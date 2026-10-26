/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#include <osgEarthPrestige/GradingLayer>
#include "Grading.h"
#include <osgEarth/ElevationPool>
#include <osgEarth/HeightFieldUtils>
#include <osgEarth/Map>
#include <osgEarth/Session>
#include <algorithm>
#include <functional>
#include <map>
#include <sstream>
#include <tuple>

using namespace osgEarth;

namespace osgEarthPrestige
{
REGISTER_OSGEARTH_LAYER_FACTORY("prestige:grading", GradingLayer);

namespace
{
    using namespace Grading;
    constexpr std::size_t maxSamples = 1000000u;
    constexpr std::size_t maxLatticeNodes = 4000000u;
    // Detected junctions: neighbor gap between seed footprints, and how far around a junction its plane is fitted.
    constexpr double junctionGap = 8.0;
    constexpr double junctionFitReach = 15.0;

    //! Checks cancellation at preparation, sampling, and rasterization boundaries.
    bool canceled(ProgressCallback* progress) { return progress && progress->isCanceled(); }

    //! Tests finite coordinates before passing user geometry to transforms or indexes.
    bool finite(const osg::Vec3d& p)
    {
        return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
    }

    //! Compares exact geometry, including components and holes; absolute profiles also compare source Z.
    //! A shared feature ID alone cannot identify duplicates when a tiled source returns clipped fragments.
    bool sameGeometry(const Geometry& a, const Geometry& b, bool includeZ)
    {
        if (a.getType() != b.getType() || a.size() != b.size()) return false;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (a[i].x() != b[i].x() || a[i].y() != b[i].y() || (includeZ && a[i].z() != b[i].z())) return false;
        if (a.getType() == Geometry::TYPE_MULTI)
        {
            const auto& ac = static_cast<const MultiGeometry&>(a).getComponents();
            const auto& bc = static_cast<const MultiGeometry&>(b).getComponents();
            if (ac.size() != bc.size()) return false;
            for (std::size_t i = 0; i < ac.size(); ++i)
                if (!sameGeometry(*ac[i], *bc[i], includeZ)) return false;
        }
        else if (a.getType() == Geometry::TYPE_POLYGON)
        {
            const auto& ah = static_cast<const Polygon&>(a).getHoles();
            const auto& bh = static_cast<const Polygon&>(b).getHoles();
            if (ah.size() != bh.size()) return false;
            for (std::size_t i = 0; i < ah.size(); ++i)
                if (!sameGeometry(*ah[i], *bh[i], includeZ)) return false;
        }
        return true;
    }

    //! Tests conservative metric bounds against the tile's ECEF sample bounds.
    bool intersects(const osg::Vec3d& a, const osg::Vec3d& b, double radius,
        const osg::Vec3d& lo, const osg::Vec3d& hi)
    {
        for (unsigned axis = 0; axis < 3; ++axis)
            if (std::max(a[axis], b[axis]) + radius < lo[axis] ||
                std::min(a[axis], b[axis]) - radius > hi[axis]) return false;
        return true;
    }

    //! Builds a ground-meter halo in geographic space before selecting source tiles or bounds.
    //! This avoids projected scale distortion and GeoExtent::expand's total-width (not per-side) convention.
    bool featureQueries(const TileKey& key, const FeatureProfile* profile, double radius, std::vector<Query>& queries)
    {
        if (!profile || !profile->getSRS()) return false;
        const auto* geo = SpatialReference::get("wgs84");
        GeoExtent extent = key.getExtent().transform(geo);
        if (!extent.isValid()) return false;
        const double dy = osg::RadiansToDegrees(radius / 6300000.0);
        const double latitude = std::max(std::abs(extent.yMin()), std::abs(extent.yMax())) + dy;
        const double dx = latitude >= 89.999 ? 180.0 : dy/std::cos(osg::DegreesToRadians(latitude));
        extent.expand(2.0*dx, 2.0*dy);
        if (profile->isTiled())
        {
            const auto* tiling = profile->getTilingProfile();
            unsigned lod = tiling->getEquivalentLOD(key.getProfile(), key.getLOD());
            if (profile->getMaxLevel() >= 0) lod = std::min(lod, static_cast<unsigned>(profile->getMaxLevel()));
            if (lod < static_cast<unsigned>(profile->getFirstLevel())) return true;
            std::vector<TileKey> keys;
            tiling->getIntersectingTiles(extent, lod, keys);
            for (const auto& tile : keys) queries.emplace_back(tile);
        }
        else
        {
            GeoExtent west, east;
            std::vector<GeoExtent> extents;
            if (extent.splitAcrossAntimeridian(west, east)) { extents.push_back(west); extents.push_back(east); }
            else extents.push_back(extent);
            for (const auto& part : extents)
            {
                const auto transformed = part.transform(profile->getSRS());
                if (!transformed.isValid()) return false;
                Query query;
                query.bounds() = transformed.bounds();
                queries.push_back(query);
            }
        }
        return true;
    }

    //! Builds exact grid locations in WGS84 and ECEF. Geographic/Mercator axes are separable, so only
    //! one diagonal needs projection and one sine/cosine pair per column; other projections use a full transform.
    bool makeGrid(const GeoExtent& extent, unsigned size, std::vector<osg::Vec3d>& geo,
        std::vector<osg::Vec3d>& world, osg::Vec3d& lo, osg::Vec3d& hi)
    {
        const auto* geographic = SpatialReference::get("wgs84");
        const auto* srs = extent.getSRS();
        const auto& ellipsoid = geographic->getEllipsoid();
        const bool separable = srs->isHorizEquivalentTo(geographic) ||
            srs->isHorizEquivalentTo(SpatialReference::get("spherical-mercator")) ||
            (srs->isMercator() && srs->getGeodeticSRS()->isHorizEquivalentTo(geographic));
        geo.resize(size*size);
        world.resize(size*size);
        if (separable)
        {
            std::vector<osg::Vec3d> diagonal;
            std::vector<double> cosine(size), sine(size);
            for (unsigned i = 0; i < size; ++i)
                diagonal.emplace_back(extent.xMin()+extent.width()*i/(size-1),
                    extent.yMin()+extent.height()*i/(size-1), 0.0);
            if (!srs->transform(diagonal, geographic)) return false;
            for (unsigned i = 0; i < size; ++i)
            {
                cosine[i] = std::cos(osg::DegreesToRadians(diagonal[i].x()));
                sine[i] = std::sin(osg::DegreesToRadians(diagonal[i].x()));
            }
            for (unsigned row = 0; row < size; ++row)
            {
                const auto meridian = ellipsoid.geodeticToGeocentric(osg::Vec3d(0.0, diagonal[row].y(), 0.0));
                for (unsigned col = 0; col < size; ++col)
                {
                    const unsigned index = row*size+col;
                    geo[index].set(diagonal[col].x(), diagonal[row].y(), 0.0);
                    world[index].set(meridian.x()*cosine[col], meridian.x()*sine[col], meridian.z());
                }
            }
        }
        else
        {
            for (unsigned row = 0; row < size; ++row)
            for (unsigned col = 0; col < size; ++col)
                geo[row*size+col].set(extent.xMin()+extent.width()*col/(size-1),
                    extent.yMin()+extent.height()*row/(size-1), 0.0);
            if (!srs->transform(geo, geographic)) return false;
            for (std::size_t i = 0; i < geo.size(); ++i)
            {
                geo[i].z() = 0.0;
                world[i] = ellipsoid.geodeticToGeocentric(geo[i]);
            }
        }
        for (const auto& point : world)
        {
            if (!finite(point)) return false;
            for (unsigned axis = 0; axis < 3; ++axis)
            {
                lo[axis] = std::min(lo[axis], point[axis]);
                hi[axis] = std::max(hi[axis], point[axis]);
            }
        }
        return true;
    }

    //! Finds a point strictly inside a local XY polygon (holes excluded): the bounds center when inside,
    //! otherwise a scan between vertex ordinates, which avoids vertex/ray degeneracy. Returns false when none
    //! is found or the job is canceled.
    bool interiorPoint(const Polygon& polygon, osg::Vec3d& anchor, ProgressCallback* progress)
    {
        anchor = polygon.getBounds().center();
        anchor.z() = 0.0;
        if (polygon.contains2D(anchor.x(), anchor.y())) return true;
        std::vector<double> ys;
        ConstGeometryIterator allRings(&polygon, true);
        while (allRings.hasMore()) for (const auto& p : *allRings.next()) ys.push_back(p.y());
        std::sort(ys.begin(), ys.end());
        ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
        for (std::size_t i = 1; i < ys.size(); ++i)
        {
            if (canceled(progress)) return false;
            const double y = 0.5*(ys[i-1]+ys[i]);
            std::vector<double> xs;
            ConstGeometryIterator edgesByRing(&polygon, true);
            while (edgesByRing.hasMore())
            {
                ConstSegmentIterator edges(edgesByRing.next(), true);
                while (edges.hasMore())
                {
                    const auto& edge = edges.next();
                    const auto& a = edge.first; const auto& b = edge.second;
                    if ((a.y() > y) != (b.y() > y))
                        xs.push_back(a.x()+(b.x()-a.x())*(y-a.y())/(b.y()-a.y()));
                }
            }
            std::sort(xs.begin(), xs.end());
            for (std::size_t j = 1; j < xs.size(); ++j)
            {
                anchor.set(0.5*(xs[j-1]+xs[j]), y, 0.0);
                if (polygon.contains2D(anchor.x(), anchor.y())) return true;
            }
        }
        return false;
    }

    //! Returns the unsigned XY area of local rings: the first ring minus the others (holes).
    double ringArea(const std::vector<std::vector<osg::Vec3d>>& rings)
    {
        double area = 0.0;
        for (std::size_t r = 0; r < rings.size(); ++r)
        {
            double twice = 0.0;
            const auto& ring = rings[r];
            for (std::size_t i = 0; i < ring.size(); ++i)
            {
                const auto& a = ring[i];
                const auto& b = ring[(i+1) % ring.size()];
                twice += a.x()*b.y()-b.x()*a.y();
            }
            area += (r == 0 ? 0.5 : -0.5)*std::abs(twice);
        }
        return std::max(0.0, area);
    }

    //! Owns one paging job's inputs and scratch; no mutable state is shared between workers.
    struct Compiler
    {
        const GradingLayer::Options& options;
        const SpatialReference* geographic = SpatialReference::get("wgs84");
        ElevationPool* pool;
        ElevationPool::WorkingSet workingSet;
        ProgressCallback* progress;
        osg::Vec3d lo, hi;
        std::size_t sampleCount = 0;
        std::shared_ptr<Junctions> junctions = std::make_shared<Junctions>();
        Surface surface;
        std::map<std::string, unsigned> profileGroups;
        struct Anchor { osg::Vec3d world; double height, ramp; };
        std::vector<NetworkLine> lines;
        // Destination surfaces per line; null is the main surface. Shared ownership keeps paired profiles alive
        // until emission, even when their footprint missed this tile.
        std::vector<std::vector<std::shared_ptr<Surface>>> destinations;
        std::map<std::tuple<unsigned, std::int64_t, std::int64_t, std::int64_t, std::int64_t, std::size_t>, unsigned> lineIndex;
        std::vector<Anchor> anchors;
        // Water bodies share one union in a tangent frame at the tile center.
        osg::Vec3d center;
        osg::Matrixd localToWorld, worldToLocal;
        std::shared_ptr<WaterShape> waterShape = std::make_shared<WaterShape>();
        struct WaterBody { osg::Vec3d anchor; double depth; };
        std::vector<WaterBody> smallBodies;

        //! Retains the caller-owned pool for this stack lifetime and selects only preceding elevation layers.
        //! The shared surface carries detected junction planes and edge-relative weighting (level sections).
        Compiler(const GradingLayer::Options& o, ElevationPool* p, const ElevationLayerVector& layers,
            ProgressCallback* callback, const osg::Vec3d& minimum, const osg::Vec3d& maximum) :
            options(o), pool(p), progress(callback), lo(minimum), hi(maximum), surface(junctions, false, true)
        {
            workingSet.setElevationLayers(layers);
            surface.setWaterBanks(options.policy().get() == "water");
            auto geo = geographic->getEllipsoid().geocentricToGeodetic((lo+hi)*0.5);
            geo.z() = 0.0;
            center = geographic->getEllipsoid().geodeticToGeocentric(geo);
            localToWorld = geographic->getEllipsoid().geocentricToLocalToWorld(center);
            worldToLocal = osg::Matrixd::inverse(localToWorld);
        }

        //! Samples a compact batch in the map SRS at a fixed physical resolution, independent of output LOD.
        bool sample(std::vector<osg::Vec3d>& positions)
        {
            if (positions.empty()) return !canceled(progress);
            if (canceled(progress) || positions.size() > maxSamples - sampleCount) return false;
            sampleCount += positions.size();
            if (!geographic->transform(positions, pool->getMapSRS())) return false;
            const int count = pool->sampleMapCoords(positions.begin(), positions.end(),
                Distance(options.elevationResolution().get(), Units::METERS), &workingSet, progress);
            return count >= 0 && !canceled(progress);
        }

        //! Copies source XY onto WGS84's zero-height ellipsoid; never changes cursor-owned features or uses source Z.
        bool position(const osg::Vec3d& input, const SpatialReference* srs, osg::Vec3d& world) const
        {
            osg::Vec3d xy(input.x(), input.y(), 0.0), geo;
            if (!finite(xy) || !srs->transform(xy, geographic, geo) || !finite(geo) || std::abs(geo.y()) > 90.0)
                return false;
            geo.z() = 0.0;
            return geographic->transformToWorld(geo, world) && finite(world);
        }

        //! Returns lon/lat on the reference ellipsoid from a chord position, discarding chord sag.
        osg::Vec3d onSurface(const osg::Vec3d& world) const
        {
            auto geo = geographic->getEllipsoid().geocentricToGeodetic(world);
            geo.z() = 0.0;
            return geo;
        }

        //! Preserves a solved piecewise-linear profile in feature-datum heights, converted to ellipsoid meters.
        //! Rejects invalid/vertical segments; keeps complete source stations and never resamples or smooths Z.
        bool absoluteRoad(const Geometry& geometry, const SpatialReference* srs, double radius, double buffer,
            unsigned group, Surface* destination, double offset)
        {
            if (!destination || geometry.size() < 2 || geometry.size() > maxSamples-sampleCount) return false;
            sampleCount += geometry.size();
            std::vector<osg::Vec3d> world, geographicPoints;
            std::vector<double> stations;
            world.reserve(geometry.size()); geographicPoints.reserve(geometry.size());
            for (const auto& input : geometry)
            {
                if (canceled(progress) || !finite(input) || !validHeight(input.z())) return false;
                osg::Vec3d p, geo;
                if (!position(input, srs, p) || !srs->transform(input, geographic, geo) || !finite(geo) ||
                    !validHeight(geo.z())) return false;
                stations.push_back(world.empty() ? 0.0 : stations.back()+(p-world.back()).length());
                geo.z() += offset;
                world.push_back(p); geographicPoints.push_back(geo);
            }
            for (std::size_t i = 1; i < world.size(); ++i)
            {
                const double length = (world[i]-world[i-1]).length();
                if (length > 100000.0) return false;
                if (length < 1e-6)
                {
                    if (std::abs(geographicPoints[i].z()-geographicPoints[i-1].z()) > 1e-6) return false;
                    continue;
                }
                const double sag = length*length / (4.0*geographic->getEllipsoid().getRadiusPolar());
                if (!intersects(world[i-1], world[i], radius+buffer+sag, lo, hi)) continue;
                Segment segment;
                segment.a = world[i-1]; segment.b = world[i];
                segment.ha = geographicPoints[i-1].z(); segment.hb = geographicPoints[i].z();
                segment.radius = radius; segment.buffer = buffer; segment.group = group;
                segment.stationA = stations[i-1]; segment.stationB = stations[i];
                segment.profileLength = stations.back(); segment.transition = options.profileTransition().get();
                destination->add(segment);
            }
            return true;
        }

        //! Registers one source polyline for the network pass. Identical lines that several features share
        //! (same group and geometry) are profiled once and emitted to every destination surface.
        bool addLine(const Geometry& geometry, const SpatialReference* srs, double radius, double buffer, unsigned group,
            const std::shared_ptr<Surface>& destination, double depth = 0.0)
        {
            if (geometry.size() < 2 || radius <= 0.0) return true;
            NetworkLine line;
            line.radius = radius; line.buffer = buffer; line.group = group; line.depth = depth;
            for (const auto& input : geometry)
            {
                if (canceled(progress)) return false;
                osg::Vec3d world;
                if (!position(input, srs, world)) return false;
                if (!line.world.empty() && (world-line.world.back()).length() < 1e-6) continue;
                if (!line.world.empty() && (world-line.world.back()).length() > 100000.0) return false;
                line.world.push_back(world);
            }
            if (line.world.size() < 2) return true;
            const auto key = std::make_tuple(group, vertexKey(line.world.front()).x, vertexKey(line.world.front()).y,
                vertexKey(line.world.back()).x, vertexKey(line.world.back()).y, line.world.size());
            const auto found = lineIndex.find(key);
            if (found != lineIndex.end())
            {
                auto& targets = destinations[found->second];
                if (std::find(targets.begin(), targets.end(), destination) == targets.end()) targets.push_back(destination);
                return true;
            }
            lineIndex.emplace(key, static_cast<unsigned>(lines.size()));
            lines.push_back(line);
            destinations.push_back(std::vector<std::shared_ptr<Surface>>(1, destination));
            return true;
        }

        //! Registers an anchor: road profiles touching it must reach ``height`` (ellipsoid meters) over ``ramp``.
        void addAnchor(const osg::Vec3d& world, double height, double ramp)
        {
            anchors.push_back(Anchor{world, height, ramp});
        }

        struct Sample { osg::Vec3d world, geo; double station; unsigned line; };
        std::vector<std::vector<PathPiece>> paths;
        std::vector<std::vector<Sample>> samples;
        // Profile heights per path sample; for water, the channel surface before bed depth.
        std::vector<std::vector<double>> heights;

        //! Splits path ``p`` into runs of adjacent valid samples and visits each; never bridges a data gap or a
        //! skipped span. Returns false as soon as ``visit`` does.
        bool runs(std::size_t p, const std::function<bool(std::size_t, std::size_t)>& visit) const
        {
            const double spacing = options.sampleSpacing().get();
            std::size_t first = 0;
            const auto& path = samples[p];
            while (first < path.size())
            {
                if (!validHeight(heights[p][first])) { ++first; continue; }
                std::size_t end = first+1;
                while (end < path.size() && validHeight(heights[p][end]) &&
                    path[end].station-path[end-1].station <= spacing*1.01) ++end;
                if (!visit(first, end)) return false;
                first = end;
            }
            return true;
        }

        //! Builds profiles for every registered line. Lines are chained through shared ends so smoothing and
        //! junction transitions are continuous; vertices are sampled per source segment, so sample positions do
        //! not depend on which neighboring lines a paging job happened to load. Water paths become the channel
        //! surface (a windowed upstream minimum, smoothed); road paths are pinned to anchors and blended into
        //! fitted junction planes (one per detected junction, fitted over its close neighbors). Every step has
        //! bounded reach, which the caller's query halo covers, so adjacent paging jobs agree.
        bool profileNetwork(double halo, bool water)
        {
            if (lines.empty()) return !canceled(progress);
            const double spacing = options.sampleSpacing().get();
            paths = chainLines(lines, water);
            samples.assign(paths.size(), std::vector<Sample>());
            std::vector<osg::Vec3d> positions;
            for (std::size_t p = 0; p < paths.size(); ++p)
            {
                double base = 0.0;
                bool first = true;
                for (const auto& piece : paths[p])
                {
                    const auto& line = lines[piece.line];
                    const std::size_t n = line.world.size();
                    for (std::size_t k = 1; k < n; ++k)
                    {
                        if (canceled(progress)) return false;
                        const auto& a = line.world[piece.reversed ? n-k : k-1];
                        const auto& b = line.world[piece.reversed ? n-k-1 : k];
                        const double length = (b-a).length();
                        const double sag = length*length / (4.0*geographic->getEllipsoid().getRadiusPolar());
                        if (intersects(a, b, halo+sag, lo, hi))
                        {
                            const unsigned count = std::max(1u, static_cast<unsigned>(std::ceil(length/spacing)));
                            for (unsigned j = first ? 0u : 1u; j <= count; ++j)
                            {
                                const double t = double(j)/count;
                                const auto geo = onSurface(a + (b-a)*t);
                                const auto world = geographic->getEllipsoid().geodeticToGeocentric(geo);
                                if (!intersects(world, world, halo, lo, hi)) continue;
                                if (positions.size() >= maxSamples || positions.size() >= maxSamples-sampleCount)
                                    return false;
                                samples[p].push_back(Sample{world, geo, base+length*t, piece.line});
                                positions.push_back(geo);
                            }
                        }
                        base += length;
                        first = false;
                    }
                }
            }
            if (!sample(positions)) return false;
            std::size_t cursor = 0;
            heights.assign(paths.size(), std::vector<double>());
            for (std::size_t p = 0; p < paths.size(); ++p)
                for (std::size_t i = 0; i < samples[p].size(); ++i)
                    heights[p].push_back(positions[cursor++].z());

            const double smoothing = options.smoothingDistance().get();
            const unsigned passes = std::max(1u, options.smoothingPasses().get());
            for (std::size_t p = 0; p < paths.size(); ++p)
            {
                std::vector<double> station;
                for (const auto& s : samples[p]) station.push_back(s.station);
                const bool ok = runs(p, [&](std::size_t first, std::size_t end)
                {
                    if (canceled(progress)) return false;
                    const double window = options.flowWindow().get();
                    std::vector<double> terrain;
                    if (water)
                    {
                        terrain = heights[p];
                        smoothProfile(station, terrain, first, end, smoothing, passes);
                        flowMinimum(station, heights[p], first, end, window);
                    }
                    smoothProfile(station, heights[p], first, end, smoothing, passes);
                    if (water)
                    {
                        flowMinimum(station, heights[p], first, end, window);
                        // Where a line strays from the terrain's valley (coarse elevation, misregistration), the
                        // upstream minimum would hold a low bed through a hillside; never cut deeper than the limit
                        // below the smoothed terrain, then round the clamped corners once more.
                        const double limit = options.channelCutLimit().get();
                        for (std::size_t i = first; i < end; ++i)
                            heights[p][i] = std::max(heights[p][i], terrain[i]-limit);
                        smoothProfile(station, heights[p], first, end, smoothing, 1u);
                        return true;
                    }
                    // Anchors pin the profile with a quintic ramp centered on the touching sample.
                    for (const auto& anchor : anchors)
                    {
                        std::size_t at = end;
                        for (std::size_t i = first; i < end; ++i)
                            if ((samples[p][i].world-anchor.world).length() <= 1.0) { at = i; break; }
                        if (at == end) continue;
                        const double delta = anchor.height-heights[p][at];
                        for (std::size_t i = first; i < end; ++i)
                        {
                            const double d = std::abs(station[i]-station[at]);
                            if (d < anchor.ramp) heights[p][i] += delta*(1.0-ease(d/anchor.ramp));
                        }
                    }
                    return true;
                });
                if (!ok) return false;
            }

            // Paired footprints carry the compiler's own junctions, and their line widths are profile search
            // widths rather than pavement, so detection applies to plain line features only.
            if (!water && options.junctionDetection().get() && !options.roadFootprints().get())
            {
                const double margin = options.junctionMargin().get();
                const double transition = options.junctionTransition().get();
                const double seedRadius = 0.5*options.maxLineWidth().get() + margin;
                // A junction's footprint spans its neighborhood: up to three seed radii plus a gap from its seed.
                const double footprintReach = 3.0*seedRadius + junctionGap;
                std::vector<JunctionSeed> seeds;
                for (const auto& seed : findJunctions(lines, margin))
                    if (intersects(seed.world, seed.world, footprintReach + transition + 2.0*seedRadius + junctionGap,
                        lo, hi)) seeds.push_back(seed);
                for (unsigned a = 0; a < seeds.size(); ++a)
                {
                    if (canceled(progress)) return false;
                    if (!intersects(seeds[a].world, seeds[a].world, footprintReach + transition, lo, hi)) continue;
                    const auto hood = neighborhood(seeds, a, junctionGap);
                    Junction j;
                    j.worldToLocal = osg::Matrixd::inverse(
                        geographic->getEllipsoid().geocentricToLocalToWorld(seeds[a].world));
                    j.transition = transition;
                    osg::Vec3d minimum(DBL_MAX, DBL_MAX, DBL_MAX), maximum(-DBL_MAX, -DBL_MAX, -DBL_MAX);
                    for (unsigned m : hood)
                    {
                        auto local = seeds[m].world*j.worldToLocal;
                        local.z() = 0.0;
                        j.capsules.push_back(Capsule{local, local, seeds[m].radius});
                        j.members.insert(j.members.end(), seeds[m].members.begin(), seeds[m].members.end());
                        for (unsigned n : hood)
                        {
                            if (n <= m || !neighbors(seeds[m], seeds[n], junctionGap)) continue;
                            auto other = seeds[n].world*j.worldToLocal;
                            other.z() = 0.0;
                            j.capsules.push_back(Capsule{local, other, std::min(seeds[m].radius, seeds[n].radius)});
                        }
                        for (unsigned axis = 0; axis < 3; ++axis)
                        {
                            minimum[axis] = std::min(minimum[axis], seeds[m].world[axis]-seeds[m].radius);
                            maximum[axis] = std::max(maximum[axis], seeds[m].world[axis]+seeds[m].radius);
                        }
                    }
                    std::sort(j.members.begin(), j.members.end());
                    j.members.erase(std::unique(j.members.begin(), j.members.end()), j.members.end());
                    // Fit the plane to the members' smoothed profiles near the junction, so arms blend in gently.
                    std::vector<osg::Vec3d> fit;
                    std::vector<double> weights;
                    for (std::size_t p = 0; p < paths.size(); ++p)
                    {
                        for (std::size_t i = 0; i < samples[p].size(); ++i)
                        {
                            const auto& s = samples[p][i];
                            if (!validHeight(heights[p][i]) ||
                                !std::binary_search(j.members.begin(), j.members.end(), lines[s.line].group)) continue;
                            auto local = s.world*j.worldToLocal;
                            local.z() = 0.0;
                            if (j.distance2(local) > junctionFitReach*junctionFitReach) continue;
                            fit.emplace_back(local.x(), local.y(), heights[p][i]);
                            weights.push_back(spacing);
                        }
                    }
                    if (!fitPlane(fit, weights, j.plane)) continue;
                    const double pad = transition + (maximum-minimum).length();
                    minimum -= osg::Vec3d(pad, pad, pad);
                    maximum += osg::Vec3d(pad, pad, pad);
                    junctions->add(j, minimum, maximum);
                }
            }

            return !canceled(progress);
        }

        //! Emits profile segments to their destination surfaces. Water channels are lowered by their bed depth and
        //! skipped inside water bodies, whose own surface governs there.
        bool emitNetwork(bool water)
        {
            std::vector<unsigned> hits;
            for (std::size_t p = 0; p < paths.size(); ++p)
            {
                const auto& path = samples[p];
                const bool ok = runs(p, [&](std::size_t first, std::size_t end)
                {
                    for (std::size_t i = first+1; i < end; ++i)
                    {
                        const auto& line = lines[path[i].line];
                        if (!intersects(path[i-1].world, path[i].world, line.radius+line.buffer, lo, hi)) continue;
                        if (water && !waterShape->empty())
                        {
                            waterShape->containing((path[i-1].world+path[i].world)*0.5*worldToLocal, hits);
                            if (!hits.empty()) continue;
                        }
                        Segment s;
                        s.a = path[i-1].world; s.b = path[i].world;
                        s.ha = heights[p][i-1] - (water ? lines[path[i-1].line].depth : 0.0);
                        s.hb = heights[p][i] - (water ? lines[path[i].line].depth : 0.0);
                        s.radius = line.radius; s.buffer = line.buffer; s.group = line.group;
                        s.extendStart = i > first+1;
                        s.extendEnd = i+1 < end;
                        for (const auto& destination : destinations[path[i].line])
                            (destination ? *destination : surface).add(s);
                    }
                    return true;
                });
                if (!ok) return false;
            }
            return !canceled(progress);
        }

        //! Adds a water body to the job's water union in the tile frame, with its outward ``bank`` and bed
        //! ``depth``; ``flowing`` bodies (river areas) follow their channel profile instead of a level. Levels
        //! come later from buildWater(). Any size is accepted; bodies reaching more than 3000 km from the tile
        //! are ignored, since the tangent frame cannot represent them.
        bool waterPolygon(const Polygon& input, const SpatialReference* srs, double bank, double depth, bool flowing)
        {
            std::vector<std::vector<osg::Vec3d>> rings;
            ConstGeometryIterator parts(&input, true);
            while (parts.hasMore())
            {
                std::vector<osg::Vec3d> ring;
                for (const auto& p : *parts.next())
                {
                    if (canceled(progress)) return false;
                    osg::Vec3d world;
                    if (!position(p, srs, world)) return false;
                    if ((world-center).length() > 3.0e6) return true;
                    auto local = world*worldToLocal;
                    local.z() = 0.0;
                    ring.push_back(local);
                }
                if (ring.size() >= 3) rings.push_back(ring);
                else if (rings.empty()) return true;
            }
            if (rings.empty()) return true;
            waterShape->add(rings, bank, depth, flowing);
            // A body that may hold no lattice node still contributes one interior sample. The rule depends on
            // the polygon alone, so every paging job that loads it agrees.
            const double lattice = std::max(options.sampleSpacing().get(), options.waterLevelRadius().get()/20.0);
            if (ringArea(rings) < 4.0*lattice*lattice)
            {
                osg::ref_ptr<Polygon> local = new Polygon();
                local->asVector() = rings.front();
                for (std::size_t r = 1; r < rings.size(); ++r)
                {
                    osg::ref_ptr<Ring> hole = new Ring();
                    hole->asVector() = rings[r];
                    local->getHoles().push_back(hole);
                }
                osg::Vec3d anchor;
                if (interiorPoint(*local, anchor, progress)) smallBodies.push_back(WaterBody{anchor, depth});
                else if (canceled(progress)) return false;
            }
            return true;
        }

        //! Builds the water level field and one pad covering the water union. Lattice rows are fixed latitudes
        //! with fixed longitudes per row, so every paging job samples identical nodes; a job samples wet nodes
        //! within the level radius of its tile and computes levels only for nodes near its tile and the water.
        //! Wet nodes of flowing bodies take the nearest channel surface (from profileNetwork()) minus their
        //! depth, so river areas fall with their river and stay level across. Antimeridian-crossing tiles use a
        //! lattice that is not continuous across 180 degrees.
        bool buildWater(const GeoExtent& tile, unsigned group)
        {
            if (waterShape->empty()) return !canceled(progress);
            const double radius = options.waterLevelRadius().get();
            // Channel surface segments in the tile frame, (x, y, height) at both ends.
            std::vector<std::pair<osg::Vec3d, osg::Vec3d>> channel;
            RTree<unsigned, double, 2> channelIndex;
            for (std::size_t p = 0; p < paths.size(); ++p)
            {
                runs(p, [&](std::size_t first, std::size_t end)
                {
                    for (std::size_t i = first+1; i < end; ++i)
                    {
                        auto a = samples[p][i-1].world*worldToLocal, b = samples[p][i].world*worldToLocal;
                        a.z() = heights[p][i-1];
                        b.z() = heights[p][i];
                        double lo2[2] = {std::min(a.x(), b.x()), std::min(a.y(), b.y())};
                        double hi2[2] = {std::max(a.x(), b.x()), std::max(a.y(), b.y())};
                        channelIndex.Insert(lo2, hi2, static_cast<unsigned>(channel.size()));
                        channel.emplace_back(a, b);
                    }
                    return true;
                });
            }
            //! Channel surface height at the nearest channel point within the level radius, or NO_DATA.
            auto channelHeight = [&](const osg::Vec3d& p)
            {
                double best = radius*radius, height = NO_DATA_VALUE;
                const double lo2[2] = {p.x()-radius, p.y()-radius}, hi2[2] = {p.x()+radius, p.y()+radius};
                channelIndex.Search(lo2, hi2, [&](unsigned i)
                {
                    const auto& a = channel[i].first;
                    const auto& b = channel[i].second;
                    const osg::Vec3d ab(b.x()-a.x(), b.y()-a.y(), 0.0), ap(p.x()-a.x(), p.y()-a.y(), 0.0);
                    const double t = ab.length2() > 0.0 ? unit((ap*ab)/ab.length2()) : 0.0;
                    const double d2 = (ap-ab*t).length2();
                    if (d2 < best) { best = d2; height = a.z()+(b.z()-a.z())*t; }
                    return RTREE_KEEP_SEARCHING;
                });
                return height;
            };
            const double spacing = std::max(options.sampleSpacing().get(), radius/20.0);
            const double bankReach = waterShape->maxBank() + 2.0*spacing;
            const double margin = std::max(radius, bankReach) + 3.0*spacing;
            const auto& ellipsoid = geographic->getEllipsoid();
            const double dy = osg::RadiansToDegrees(margin/6300000.0);
            const double latMin = std::max(-89.9, tile.yMin()-dy), latMax = std::min(89.9, tile.yMax()+dy);
            const double dlat = spacing/111320.0;
            std::vector<osg::Vec3d> wetLocal, wetGeo, leveled, flowingWet;
            std::vector<double> wetDepth;
            std::vector<unsigned> hits;
            std::size_t nodes = 0;
            for (long long row = static_cast<long long>(std::ceil(latMin/dlat));
                row <= static_cast<long long>(std::floor(latMax/dlat)); ++row)
            {
                if (canceled(progress)) return false;
                const double lat = row*dlat;
                const double dlon = dlat/std::max(0.05, std::cos(osg::DegreesToRadians(lat)));
                const double reachLat = std::min(89.9, std::max(std::abs(tile.yMin()), std::abs(tile.yMax()))+dy);
                const double dx = dy/std::max(0.05, std::cos(osg::DegreesToRadians(reachLat)));
                for (long long col = static_cast<long long>(std::ceil((tile.xMin()-dx)/dlon));
                    col <= static_cast<long long>(std::floor((tile.xMax()+dx)/dlon)); ++col)
                {
                    const osg::Vec3d geo(col*dlon, lat, 0.0);
                    const auto world = ellipsoid.geodeticToGeocentric(geo);
                    if (!intersects(world, world, margin, lo, hi)) continue;
                    if (++nodes > maxLatticeNodes) return false;
                    auto local = world*worldToLocal;
                    local.z() = 0.0;
                    waterShape->containing(local, hits);
                    if (!hits.empty())
                    {
                        const double channelSurface = waterShape->flowing(hits) ? channelHeight(local) : NO_DATA_VALUE;
                        if (validHeight(channelSurface))
                            flowingWet.emplace_back(local.x(), local.y(), channelSurface-waterShape->depth(hits));
                        else
                        {
                            wetLocal.push_back(local);
                            wetGeo.push_back(geo);
                            wetDepth.push_back(waterShape->depth(hits));
                        }
                    }
                    if (intersects(world, world, 3.0*spacing, lo, hi) &&
                        (!hits.empty() || waterShape->near(local, bankReach)))
                        leveled.push_back(local);
                }
            }
            for (const auto& body : smallBodies)
            {
                const auto geo = onSurface(body.anchor*localToWorld);
                const auto world = ellipsoid.geodeticToGeocentric(geo);
                if (!intersects(world, world, margin, lo, hi)) continue;
                wetLocal.push_back(body.anchor);
                wetGeo.push_back(geo);
                wetDepth.push_back(body.depth);
            }
            if (leveled.empty() || (wetGeo.empty() && flowingWet.empty())) return !canceled(progress);
            if (!sample(wetGeo)) return false;
            std::vector<osg::Vec3d> wet;
            for (std::size_t i = 0; i < wetGeo.size(); ++i)
                if (validHeight(wetGeo[i].z()))
                    wet.emplace_back(wetLocal[i].x(), wetLocal[i].y(), wetGeo[i].z()-wetDepth[i]);
            if (wet.empty() && flowingWet.empty()) return true;
            Pad pad;
            pad.shape = waterShape;
            pad.worldToLocal = worldToLocal;
            pad.group = group;
            pad.level = std::make_shared<WaterField>(wet, flowingWet, leveled, radius, spacing, 0.2);
            // Pixels exist only inside the tile; one spacing of slack covers chord sag of the bounds.
            surface.add(pad, lo-osg::Vec3d(spacing, spacing, spacing), hi+osg::Vec3d(spacing, spacing, spacing));
            return !canceled(progress);
        }

        //! Fits one common junction plane from its own centerline terrain samples, independently of output tiles.
        //! Only explicit member profiles acquire the transition; disconnected nearby roads remain unaffected.
        bool junction(const std::vector<const Geometry*>& lines, const std::vector<const Polygon*>& polygons,
            const SpatialReference* srs, const std::vector<unsigned>& members)
        {
            if (polygons.empty() || polygons.front()->empty()) return false;
            osg::Vec3d origin;
            if (!position(polygons.front()->front(), srs, origin)) return false;
            Junction j;
            j.worldToLocal = osg::Matrixd::inverse(geographic->getEllipsoid().geocentricToLocalToWorld(origin));
            j.members = members;
            std::sort(j.members.begin(), j.members.end());
            j.members.erase(std::unique(j.members.begin(), j.members.end()), j.members.end());
            j.transition = options.junctionTransition().get();
            osg::Vec3d minimum(DBL_MAX, DBL_MAX, DBL_MAX), maximum(-DBL_MAX, -DBL_MAX, -DBL_MAX);
            for (const auto* input : polygons)
            {
                osg::ref_ptr<Polygon> polygon = static_cast<Polygon*>(input->clone());
                GeometryIterator rings(polygon.get(), true);
                while (rings.hasMore()) for (auto& p : *rings.next())
                {
                    osg::Vec3d world;
                    if (!position(p, srs, world) || (world-origin).length() > 100000.0) return false;
                    for (unsigned axis = 0; axis < 3; ++axis)
                    {
                        minimum[axis] = std::min(minimum[axis], world[axis]);
                        maximum[axis] = std::max(maximum[axis], world[axis]);
                    }
                    p = world*j.worldToLocal;
                    p.z() = 0.0;
                }
                j.footprints.push_back(polygon);
            }
            const double span = (maximum-minimum).length();
            const double pad = j.transition + span*span/geographic->getEllipsoid().getRadiusPolar();
            minimum -= osg::Vec3d(pad, pad, pad);
            maximum += osg::Vec3d(pad, pad, pad);
            if (!intersects(minimum, maximum, 0.0, lo, hi)) return true;

            std::vector<osg::Vec3d> samples, local;
            std::vector<double> weights;
            // Sampling is anchored to complete source segments, never to a tile's clipped portion.
            const double spacing = std::min(5.0, options.sampleSpacing().get());
            for (const auto* line : lines)
            {
                for (std::size_t i = 1; i < line->size(); ++i)
                {
                    if (canceled(progress)) return false;
                    osg::Vec3d a, b;
                    if (!position((*line)[i-1], srs, a) || !position((*line)[i], srs, b)) return false;
                    const double length = (b-a).length();
                    if (length > 100000.0) return false;
                    if (length < 1e-6 || !intersects(a, b, 2.0, minimum, maximum)) continue;
                    const unsigned count = static_cast<unsigned>(std::ceil(length/spacing));
                    for (unsigned k = 0; k <= count; ++k)
                    {
                        const auto geo = onSurface(a+(b-a)*(double(k)/count));
                        const auto p = geographic->getEllipsoid().geodeticToGeocentric(geo)*j.worldToLocal;
                        if (j.distance2(p) > 4.0) continue;
                        if (samples.size() >= maxSamples) return false;
                        samples.push_back(geo);
                        local.push_back(p);
                        weights.push_back(length/count*(k == 0 || k == count ? 0.5 : 1.0));
                    }
                }
            }
            if (!sample(samples)) return false;
            double weight = 0.0, x = 0.0, y = 0.0, h = 0.0;
            unsigned valid = 0;
            for (std::size_t i = 0; i < samples.size(); ++i) if (validHeight(samples[i].z()))
            {
                weight += weights[i];
                x += weights[i]*local[i].x(); y += weights[i]*local[i].y(); h += weights[i]*samples[i].z();
                ++valid;
            }
            // Never fill a missing elevation region from an absent or single terrain sample.
            if (valid < 3 || weight <= 0.0) return true;
            x /= weight; y /= weight; h /= weight;
            double xx = 0.0, xy = 0.0, yy = 0.0, xh = 0.0, yh = 0.0;
            for (std::size_t i = 0; i < samples.size(); ++i) if (validHeight(samples[i].z()))
            {
                const double dx = local[i].x()-x, dy = local[i].y()-y, dh = samples[i].z()-h;
                xx += weights[i]*dx*dx; xy += weights[i]*dx*dy; yy += weights[i]*dy*dy;
                xh += weights[i]*dx*dh; yh += weights[i]*dy*dh;
            }
            // A tiny ridge term stabilizes collinear approaches without inventing an unconstrained cross-slope.
            const double regularization = std::max(1e-8, (xx+yy)*1e-8);
            xx += regularization; yy += regularization;
            const double determinant = xx*yy-xy*xy;
            if (!(determinant > 0.0)) return true;
            j.plane.x() = (xh*yy-yh*xy)/determinant;
            j.plane.y() = (yh*xx-xh*xy)/determinant;
            j.plane.z() = h-j.plane.x()*x-j.plane.y()*y;
            if (!finite(j.plane)) return false;
            junctions->add(j, minimum, maximum);
            return true;
        }

        //! Creates a round level pad for each point, sampling its center once.
        bool points(const Geometry& geometry, const SpatialReference* srs, double radius, double buffer, unsigned& group)
        {
            for (const auto& point : geometry)
            {
                osg::Vec3d world;
                if (!position(point, srs, world)) return false;
                if (!intersects(world, world, radius+buffer, lo, hi)) continue;
                std::vector<osg::Vec3d> samples(1, onSurface(world));
                if (!sample(samples)) return false;
                Segment segment;
                segment.a = segment.b = world;
                segment.ha = segment.hb = samples.front().z();
                segment.radius = radius; segment.buffer = buffer; segment.group = group++;
                surface.add(segment);
            }
            return true;
        }

        //! Masks a supplied centerline profile, or flattens to one interior sample when no profile is supplied.
        //! Edges are local tangent-plane segments; polygons spanning more than 100 km are unsupported.
        bool polygon(const Polygon& input, const SpatialReference* srs, double buffer, unsigned group,
            std::shared_ptr<const Surface> profile = nullptr, double inset = 0.0)
        {
            if (input.size() < 3) return true;
            osg::Vec3d origin;
            if (!position(input.front(), srs, origin)) return false;
            const auto localToWorld = geographic->getEllipsoid().geocentricToLocalToWorld(origin);
            const auto worldToLocal = osg::Matrixd::inverse(localToWorld);
            osg::ref_ptr<Polygon> polygon = static_cast<Polygon*>(input.clone());
            GeometryIterator rings(polygon.get(), true);
            osg::Vec3d minimum(DBL_MAX, DBL_MAX, DBL_MAX), maximum(-DBL_MAX, -DBL_MAX, -DBL_MAX);
            while (rings.hasMore())
            {
                auto* ring = rings.next();
                for (auto& p : *ring)
                {
                    osg::Vec3d world;
                    if (!position(p, srs, world) || (world-origin).length() > 100000.0) return false;
                    for (unsigned axis = 0; axis < 3; ++axis)
                    {
                        minimum[axis] = std::min(minimum[axis], world[axis]);
                        maximum[axis] = std::max(maximum[axis], world[axis]);
                    }
                    p = world * worldToLocal;
                    p.z() = 0.0;
                }
            }
            const double span = (maximum-minimum).length();
            const double padding = buffer + span*span/geographic->getEllipsoid().getRadiusPolar();
            minimum -= osg::Vec3d(padding, padding, padding);
            maximum += osg::Vec3d(padding, padding, padding);
            if (!intersects(minimum, maximum, 0.0, lo, hi)) return true;

            if (profile)
            {
                Pad pad;
                pad.polygon = polygon; pad.worldToLocal = worldToLocal;
                pad.buffer = buffer; pad.group = group; pad.profile = profile;
                pad.inset = inset;
                pad.useProfileAlpha = options.profileElevation().get() == "absolute";
                surface.add(pad, minimum, maximum);
                return true;
            }

            osg::Vec3d anchor;
            if (!interiorPoint(*polygon, anchor, progress)) return !canceled(progress);
            std::vector<osg::Vec3d> samples(1, onSurface(anchor*localToWorld));
            if (!sample(samples)) return false;
            Pad pad;
            pad.polygon = polygon; pad.worldToLocal = worldToLocal;
            pad.height = samples.front().z(); pad.buffer = buffer; pad.group = group;
            surface.add(pad, minimum, maximum);
            return true;
        }
    };
}

Config GradingLayer::Options::getConfig() const
{
    Config conf = ElevationLayer::Options::getConfig();
    featureSource().set(conf, "features");
    conf.set("line_width", lineWidth());
    conf.set("buffer_width", bufferWidth());
    conf.set("max_line_width", maxLineWidth());
    conf.set("max_buffer_width", maxBufferWidth());
    conf.set("sample_spacing", sampleSpacing());
    conf.set("smoothing_distance", smoothingDistance());
    conf.set("elevation_resolution", elevationResolution());
    conf.set("road_footprints", roadFootprints());
    conf.set("junction", junction());
    conf.set("profile_id_attribute", profileIdAttribute());
    conf.set("junction_transition", junctionTransition());
    conf.set("profile_elevation", profileElevation());
    conf.set("profile_transition", profileTransition());
    conf.set("profile_offset", profileOffset());
    conf.set("footprint_inset", footprintInset());
    conf.set("cut_only", cutOnly());
    conf.set("policy", policy());
    conf.set("smoothing_passes", smoothingPasses());
    conf.set("junction_detection", junctionDetection());
    conf.set("junction_margin", junctionMargin());
    conf.set("depth", depth());
    conf.set("anchor_height", anchorHeight());
    conf.set("anchor_ramp", anchorRamp());
    conf.set("max_anchor_ramp", maxAnchorRamp());
    conf.set("flow_window", flowWindow());
    conf.set("water_level_radius", waterLevelRadius());
    conf.set("bank_slope", bankSlope());
    conf.set("flowing", flowing());
    conf.set("channel_cut_limit", channelCutLimit());
    return conf;
}

void GradingLayer::Options::fromConfig(const Config& conf)
{
    featureSource().get(conf, "features");
    conf.get("line_width", lineWidth());
    conf.get("buffer_width", bufferWidth());
    conf.get("max_line_width", maxLineWidth());
    conf.get("max_buffer_width", maxBufferWidth());
    conf.get("sample_spacing", sampleSpacing());
    conf.get("smoothing_distance", smoothingDistance());
    conf.get("elevation_resolution", elevationResolution());
    conf.get("road_footprints", roadFootprints());
    conf.get("junction", junction());
    conf.get("profile_id_attribute", profileIdAttribute());
    conf.get("junction_transition", junctionTransition());
    conf.get("profile_elevation", profileElevation());
    conf.get("profile_transition", profileTransition());
    conf.get("profile_offset", profileOffset());
    conf.get("footprint_inset", footprintInset());
    conf.get("cut_only", cutOnly());
    conf.get("policy", policy());
    conf.get("smoothing_passes", smoothingPasses());
    conf.get("junction_detection", junctionDetection());
    conf.get("junction_margin", junctionMargin());
    conf.get("depth", depth());
    conf.get("anchor_height", anchorHeight());
    conf.get("anchor_ramp", anchorRamp());
    conf.get("max_anchor_ramp", maxAnchorRamp());
    conf.get("flow_window", flowWindow());
    conf.get("water_level_radius", waterLevelRadius());
    conf.get("bank_slope", bankSlope());
    conf.get("flowing", flowing());
    conf.get("channel_cut_limit", channelCutLimit());
}

void GradingLayer::setFeatureSource(FeatureSource* source) { options().featureSource().setLayer(source); }
FeatureSource* GradingLayer::getFeatureSource() const { return options().featureSource().getLayer(); }

void GradingLayer::init()
{
    ElevationLayer::init();
    // Output depends on preceding elevation, whose revision is not part of this layer's persistent cache key.
    options().cachePolicy() = CachePolicy::NO_CACHE;
    options().l2CacheSize() = 0u;
}

Status GradingLayer::openImplementation()
{
    const auto& o = options();
    if ((o.profileElevation().get() != "terrain" && o.profileElevation().get() != "absolute") ||
        (o.profileElevation().get() == "absolute" && !o.roadFootprints().get()))
        return Status(Status::ConfigurationError, "Absolute profile elevation requires paired road footprints");
    if (!std::isfinite(o.maxLineWidth().get()) || o.maxLineWidth().get() <= 0.0 ||
        o.maxLineWidth().get() > 10000.0 || !std::isfinite(o.maxBufferWidth().get()) ||
        o.maxBufferWidth().get() < 0.0 || o.maxBufferWidth().get() > 10000.0 ||
        !std::isfinite(o.sampleSpacing().get()) || o.sampleSpacing().get() < 1.0 || o.sampleSpacing().get() > 100.0 ||
        !std::isfinite(o.smoothingDistance().get()) || o.smoothingDistance().get() < 0.0 ||
        o.smoothingDistance().get() > 1000.0 || !std::isfinite(o.elevationResolution().get()) ||
        o.elevationResolution().get() <= 0.0 || !std::isfinite(o.junctionTransition().get()) ||
        o.junctionTransition().get() <= 0.0 || o.junctionTransition().get() > 1000.0 ||
        !std::isfinite(o.profileTransition().get()) || o.profileTransition().get() < 0.0 ||
        o.profileTransition().get() > 1000.0 ||
        (o.policy().get() != "road" && o.policy().get() != "water") ||
        (o.policy().get() == "water" && o.roadFootprints().get()) ||
        o.smoothingPasses().get() < 1u || o.smoothingPasses().get() > 8u ||
        !std::isfinite(o.junctionMargin().get()) || o.junctionMargin().get() < 0.0 || o.junctionMargin().get() > 100.0 ||
        !std::isfinite(o.maxAnchorRamp().get()) || o.maxAnchorRamp().get() < 0.0 || o.maxAnchorRamp().get() > 2000.0 ||
        !std::isfinite(o.flowWindow().get()) || o.flowWindow().get() < 0.0 || o.flowWindow().get() > 5000.0 ||
        !std::isfinite(o.waterLevelRadius().get()) || o.waterLevelRadius().get() < 10.0 ||
        o.waterLevelRadius().get() > 2000.0 ||
        !std::isfinite(o.bankSlope().get()) || o.bankSlope().get() < 0.01 || o.bankSlope().get() > 10.0 ||
        !std::isfinite(o.channelCutLimit().get()) || o.channelCutLimit().get() < 0.0 ||
        o.channelCutLimit().get() > 1000.0 ||
        getTileSize() < 2u || getTileSize() > 1024u)
        return Status(Status::ConfigurationError, "Invalid grading widths, spacing, smoothing, resolution, or tile size");
    Status status = ElevationLayer::openImplementation();
    if (status.isError()) return status;
    if (!options().featureSource().isSet())
        return Status(Status::ConfigurationError, "Grading requires a feature source");
    status = options().featureSource().open(getReadOptions());
    if (status.isError()) return status;
    if (!getProfile())
    {
        setProfile(Profile::create(Profile::GLOBAL_GEODETIC));
        // Avoid enumerating continent-scale road datasets for grids too coarse to represent road corridors.
        // Applications with a regional profile keep their own LOD convention.
        if (!options().minLevel().isSet()) options().minLevel().init(12u);
    }
    if (getProfile()->getSRS()->getVerticalDatum())
        return Status(Status::ConfigurationError, "Grading output uses ellipsoid heights; omit the output vdatum");
    addDataExtent(DataExtent(getProfile()->getExtent(), getMinLevel(), getMaxDataLevel()));
    return Status::NoError;
}

void GradingLayer::addedToMap(const Map* map)
{
    ElevationLayer::addedToMap(map);
    options().featureSource().addedToMap(map);
    if (!getFeatureSource() || !getFeatureSource()->isOpen())
        setStatus(Status(Status::ConfigurationError, "Grading feature source could not be opened or resolved"));
    _map = map;
}

void GradingLayer::removedFromMap(const Map* map)
{
    _map = nullptr;
    options().featureSource().removedFromMap(map);
    ElevationLayer::removedFromMap(map);
}

GeoHeightField GradingLayer::createHeightFieldImplementation(const TileKey& key, ProgressCallback* progress) const
{
    osg::ref_ptr<const Map> map;
    auto* source = getFeatureSource();
    if (!_map.lock(map) || !source || !source->isOpen() || canceled(progress))
        return GeoHeightField::INVALID;
    ElevationLayerVector layers, preceding;
    map->getLayers(layers);
    bool found = false;
    for (const auto& layer : layers)
    {
        if (layer.get() == this) { found = true; break; }
        if (layer->isOpen()) preceding.push_back(layer);
    }
    // An empty WorkingSet selects ALL map layers in ElevationPool, which would recursively call us.
    if (!found || preceding.empty()) return GeoHeightField::INVALID;

    const auto& extent = key.getExtent();
    const unsigned size = getTileSize();
    std::vector<osg::Vec3d> grid, world;
    osg::Vec3d lo(DBL_MAX, DBL_MAX, DBL_MAX), hi(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    std::unique_ptr<Compiler> compiler;
    std::vector<Query> queries;
    const bool water = options().policy().get() == "water";
    // Water surfaces are set exactly inside channels and bodies (filling pits so water lies level); banks are
    // slope limited both ways. The explicit cut-only option never fills anywhere.
    const bool cutOnly = options().cutOnly().get();
    const double bankSlope = options().bankSlope().get();
    // Everything that can change a pixel here must be loaded identically by every neighboring tile, so shared
    // tile edges match exactly. Profiles: segment reach plus the smoothing window (and flow windows for water).
    // Detected junctions: a footprint spans its seed neighborhood (3 seed radii + gap), whose seeds see their
    // own neighbors (2 more radii + gap), plus transition and fit reach. Water bodies: the level radius.
    const double reach = 0.5*options().maxLineWidth().get()+options().maxBufferWidth().get();
    const double profile = options().smoothingDistance().get()*std::max(1u, options().smoothingPasses().get()) +
        2.0*options().sampleSpacing().get();
    double halo = reach + profile;
    if (water)
    {
        const double radius = options().waterLevelRadius().get();
        const double lattice = std::max(options().sampleSpacing().get(), radius/20.0);
        // Wet nodes lie within the lattice margin of the tile; flowing ones read channel surfaces within the
        // level radius of themselves, and those surfaces depend on two flow windows and the smoothing window.
        const double margin = std::max(radius, options().maxBufferWidth().get() + 2.0*lattice) + 3.0*lattice;
        halo = std::max(reach, margin + radius) + 2.0*options().flowWindow().get() + profile +
            options().smoothingDistance().get();
    }
    else if (options().roadFootprints().get())
        halo += options().junctionTransition().get();
    else
    {
        const double seed = 0.5*options().maxLineWidth().get() + options().junctionMargin().get();
        if (options().junctionDetection().get())
            halo = std::max(reach, 6.0*seed + 2.0*junctionGap + options().junctionTransition().get() +
                junctionFitReach) + profile;
        if (options().anchorHeight().isSet()) halo += options().maxAnchorRamp().get();
    }
    if (!featureQueries(key, source->getFeatureProfile(), halo, queries))
        return GeoHeightField(Status(Status::ConfigurationError, "Cannot transform flattening feature query"));
    osg::ref_ptr<Session> session = new Session(map);
    FilterContext context(session.get());
    // Retain read-only features for this job so their geometry survives each source cursor.
    std::map<FeatureID, std::vector<osg::ref_ptr<const Feature>>> seen;
    const bool absolute = options().profileElevation().get() == "absolute";
    unsigned group = 0;
    for (const auto& query : queries)
    {
        auto cursor = source->createFeatureCursor(query, {}, &context, progress);
        if (!cursor) return GeoHeightField::INVALID;
        while (cursor->hasMore())
        {
            if (canceled(progress)) return GeoHeightField::INVALID;
            auto* feature = cursor->nextFeature();
            if (!feature || !feature->getGeometry() || !feature->getSRS()) continue;
            // Keep distinct fragments with a shared ID; discard only matching geometry in the same SRS.
            if (feature->getFID() != 0)
            {
                auto& previous = seen[feature->getFID()];
                const bool duplicate = std::any_of(previous.begin(), previous.end(), [&](const auto& other)
                {
                    return (absolute ? feature->getSRS()->isEquivalentTo(other->getSRS()) :
                        feature->getSRS()->isHorizEquivalentTo(other->getSRS())) &&
                        sameGeometry(*feature->getGeometry(), *other->getGeometry(), absolute);
                });
                if (duplicate) continue;
                previous.emplace_back(feature);
            }
            const double width = options().lineWidth()->eval(feature, context);
            const double buffer = options().bufferWidth()->eval(feature, context);
            const double inset = options().footprintInset()->eval(feature, context);
            if (!std::isfinite(inset) || inset < 0.0 || inset > 1000.0 ||
                (inset > 0.0 && (!options().roadFootprints().get() || buffer != 0.0)))
                return GeoHeightField(Status(Status::ConfigurationError,
                    "Footprint inset requires paired profiles, no outer shoulder, and a width from zero to 1000 meters"));
            if (!std::isfinite(width) || !std::isfinite(buffer) || width < 0.0 || buffer < 0.0) continue;
            if (!compiler)
            {
                if (!makeGrid(extent, size, grid, world, lo, hi)) return GeoHeightField::INVALID;
                compiler.reset(new Compiler(options(), map->getElevationPool(), preceding, progress, lo, hi));
            }
            const double radius = 0.5*std::min(width, options().maxLineWidth().get());
            const double shoulder = std::min(buffer, options().maxBufferWidth().get());
            if (options().roadFootprints().get())
            {
                // Each feature explicitly pairs its complete centerlines with its pavement footprint.
                // Search width bounds profile support only; it never expands the pavement mask.
                auto profile = std::make_shared<Surface>(absolute ? nullptr : compiler->junctions, absolute);
                std::vector<const Polygon*> footprints;
                std::vector<const Geometry*> lines;
                ConstGeometryIterator components(feature->getGeometry(), false);
                while (components.hasMore())
                {
                    const auto* component = components.next();
                    if (component->getType() == Geometry::TYPE_POLYGON)
                        footprints.push_back(static_cast<const Polygon*>(component));
                    else if (component->getType() == Geometry::TYPE_LINESTRING)
                        lines.push_back(component);
                    else return GeoHeightField(Status(Status::ConfigurationError,
                        "Road footprints require only centerlines and polygons in each feature"));
                }
                if (lines.empty() || footprints.empty() || radius <= 0.0)
                    return GeoHeightField(Status(Status::ConfigurationError,
                        "Road footprints require centerlines, pavement polygons, and positive profile search width"));
                std::vector<std::string> ids;
                if (!options().profileIdAttribute()->empty())
                {
                    std::istringstream input(feature->getString(options().profileIdAttribute().get()));
                    std::string id;
                    while (std::getline(input, id, ',')) ids.push_back(id);
                    if (ids.size() != lines.size() || std::any_of(ids.begin(), ids.end(),
                        [](const std::string& value) { return value.empty(); }))
                        return GeoHeightField(Status(Status::ConfigurationError,
                            "Profile IDs must pair one nonempty ID with each complete centerline"));
                }
                const bool commonJunction = options().junction()->eval(feature, context);
                if (absolute && (commonJunction || lines.size() != 1))
                    return GeoHeightField(Status(Status::ConfigurationError,
                        "Absolute road footprints require one continuous solved centerline and no junction refit"));
                if (commonJunction && ids.empty())
                    return GeoHeightField(Status(Status::ConfigurationError, "Junction grading requires profile IDs"));
                std::vector<unsigned> members;
                const double offset = absolute ? options().profileOffset()->eval(feature, context) : 0.0;
                if (!std::isfinite(offset) || std::abs(offset) > 10000.0)
                    return GeoHeightField(Status(Status::ConfigurationError, "Invalid absolute profile offset"));
                for (std::size_t i = 0; i < lines.size(); ++i)
                {
                    unsigned member = group++;
                    if (!ids.empty())
                        member = compiler->profileGroups.emplace(ids[i], member).first->second;
                    members.push_back(member);
                    const bool ok = absolute ?
                        compiler->absoluteRoad(*lines[i], feature->getSRS(), radius, shoulder, member, profile.get(), offset) :
                        compiler->addLine(*lines[i], feature->getSRS(), radius, shoulder, member, profile);
                    if (!ok)
                        return GeoHeightField::INVALID;
                }
                if (commonJunction && !compiler->junction(lines, footprints, feature->getSRS(), members))
                    return GeoHeightField::INVALID;
                const unsigned owner = group++;
                for (const auto* footprint : footprints)
                    if (!compiler->polygon(*footprint, feature->getSRS(), shoulder, owner, profile, inset))
                        return GeoHeightField::INVALID;
                continue;
            }
            const double depth = water ? options().depth()->eval(feature, context) : 0.0;
            if (!std::isfinite(depth) || depth < 0.0 || depth > 100.0)
                return GeoHeightField(Status(Status::ConfigurationError, "Water depth must be from 0 to 100 meters"));
            double anchorHeight = NO_DATA_VALUE, anchorRamp = 0.0;
            if (options().anchorHeight().isSet())
            {
                anchorHeight = options().anchorHeight()->eval(feature, context);
                anchorRamp = std::min(options().anchorRamp()->eval(feature, context), options().maxAnchorRamp().get());
            }
            ConstGeometryIterator parts(feature->getGeometry(), false);
            while (parts.hasMore())
            {
                const auto* part = parts.next();
                bool ok = true;
                const bool point = part->getType() == Geometry::TYPE_POINTSET || part->getType() == Geometry::TYPE_POINT;
                if (part->getType() == Geometry::TYPE_LINESTRING)
                    ok = compiler->addLine(*part, feature->getSRS(), radius, shoulder, group++, nullptr, depth);
                else if (point && Grading::validHeight(anchorHeight) && anchorRamp > 0.0)
                {
                    for (const auto& p : *part)
                    {
                        osg::Vec3d world;
                        if (!compiler->position(p, feature->getSRS(), world)) { ok = false; break; }
                        compiler->addAnchor(world, anchorHeight, anchorRamp);
                    }
                }
                else if (point && !water)
                    ok = compiler->points(*part, feature->getSRS(), radius, shoulder, group);
                else if (part->getType() == Geometry::TYPE_POLYGON || part->getType() == Geometry::TYPE_RING)
                {
                    osg::ref_ptr<Geometry> owned = part->getType() == Geometry::TYPE_RING ?
                        part->cloneAs(Geometry::TYPE_POLYGON) : nullptr;
                    const auto* polygon = owned.valid() ? static_cast<const Polygon*>(owned.get()) :
                        static_cast<const Polygon*>(part);
                    ok = water ? compiler->waterPolygon(*polygon, feature->getSRS(), shoulder, depth,
                        options().flowing()->eval(feature, context)) :
                        compiler->polygon(*polygon, feature->getSRS(), shoulder, group++);
                }
                if (!ok) return canceled(progress) ? GeoHeightField::INVALID :
                    GeoHeightField(Status(Status::ResourceUnavailable,
                        "Grading input transform/sample failed, geometry exceeds 100 km, or sample budget exceeded"));
            }
        }
    }
    if (compiler && !compiler->profileNetwork(halo, water))
        return canceled(progress) ? GeoHeightField::INVALID : GeoHeightField(Status(Status::ResourceUnavailable,
            "Grading network transform/sample failed, geometry exceeds 100 km, or sample budget exceeded"));
    if (compiler && water && !compiler->buildWater(extent.transform(SpatialReference::get("wgs84")), group++))
        return canceled(progress) ? GeoHeightField::INVALID : GeoHeightField(Status(Status::ResourceUnavailable,
            "Water level lattice exceeds its node or sample budget, or sampling failed"));
    if (compiler && !compiler->emitNetwork(water)) return GeoHeightField::INVALID;
    if (!compiler || compiler->surface.empty() || canceled(progress)) return GeoHeightField::INVALID;
    std::vector<Grading::Constraint> constraints;
    std::vector<Grading::Candidate> scratch;
    scratch.reserve(8);
    std::vector<osg::Vec3d> natural;
    std::vector<std::size_t> pixels;
    auto hf = HeightFieldUtils::createReferenceHeightField(extent, size, size, 0u, false, NO_DATA_VALUE);
    bool wrote = false;
    for (unsigned row = 0; row < size; row += 8)
    for (unsigned col = 0; col < size; col += 8)
    {
        if (canceled(progress)) return GeoHeightField::INVALID;
        osg::Vec3d minimum(DBL_MAX, DBL_MAX, DBL_MAX), maximum(-DBL_MAX, -DBL_MAX, -DBL_MAX);
        const unsigned endRow = std::min(row+8, size), endCol = std::min(col+8, size);
        // Bound actual grid samples, rather than projecting just the corners of a nonlinear map projection.
        for (unsigned r = row; r < endRow; ++r)
        for (unsigned c = col; c < endCol; ++c)
        {
            const auto& p = world[r*size+c];
            for (unsigned axis = 0; axis < 3; ++axis)
            {
                minimum[axis] = std::min(minimum[axis], p[axis]);
                maximum[axis] = std::max(maximum[axis], p[axis]);
            }
        }
        if (!compiler->surface.intersects(minimum, maximum)) continue;
        for (unsigned r = row; r < endRow; ++r)
        for (unsigned c = col; c < endCol; ++c)
        {
            const std::size_t i = r*size+c;
            const auto constraint = compiler->surface.query(world[i], scratch);
            const bool exact = water ? constraint.bank == 0.0 : constraint.alpha == 1.0;
            if (exact && Grading::validHeight(constraint.height) && !cutOnly)
            {
                double value = constraint.height;
                // A channel bed is also bounded by nearby water bodies' bank cones.
                for (unsigned t = 0; t < constraint.termCount; ++t)
                {
                    const auto& term = constraint.terms[t];
                    value = std::min(value, value + term.alpha*(term.height+bankSlope*term.bank-value));
                }
                hf->getHeightList()[i] = static_cast<float>(value);
                wrote = true;
            }
            else if (constraint.alpha > 0.0)
            {
                pixels.push_back(i);
                natural.push_back(grid[i]);
                constraints.push_back(constraint);
            }
        }
    }
    if (!compiler->sample(natural)) return GeoHeightField::INVALID;
    for (std::size_t j = 0; j < pixels.size(); ++j)
    {
        const auto& c = constraints[j];
        if (Grading::validHeight(natural[j].z()) && Grading::validHeight(c.height))
        {
            const double h = natural[j].z();
            double result = h + c.alpha*((cutOnly ? std::min(c.height, h) : c.height)-h);
            if (water && c.termCount)
            {
                // Each nearby water feature bounds the bank by its own slope cone: the deepest cut and the highest
                // fill win, and both combine additively so the result stays continuous where features meet.
                double upper = h, lower = h;
                for (unsigned t = 0; t < c.termCount; ++t)
                {
                    const auto& term = c.terms[t];
                    upper = std::min(upper, h + term.alpha*(std::min(h, term.height+bankSlope*term.bank)-h));
                    lower = std::max(lower, h + term.alpha*(std::max(h, term.height-bankSlope*term.bank)-h));
                }
                result = cutOnly ? upper : upper+lower-h;
            }
            hf->getHeightList()[pixels[j]] = static_cast<float>(result);
            wrote = true;
        }
    }
    return wrote && !canceled(progress) ? GeoHeightField(hf, extent) : GeoHeightField::INVALID;
}
}
