/* osgEarth
 * Copyright 2026 Pelican Mapping
 * MIT License
 */
#pragma once

#include <osgEarth/Geometry>
#include <osgEarth/rtree.h>
#include <osg/Matrixd>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace osgEarthPrestige { namespace Grading
{
    //! Accepts finite real elevations only; never interpolates the no-data sentinel.
    inline bool validHeight(double h)
    {
        return std::isfinite(h) && h != NO_DATA_VALUE;
    }

    //! Clamps a fraction without extrapolating.
    inline double unit(double t) { return std::max(0.0, std::min(1.0, t)); }

    //! Quintic smoothstep on [0, 1] with zero first and second derivatives at both ends.
    inline double ease(double t)
    {
        t = unit(t);
        return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
    }

    //! Returns compact quintic influence with zero first/second derivative at the shoulder edges.
    inline double influence(double distance, double radius, double buffer)
    {
        if (distance <= radius) return 1.0;
        if (buffer <= 0.0 || distance >= radius + buffer) return 0.0;
        return 1.0 - ease((distance - radius) / buffer);
    }

    //! Water bank influence at ``outside`` meters beyond the water's edge: full over the inner 60% of ``reach``,
    //! then a quintic fade to zero, so a slope-limited bank that would extend past the reach meets natural
    //! terrain smoothly instead of ending in a step.
    inline double bankInfluence(double outside, double reach)
    {
        if (reach <= 0.0) return outside <= 0.0 ? 1.0 : 0.0;
        return influence(outside, 0.6*reach, 0.4*reach);
    }

    struct Segment
    {
        osg::Vec3d a, b;
        double ha = NO_DATA_VALUE, hb = NO_DATA_VALUE;
        double radius = 0.0, buffer = 0.0;
        unsigned group = 0;
        bool extendStart = false, extendEnd = false;
        double stationA = 0.0, stationB = 0.0, profileLength = 0.0, transition = 0.0;
    };

    class Surface;

    //! Returns the squared local XY distance from p to the segment ab (Z ignored).
    inline double segmentDistance2(const osg::Vec3d& a, const osg::Vec3d& b, const osg::Vec3d& p)
    {
        const osg::Vec3d ab(b.x()-a.x(), b.y()-a.y(), 0.0), ap(p.x()-a.x(), p.y()-a.y(), 0.0);
        const double t = ab.length2() > 0.0 ? unit((ap*ab)/ab.length2()) : 0.0;
        return (ap-ab*t).length2();
    }

    //! Returns squared distance to every local XY boundary, including island holes, from either side.
    inline double boundaryDistance2(const osgEarth::Polygon* polygon, const osg::Vec3d& p)
    {
        double d2 = std::numeric_limits<double>::max();
        osgEarth::ConstGeometryIterator rings(polygon, true);
        while (rings.hasMore())
        {
            osgEarth::ConstSegmentIterator edges(rings.next(), true);
            while (edges.hasMore())
            {
                const auto& edge = edges.next();
                d2 = std::min(d2, segmentDistance2(edge.first, edge.second, p));
            }
        }
        return d2;
    }

    //! Returns squared distance outside a local XY footprint, including its holes; interiors return zero.
    inline double footprintDistance2(const osgEarth::Polygon* polygon, const osg::Vec3d& p)
    {
        return polygon->contains2D(p.x(), p.y()) ? 0.0 : boundaryDistance2(polygon, p);
    }

    //! A disk swept along a local XY segment; a disk when a == b. Used for detected junction footprints.
    struct Capsule
    {
        osg::Vec3d a, b;
        double radius = 0.0;

        //! Returns the XY distance from p to the capsule, zero inside.
        double distance(const osg::Vec3d& p) const
        {
            const osg::Vec3d ab(b.x()-a.x(), b.y()-a.y(), 0.0), ap(p.x()-a.x(), p.y()-a.y(), 0.0);
            const double t = ab.length2() > 0.0 ? unit((ap*ab)/ab.length2()) : 0.0;
            return std::max(0.0, (ap-ab*t).length() - radius);
        }
    };

    struct Junction
    {
        std::vector<osg::ref_ptr<const osgEarth::Polygon>> footprints;
        std::vector<Capsule> capsules;
        osg::Matrixd worldToLocal;
        osg::Vec3d plane; // Ellipsoid height = plane.z + plane.x*local.x + plane.y*local.y.
        std::vector<unsigned> members;
        double transition = 30.0;

        //! Returns squared distance to the union of local footprints and capsules, respecting polygon holes.
        double distance2(const osg::Vec3d& local) const
        {
            double result = std::numeric_limits<double>::max();
            for (const auto& polygon : footprints)
                result = std::min(result, footprintDistance2(polygon.get(), local));
            for (const auto& capsule : capsules)
            {
                const double d = capsule.distance(local);
                result = std::min(result, d*d);
            }
            return result;
        }
    };

    //! Shared per-job junction planes, built before querying and used only by explicitly identified profiles.
    class Junctions
    {
    public:
        //! Retains a fully fitted plane and conservative ECEF bounds, including its transition region.
        void add(const Junction& junction, const osg::Vec3d& lo, const osg::Vec3d& hi)
        {
            const unsigned index = static_cast<unsigned>(_junctions.size());
            _junctions.push_back(junction);
            _index.Insert(lo.ptr(), hi.ptr(), index);
        }

        //! Blends a member profile into a common plane with zero derivative at both ends of the transition.
        //! This changes height only; pavement masks and terrain shoulders retain exclusive coverage control.
        double apply(unsigned member, const osg::Vec3d& world, double original) const
        {
            double sum = 0.0, weight = 0.0, alpha = 0.0;
            _index.Search(world.ptr(), world.ptr(), [&](unsigned index)
            {
                const auto& j = _junctions[index];
                if (!std::binary_search(j.members.begin(), j.members.end(), member)) return RTREE_KEEP_SEARCHING;
                const auto p = world*j.worldToLocal;
                const double a = influence(std::sqrt(j.distance2(p)), 0.0, j.transition);
                if (a > 0.0)
                {
                    // An interior plane dominates neighboring transition bands without a discrete priority switch.
                    const double w = a / (1.0-a+1e-8);
                    sum += w*(j.plane.z()+j.plane.x()*p.x()+j.plane.y()*p.y());
                    weight += w;
                    alpha = std::max(alpha, a);
                }
                return RTREE_KEEP_SEARCHING;
            });
            return weight > 0.0 ? original + alpha*(sum/weight-original) : original;
        }

        //! Number of retained junctions (diagnostics and tests).
        std::size_t size() const { return _junctions.size(); }

    private:
        std::vector<Junction> _junctions;
        RTree<unsigned, double, 3> _index;
    };

    //! Returns the value at quantile q (0..1) of ``values``, sorting them in place; empty input returns NO_DATA.
    inline double quantile(std::vector<double>& values, double q)
    {
        if (values.empty()) return NO_DATA_VALUE;
        std::sort(values.begin(), values.end());
        const double position = unit(q)*(values.size()-1);
        const std::size_t i = static_cast<std::size_t>(position);
        const double t = position-i;
        return i+1 < values.size() ? values[i]*(1.0-t)+values[i+1]*t : values[i];
    }

    //! Union of water polygons in one local tangent frame. An edge index keeps containment and bank distance
    //! fast for very large or densely digitized bodies. Overlapping polygons act as their union, so fragments
    //! that a tiled source clipped (with or without a buffer) behave exactly like the complete body.
    class WaterShape
    {
    public:
        //! Adds one polygon as local XY rings (outer ring and holes, closed or open) with its outward bank width
        //! and bed depth. Containment is even-odd per polygon, so holes (islands) stay dry.
        void add(const std::vector<std::vector<osg::Vec3d>>& rings, double bank, double depth, bool flowing = false)
        {
            const unsigned polygon = static_cast<unsigned>(_banks.size());
            _banks.push_back(bank);
            _depths.push_back(depth);
            _flowing.push_back(flowing ? 1 : 0);
            _maxBank = std::max(_maxBank, bank);
            for (const auto& ring : rings)
            {
                for (std::size_t i = 0; i < ring.size(); ++i)
                {
                    const auto& a = ring[i];
                    const auto& b = ring[(i+1) % ring.size()];
                    if (a.x() == b.x() && a.y() == b.y()) continue;
                    const unsigned index = static_cast<unsigned>(_edges.size());
                    _edges.push_back(Edge{osg::Vec3d(a.x(), a.y(), 0.0), osg::Vec3d(b.x(), b.y(), 0.0), polygon});
                    double lo[2] = {std::min(a.x(), b.x()), std::min(a.y(), b.y())};
                    double hi[2] = {std::max(a.x(), b.x()), std::max(a.y(), b.y())};
                    _index.Insert(lo, hi, index);
                }
            }
        }

        //! Fills ``hits`` with the sorted indices of polygons containing a local point (half-open ray rule).
        void containing(const osg::Vec3d& p, std::vector<unsigned>& hits) const
        {
            hits.clear();
            const double lo[2] = {p.x(), p.y()}, hi[2] = {std::numeric_limits<double>::max(), p.y()};
            _index.Search(lo, hi, [&](unsigned index)
            {
                const auto& e = _edges[index];
                if ((e.a.y() > p.y()) != (e.b.y() > p.y()) &&
                    e.a.x() + (p.y()-e.a.y())*(e.b.x()-e.a.x())/(e.b.y()-e.a.y()) > p.x())
                    hits.push_back(e.polygon);
                return RTREE_KEEP_SEARCHING;
            });
            std::sort(hits.begin(), hits.end());
            // Keep the polygons whose ray crossing count is odd.
            std::size_t out = 0;
            for (std::size_t i = 0; i < hits.size();)
            {
                std::size_t j = i;
                while (j < hits.size() && hits[j] == hits[i]) ++j;
                if ((j-i) % 2u) hits[out++] = hits[i];
                i = j;
            }
            hits.resize(out);
        }

        //! Bank influence at a local point: one inside any polygon, otherwise the strongest bankInfluence()
        //! among polygons within their bank reach. ``distance2`` receives the squared distance to water.
        double coverage(const osg::Vec3d& p, std::vector<unsigned>& scratch, double& distance2) const
        {
            containing(p, scratch);
            distance2 = 0.0;
            if (!scratch.empty()) return 1.0;
            distance2 = std::numeric_limits<double>::max();
            double alpha = 0.0;
            if (_maxBank <= 0.0) return alpha;
            const double lo[2] = {p.x()-_maxBank, p.y()-_maxBank}, hi[2] = {p.x()+_maxBank, p.y()+_maxBank};
            _index.Search(lo, hi, [&](unsigned index)
            {
                const auto& e = _edges[index];
                const double d2 = segmentDistance2(e.a, e.b, p);
                distance2 = std::min(distance2, d2);
                alpha = std::max(alpha, bankInfluence(std::sqrt(d2), _banks[e.polygon]));
                return RTREE_KEEP_SEARCHING;
            });
            return alpha;
        }

        //! Reports whether any water edge lies within ``limit`` of a local point.
        bool near(const osg::Vec3d& p, double limit) const
        {
            bool hit = false;
            const double lo[2] = {p.x()-limit, p.y()-limit}, hi[2] = {p.x()+limit, p.y()+limit};
            _index.Search(lo, hi, [&](unsigned index)
            {
                const auto& e = _edges[index];
                if (segmentDistance2(e.a, e.b, p) <= limit*limit) { hit = true; return RTREE_STOP_SEARCHING; }
                return RTREE_KEEP_SEARCHING;
            });
            return hit;
        }

        //! Deepest bed depth among the given polygons (from containing()); zero when none.
        double depth(const std::vector<unsigned>& polygons) const
        {
            double result = 0.0;
            for (unsigned i : polygons) result = std::max(result, _depths[i]);
            return result;
        }

        //! True when any of the given polygons (from containing()) is flowing water.
        bool flowing(const std::vector<unsigned>& polygons) const
        {
            for (unsigned i : polygons) if (_flowing[i]) return true;
            return false;
        }

        //! Widest bank of any polygon.
        double maxBank() const { return _maxBank; }
        //! Number of polygons added.
        std::size_t size() const { return _banks.size(); }
        //! True when no polygon edges were added.
        bool empty() const { return _edges.empty(); }

    private:
        struct Edge { osg::Vec3d a, b; unsigned polygon; };
        std::vector<Edge> _edges;
        std::vector<double> _banks, _depths;
        std::vector<char> _flowing;
        double _maxBank = 0.0;
        RTree<unsigned, double, 2> _index;
    };

    //! Water surface over lattice nodes in one local frame. Still water: a node's level is a low quantile of the
    //! still samples (terrain minus bed depth) within ``radius``, so bodies narrower than the radius come out
    //! exactly level and larger ones follow their surface gently. Flowing water: a node whose nearest wet sample
    //! is flowing takes that sample's value (channel surface minus depth). A point blends nearby node levels with
    //! a compact smooth kernel. Nothing depends on how a source split the polygons.
    class WaterField
    {
    public:
        //! ``still`` and ``flowing`` hold local (x, y, value) samples; ``nodes`` holds local lattice positions that
        //! need a level. ``spacing`` is the lattice spacing; ``q`` the still-water quantile (0..1).
        WaterField(const std::vector<osg::Vec3d>& still, const std::vector<osg::Vec3d>& flowing,
            const std::vector<osg::Vec3d>& nodes, double radius, double spacing, double q) : _reach(2.0*spacing)
        {
            RTree<unsigned, double, 2> stillIndex, flowingIndex;
            for (unsigned i = 0; i < still.size(); ++i)
            {
                double lo[2] = {still[i].x(), still[i].y()}, hi[2] = {lo[0], lo[1]};
                stillIndex.Insert(lo, hi, i);
            }
            for (unsigned i = 0; i < flowing.size(); ++i)
            {
                double lo[2] = {flowing[i].x(), flowing[i].y()}, hi[2] = {lo[0], lo[1]};
                flowingIndex.Insert(lo, hi, i);
            }
            std::vector<double> values;
            for (const auto& node : nodes)
            {
                values.clear();
                double nearestStill = radius*radius, nearestFlowing = radius*radius, flowingValue = NO_DATA_VALUE;
                const double lo[2] = {node.x()-radius, node.y()-radius}, hi[2] = {node.x()+radius, node.y()+radius};
                stillIndex.Search(lo, hi, [&](unsigned i)
                {
                    const double dx = still[i].x()-node.x(), dy = still[i].y()-node.y(), d2 = dx*dx+dy*dy;
                    if (d2 <= radius*radius) { values.push_back(still[i].z()); nearestStill = std::min(nearestStill, d2); }
                    return RTREE_KEEP_SEARCHING;
                });
                flowingIndex.Search(lo, hi, [&](unsigned i)
                {
                    const double dx = flowing[i].x()-node.x(), dy = flowing[i].y()-node.y(), d2 = dx*dx+dy*dy;
                    if (d2 < nearestFlowing || (d2 == nearestFlowing && flowing[i].z() < flowingValue))
                    {
                        nearestFlowing = d2;
                        flowingValue = flowing[i].z();
                    }
                    return RTREE_KEEP_SEARCHING;
                });
                double level = NO_DATA_VALUE;
                if (validHeight(flowingValue) && (values.empty() || nearestFlowing <= nearestStill)) level = flowingValue;
                else if (!values.empty()) level = quantile(values, q);
                else continue;
                const unsigned index = static_cast<unsigned>(_levels.size());
                _levels.emplace_back(node.x(), node.y(), level);
                double at[2] = {node.x(), node.y()};
                _index.Insert(at, at, index);
            }
        }

        //! Water surface height at a local point, or NO_DATA where no node level lies within two spacings.
        double height(const osg::Vec3d& p) const
        {
            double sum = 0.0, weight = 0.0;
            const double lo[2] = {p.x()-_reach, p.y()-_reach}, hi[2] = {p.x()+_reach, p.y()+_reach};
            _index.Search(lo, hi, [&](unsigned i)
            {
                const double dx = _levels[i].x()-p.x(), dy = _levels[i].y()-p.y();
                const double t = (dx*dx+dy*dy)/(_reach*_reach);
                if (t < 1.0)
                {
                    const double w = (1.0-t)*(1.0-t);
                    sum += w*_levels[i].z();
                    weight += w;
                }
                return RTREE_KEEP_SEARCHING;
            });
            return weight > 0.0 ? sum/weight : NO_DATA_VALUE;
        }

        //! True when no node received a level.
        bool empty() const { return _levels.empty(); }

    private:
        double _reach;
        std::vector<osg::Vec3d> _levels;
        RTree<unsigned, double, 2> _index;
    };

    struct Pad
    {
        osg::ref_ptr<const osgEarth::Polygon> polygon;
        osg::Matrixd worldToLocal;
        double height = NO_DATA_VALUE, buffer = 0.0;
        double inset = 0.0;
        unsigned group = 0;
        // Optional worker-owned centerline profiles; the polygon supplies coverage, not a constant height.
        std::shared_ptr<const Surface> profile;
        bool useProfileAlpha = false;
        // Optional water union (water policy): replaces ``polygon`` coverage with union containment plus
        // per-polygon outward banks. ``worldToLocal`` is then the shape's frame.
        std::shared_ptr<const WaterShape> shape;
        // Optional water level field (water policy); replaces the constant height.
        std::shared_ptr<const WaterField> level;
    };

    struct Candidate
    {
        unsigned group;
        double distance2, height, alpha;
        double weight = 0.0;
        double distanceScale = 1.0;
        double bank = 0.0; // water banks: meters beyond the water's edge (zero inside)
        bool area = false; // water banks: a water body rather than a channel
    };

    //! One water feature's contribution to a bank pixel: its height at the edge, distance beyond the edge, and
    //! bank influence.
    struct WaterTerm
    {
        double height = NO_DATA_VALUE, bank = 0.0, alpha = 0.0;
    };

    struct Constraint
    {
        double height = NO_DATA_VALUE, alpha = 0.0;
        // Water banks only: meters beyond the nearest water edge (zero inside water); negative elsewhere.
        double bank = -1.0;
        // Water banks only: the nearest features, each bounding the bank (or a channel bed) by its own slope
        // cone, so a stream bed descends to the lake it enters instead of stepping into it.
        static const unsigned maxTerms = 4u;
        WaterTerm terms[maxTerms];
        unsigned termCount = 0;
    };

    //! Worker-owned immutable-after-build geometry. Query scratch belongs to the caller.
    class Surface
    {
    public:
        //! Retains optional shared junction controls; callers finish all mutations before querying this surface.
        //! With ``edgeWeights``, a segment's blending weight uses the distance beyond its pavement half-width, so
        //! inside a road its own level cross-section dominates any neighbor's shoulder.
        explicit Surface(std::shared_ptr<const Junctions> junctions = nullptr, bool nearestSegments = false,
            bool edgeWeights = false) :
            _nearestSegments(nearestSegments), _edgeWeights(edgeWeights), _junctions(junctions) { }

        //! Switches to water banks: outside a channel's half-width or a water body, candidates keep the water
        //! height and report their distance to the edge, so the caller can build a slope-limited bank; influence
        //! uses bankInfluence() over the buffer (bank reach) instead of a shoulder fade.
        void setWaterBanks(bool enabled) { _water = enabled; }

        //! Adds a metric ECEF segment, including a zero-length segment for a point pad.
        void add(const Segment& segment)
        {
            if (!validHeight(segment.ha) || !validHeight(segment.hb)) return;
            const unsigned index = static_cast<unsigned>(_segments.size());
            IndexedSegment indexed;
            indexed.segment = segment;
            indexed.ab = segment.b - segment.a;
            indexed.inverseLength2 = indexed.ab.length2() > 0.0 ? 1.0/indexed.ab.length2() : 0.0;
            indexed.outer2 = (segment.radius+segment.buffer)*(segment.radius+segment.buffer);
            _segments.push_back(indexed);
            double lo[3], hi[3];
            const double r = segment.radius + segment.buffer;
            for (unsigned axis = 0; axis < 3; ++axis)
            {
                lo[axis] = std::min(segment.a[axis], segment.b[axis]) - r;
                hi[axis] = std::max(segment.a[axis], segment.b[axis]) + r;
            }
            _index.Insert(lo, hi, index);
        }

        //! Adds an owned local polygon and conservative ECEF bounds; holes remain unflattened interiors.
        //! A profile may still be empty (its segments are emitted later); queries then skip it.
        void add(const Pad& pad, const osg::Vec3d& lo, const osg::Vec3d& hi)
        {
            if (!pad.profile && (pad.level ? pad.level->empty() : !validHeight(pad.height))) return;
            const unsigned index = static_cast<unsigned>(_pads.size());
            _pads.push_back(pad);
            _padIndex.Insert(lo.ptr(), hi.ptr(), index);
        }

        //! Reports whether any valid terrain constraints were compiled.
        bool empty() const { return _segments.empty() && _pads.empty(); }

        //! Conservatively rejects whole raster blocks before doing individual point queries.
        bool intersects(const osg::Vec3d& lo, const osg::Vec3d& hi) const
        {
            bool hit = false;
            _index.Search(lo.ptr(), hi.ptr(), [&hit](unsigned) { hit = true; return RTREE_STOP_SEARCHING; });
            if (!hit && !_pads.empty())
                _padIndex.Search(lo.ptr(), hi.ptr(), [&hit](unsigned) { hit = true; return RTREE_STOP_SEARCHING; });
            return hit;
        }

        //! Queries intersecting bounds and continuously blends nearby segments within each geometry part.
        //! Shared endpoints are never deleted. Missing endpoint pairs are omitted when building.
        Constraint query(const osg::Vec3d& point, std::vector<Candidate>& scratch) const
        {
            scratch.clear();
            _index.Search(point.ptr(), point.ptr(), [&](unsigned index)
            {
                const auto& indexed = _segments[index];
                const auto& s = indexed.segment;
                const double projected = ((point-s.a)*indexed.ab)*indexed.inverseLength2;
                const double t = unit(projected);
                // Continue internal segment grades through joins. Collinear segments then agree exactly,
                // while weighted corner contributions avoid jumps across the nearest-segment bisector.
                const double heightT = (projected < 0.0 && s.extendStart) || (projected > 1.0 && s.extendEnd) ?
                    projected : t;
                const double d2 = (point - (s.a + indexed.ab*t)).length2();
                if (d2 <= indexed.outer2)
                {
                    const double d = std::sqrt(d2);
                    double alpha = influence(d, s.radius, s.buffer);
                    const double outside = std::max(0.0, d - s.radius);
                    if (_water) alpha = bankInfluence(outside, s.buffer);
                    if (s.transition > 0.0)
                    {
                        const double station = s.stationA + (s.stationB-s.stationA)*t;
                        alpha *= 1.0-influence(std::min(station, s.profileLength-station), 0.0, s.transition);
                    }
                    Candidate value{s.group, d2, s.ha + (s.hb-s.ha)*heightT, alpha};
                    value.bank = outside;
                    if (_edgeWeights)
                    {
                        // Distance beyond the pavement edge: zero across the whole road, so its level section
                        // wins inside, while a 0.1 m scale keeps the blend continuous just outside the edge.
                        value.distance2 = outside*outside;
                        value.distanceScale = 100.0;
                    }
                    if (_nearestSegments) accumulateNearest(scratch, value);
                    else accumulate(scratch, value);
                }
                return RTREE_KEEP_SEARCHING;
            });
            if (!_pads.empty()) _padIndex.Search(point.ptr(), point.ptr(), [&](unsigned index)
            {
                const auto& pad = _pads[index];
                const auto p = point * pad.worldToLocal;
                double d2 = 0.0, alpha = 0.0;
                if (pad.shape)
                {
                    std::vector<unsigned> hits;
                    alpha = pad.shape->coverage(p, hits, d2);
                }
                else
                {
                    d2 = footprintDistance2(pad.polygon.get(), p);
                    alpha = influence(std::sqrt(d2), 0.0, pad.buffer);
                    if (pad.inset > 0.0)
                        alpha = pad.polygon->contains2D(p.x(), p.y()) ?
                            1.0-influence(std::sqrt(boundaryDistance2(pad.polygon.get(), p)), 0.0, pad.inset) : 0.0;
                }
                double height = pad.height;
                if (alpha > 0.0 && pad.level)
                    height = pad.level->height(p);
                if (alpha > 0.0 && pad.profile)
                {
                    std::vector<Candidate> profileScratch;
                    const auto profile = pad.profile->query(point, profileScratch);
                    height = profile.height;
                    if (pad.useProfileAlpha) alpha *= profile.alpha;
                }
                if (validHeight(height))
                {
                    Candidate candidate{pad.group, d2, height, alpha};
                    candidate.bank = std::sqrt(d2);
                    candidate.area = static_cast<bool>(pad.shape);
                    // Solved road footprints own their interior elevation. A neighboring shoulder must
                    // approach zero weight at that boundary, rather than lifting the road by decimeters.
                    // A one-centimeter distance regularizer keeps the blend continuous outside the mask;
                    // overlapping interiors still share equal weight. Terrain-derived policies are unchanged.
                    if (pad.useProfileAlpha) candidate.distanceScale = 10000.0;
                    accumulate(scratch, candidate);
                }
                return RTREE_KEEP_SEARCHING;
            });
            Constraint result;
            if (_water)
            {
                // Inside a water body only the bodies count, so channel beds cannot dimple a lake. Inside a channel
                // only channel cores set the bed, while other features' banks still bound it (see Constraint).
                double areaSum = 0.0, areaWeight = 0.0, coreSum = 0.0, coreWeight = 0.0;
                for (const auto& candidate : scratch)
                {
                    if (candidate.alpha <= 0.0 || candidate.bank > 0.0) continue;
                    const double w = candidate.alpha / (1.0 + candidate.distanceScale*candidate.distance2);
                    (candidate.area ? areaSum : coreSum) += w*candidate.height;
                    (candidate.area ? areaWeight : coreWeight) += w;
                }
                if (areaWeight > 0.0)
                {
                    result.height = areaSum/areaWeight;
                    result.alpha = 1.0;
                    result.bank = 0.0;
                    return result;
                }
                // Banks keep each feature separate; the caller bounds the terrain by every feature's cone.
                std::sort(scratch.begin(), scratch.end(), [](const Candidate& a, const Candidate& b)
                {
                    return std::make_tuple(a.bank, a.group) < std::make_tuple(b.bank, b.group);
                });
                for (const auto& candidate : scratch)
                {
                    if (candidate.alpha <= 0.0 || candidate.bank <= 0.0 || result.termCount == Constraint::maxTerms)
                        continue;
                    auto& term = result.terms[result.termCount++];
                    term.height = candidate.height;
                    term.bank = candidate.bank;
                    term.alpha = candidate.alpha;
                    result.alpha = std::max(result.alpha, candidate.alpha);
                }
                if (coreWeight > 0.0)
                {
                    result.height = coreSum/coreWeight;
                    result.alpha = 1.0;
                    result.bank = 0.0;
                }
                else if (result.termCount)
                {
                    result.height = result.terms[0].height;
                    result.bank = result.terms[0].bank;
                }
                return result;
            }
            double sum = 0.0, weight = 0.0;
            // Normalize per part so digitization density does not give one road priority at a junction.
            for (const auto& candidate : scratch)
            {
                const double w = candidate.alpha / (1.0 + candidate.distanceScale*candidate.distance2);
                const double h = _junctions ? _junctions->apply(candidate.group, point, candidate.height) : candidate.height;
                sum += w * h;
                weight += w;
                result.alpha = std::max(result.alpha, candidate.alpha);
            }
            if (weight > 0.0) result.height = sum / weight;
            return result;
        }

    private:
        //! Solved profiles use the nearest station instead of smoothing across unrelated station heights.
        //! Ties retain the first segment; valid continuous paths agree at their shared endpoint.
        static void accumulateNearest(std::vector<Candidate>& candidates, const Candidate& value)
        {
            // Retain zero-weight endpoints so a more distant positive segment cannot defeat the fade.
            for (auto& candidate : candidates)
            {
                if (candidate.group != value.group) continue;
                if (value.distance2 < candidate.distance2) candidate = value;
                return;
            }
            candidates.push_back(value);
        }

        //! Blends compact segment contributions within a part and retains its maximum flattening influence.
        static void accumulate(std::vector<Candidate>& candidates, const Candidate& value)
        {
            if (value.alpha <= 0.0) return;
            const double weight = value.alpha / (1.0 + value.distanceScale*value.distance2);
            for (auto& candidate : candidates)
            {
                if (candidate.group != value.group) continue;
                candidate.height += (value.height-candidate.height)*weight/(candidate.weight+weight);
                candidate.weight += weight;
                candidate.distance2 = std::min(candidate.distance2, value.distance2);
                candidate.alpha = std::max(candidate.alpha, value.alpha);
                candidate.bank = std::min(candidate.bank, value.bank);
                return;
            }
            candidates.push_back(value);
            candidates.back().weight = weight;
        }

        struct IndexedSegment
        {
            Segment segment;
            osg::Vec3d ab;
            double inverseLength2, outer2;
        };
        std::vector<IndexedSegment> _segments;
        bool _nearestSegments = false;
        bool _edgeWeights = false;
        bool _water = false;
        std::shared_ptr<const Junctions> _junctions;
        std::vector<Pad> _pads;
        RTree<unsigned, double, 3> _index, _padIndex;
    };

    // ------------------------------------------------------------------------------------------------------------
    // Network algorithms: pure functions over world-space polylines, shared by the layer and its tests.

    //! One source polyline in ECEF on the reference ellipsoid, with its grading parameters.
    struct NetworkLine
    {
        std::vector<osg::Vec3d> world;
        double radius = 0.0, buffer = 0.0;
        unsigned group = 0;
        double depth = 0.0; // water channels: bed depth below the downhill-filtered profile
    };

    //! A line inside a chained path; ``reversed`` walks the line from its last vertex to its first.
    struct PathPiece
    {
        unsigned line = 0;
        bool reversed = false;
    };

    //! Integer key of a world position on a 1 cm grid; coincident source nodes produce the same key.
    struct VertexKey
    {
        std::int64_t x, y, z;
        bool operator==(const VertexKey& o) const { return x == o.x && y == o.y && z == o.z; }
    };

    //! Hash for VertexKey.
    struct VertexKeyHash
    {
        std::size_t operator()(const VertexKey& k) const
        {
            std::size_t h = std::hash<std::int64_t>()(k.x);
            h ^= std::hash<std::int64_t>()(k.y) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            h ^= std::hash<std::int64_t>()(k.z) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            return h;
        }
    };

    //! Quantizes an ECEF position to the 1 cm vertex grid.
    inline VertexKey vertexKey(const osg::Vec3d& p)
    {
        return VertexKey{std::llround(p.x()*100.0), std::llround(p.y()*100.0), std::llround(p.z()*100.0)};
    }

    //! Counts arms per vertex: an interior occurrence contributes two, an open line end one.
    inline std::unordered_map<VertexKey, unsigned, VertexKeyHash> vertexDegrees(const std::vector<NetworkLine>& lines)
    {
        std::unordered_map<VertexKey, unsigned, VertexKeyHash> degree;
        for (const auto& line : lines)
        {
            const std::size_t n = line.world.size();
            for (std::size_t i = 0; i < n; ++i)
                degree[vertexKey(line.world[i])] += (i > 0 ? 1u : 0u) + (i+1 < n ? 1u : 0u);
        }
        return degree;
    }

    //! Chains lines end to end through vertices where exactly two line ends meet and nothing else passes.
    //! With ``directed`` (water), only an end can join a start, so a chain keeps the digitized flow direction.
    //! Every line appears in exactly one path; closed rings and loops become single paths.
    inline std::vector<std::vector<PathPiece>> chainLines(const std::vector<NetworkLine>& lines, bool directed)
    {
        const auto degree = vertexDegrees(lines);
        // Line ends by vertex: (line index, 0 = start / 1 = end).
        std::unordered_map<VertexKey, std::vector<std::pair<unsigned, int>>, VertexKeyHash> ends;
        for (unsigned i = 0; i < lines.size(); ++i)
        {
            if (lines[i].world.size() < 2) continue;
            ends[vertexKey(lines[i].world.front())].emplace_back(i, 0);
            ends[vertexKey(lines[i].world.back())].emplace_back(i, 1);
        }
        //! Finds the line end joined to (line, end) at its vertex, or (-1, -1) when the vertex does not join.
        auto partner = [&](unsigned line, int end) -> std::pair<int, int>
        {
            const auto& p = end == 0 ? lines[line].world.front() : lines[line].world.back();
            const auto key = vertexKey(p);
            const auto d = degree.find(key);
            const auto e = ends.find(key);
            if (d == degree.end() || e == ends.end() || d->second != 2u || e->second.size() != 2u)
                return std::make_pair(-1, -1);
            const auto& other = e->second[0].first == line && e->second[0].second == end ? e->second[1] : e->second[0];
            if (other.first == line && other.second == end) return std::make_pair(-1, -1);
            if (directed && other.second == end) return std::make_pair(-1, -1); // two starts or two ends meet
            return std::make_pair(int(other.first), other.second);
        };
        std::vector<std::vector<PathPiece>> paths;
        std::vector<char> used(lines.size(), 0);
        for (unsigned seed = 0; seed < lines.size(); ++seed)
        {
            if (used[seed] || lines[seed].world.size() < 2) continue;
            // Walk backward from the seed's start to the head of its chain (or around a loop back to the seed).
            unsigned head = seed;
            bool headReversed = false;
            for (std::size_t guard = 0; guard <= lines.size(); ++guard)
            {
                const int entryEnd = headReversed ? 1 : 0;
                const auto previous = partner(head, entryEnd);
                if (previous.first < 0 || unsigned(previous.first) == seed || used[previous.first]) break;
                // The previous line must exit at this vertex: forward if its end joins, reversed if its start does.
                head = unsigned(previous.first);
                headReversed = previous.second == 0;
            }
            std::vector<PathPiece> path;
            unsigned current = head;
            bool reversed = headReversed;
            for (std::size_t guard = 0; guard <= lines.size(); ++guard)
            {
                used[current] = 1;
                path.push_back(PathPiece{current, reversed});
                const auto next = partner(current, reversed ? 0 : 1);
                if (next.first < 0 || used[next.first]) break;
                current = unsigned(next.first);
                reversed = next.second == 1;
            }
            paths.push_back(path);
        }
        return paths;
    }

    //! Symmetric triangular filter over one run of stations, integration-weighted so dense digitization does not
    //! dominate. ``passes`` > 1 repeats it, which approximates a smoother (B-spline-like) vertical curve.
    inline void smoothProfile(const std::vector<double>& station, std::vector<double>& heights,
        std::size_t first, std::size_t end, double distance, unsigned passes)
    {
        if (end <= first+1 || distance <= 0.0) return;
        std::vector<double> source;
        for (unsigned pass = 0; pass < passes; ++pass)
        {
            source.assign(heights.begin()+first, heights.begin()+end);
            for (std::size_t i = first; i < end; ++i)
            {
                double sum = 0.0, weight = 0.0;
                auto beginIt = std::lower_bound(station.begin()+first, station.begin()+end, station[i]-distance);
                for (std::size_t j = beginIt-station.begin(); j < end; ++j)
                {
                    const double d = std::abs(station[j]-station[i]);
                    if (station[j] > station[i]+distance) break;
                    const double w = std::max(0.0, 1.0-d/distance);
                    const double cell = 0.5*(station[std::min(j+1, end-1)]-station[j > first ? j-1 : first]);
                    sum += source[j-first]*w*cell;
                    weight += w*cell;
                }
                heights[i] = weight > 0.0 ? sum/weight : source[i-first];
            }
        }
    }

    //! Water flows from the first sample to the last. Replaces each height by the minimum over the preceding
    //! ``window`` meters of station, so bumps shorter than the window (bridge decks, culverts, canopy in a surface
    //! model) vanish and the bed cannot rise within a window. Unlike a running minimum from the source, the
    //! result is local: paging jobs that load different upstream extents still agree.
    inline void flowMinimum(const std::vector<double>& station, std::vector<double>& heights,
        std::size_t first, std::size_t end, double window)
    {
        if (end <= first+1 || window <= 0.0) return;
        const std::vector<double> source(heights.begin()+first, heights.begin()+end);
        std::deque<std::size_t> queue; // candidate indices with increasing source heights
        for (std::size_t i = first; i < end; ++i)
        {
            while (!queue.empty() && source[queue.back()-first] >= source[i-first]) queue.pop_back();
            queue.push_back(i);
            while (station[queue.front()] < station[i]-window) queue.pop_front();
            heights[i] = source[queue.front()-first];
        }
    }

    //! Detected junction: the shared vertex, its member groups and the radius it paves.
    struct JunctionSeed
    {
        osg::Vec3d world;
        std::vector<unsigned> members;
        double radius = 0.0;
    };

    //! Finds vertices where three or more road arms meet, with their member groups and a paving radius
    //! (widest member half-width plus ``margin``). Only shared source vertices create junctions, so roads that
    //! merely cross in plan (overpasses, underpasses) never do.
    inline std::vector<JunctionSeed> findJunctions(const std::vector<NetworkLine>& lines, double margin)
    {
        const auto degree = vertexDegrees(lines);
        std::unordered_map<VertexKey, JunctionSeed, VertexKeyHash> seeds;
        for (const auto& line : lines)
        {
            for (const auto& p : line.world)
            {
                const auto key = vertexKey(p);
                const auto d = degree.find(key);
                if (d == degree.end() || d->second < 3u) continue;
                auto& seed = seeds[key];
                seed.world = p;
                seed.members.push_back(line.group);
                seed.radius = std::max(seed.radius, line.radius + margin);
            }
        }
        std::vector<JunctionSeed> result;
        result.reserve(seeds.size());
        for (auto& entry : seeds)
        {
            auto& seed = entry.second;
            std::sort(seed.members.begin(), seed.members.end());
            seed.members.erase(std::unique(seed.members.begin(), seed.members.end()), seed.members.end());
            result.push_back(seed);
        }
        // Deterministic order independent of hash iteration.
        std::sort(result.begin(), result.end(), [](const JunctionSeed& a, const JunctionSeed& b)
        {
            return std::make_tuple(a.world.x(), a.world.y(), a.world.z()) <
                std::make_tuple(b.world.x(), b.world.y(), b.world.z());
        });
        return result;
    }

    //! Reports whether two junction seeds are close enough (sum of radii plus ``gap``) to share a plane,
    //! e.g. a road crossing both carriageways of a divided road.
    inline bool neighbors(const JunctionSeed& a, const JunctionSeed& b, double gap)
    {
        return (a.world-b.world).length() < a.radius+b.radius+gap;
    }

    //! Returns seed ``index`` and its neighbors, ascending. Each seed's junction is built from its own
    //! neighborhood rather than a transitive cluster: a chain of close junctions cannot grow without bound,
    //! so the result depends only on nearby geometry and adjacent paging jobs agree. Neighborhoods that are
    //! complete (every pair close) are identical for all their members and yield exactly one shared plane.
    inline std::vector<unsigned> neighborhood(const std::vector<JunctionSeed>& seeds, unsigned index, double gap)
    {
        std::vector<unsigned> result;
        for (unsigned i = 0; i < seeds.size(); ++i)
            if (i == index || neighbors(seeds[index], seeds[i], gap)) result.push_back(i);
        return result;
    }

    //! Weighted least-squares plane h = z + x*X + y*Y through local samples (x, y, h); a tiny ridge term keeps
    //! collinear inputs stable without inventing a cross-slope. Returns false when there is nothing to fit.
    inline bool fitPlane(const std::vector<osg::Vec3d>& samples, const std::vector<double>& weights, osg::Vec3d& plane)
    {
        double weight = 0.0, x = 0.0, y = 0.0, h = 0.0;
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            weight += weights[i];
            x += weights[i]*samples[i].x(); y += weights[i]*samples[i].y(); h += weights[i]*samples[i].z();
        }
        if (samples.size() < 3 || weight <= 0.0) return false;
        x /= weight; y /= weight; h /= weight;
        double xx = 0.0, xy = 0.0, yy = 0.0, xh = 0.0, yh = 0.0;
        for (std::size_t i = 0; i < samples.size(); ++i)
        {
            const double dx = samples[i].x()-x, dy = samples[i].y()-y, dh = samples[i].z()-h;
            xx += weights[i]*dx*dx; xy += weights[i]*dx*dy; yy += weights[i]*dy*dy;
            xh += weights[i]*dx*dh; yh += weights[i]*dy*dh;
        }
        const double regularization = std::max(1e-8, (xx+yy)*1e-8);
        xx += regularization; yy += regularization;
        const double determinant = xx*yy-xy*xy;
        if (!(determinant > 0.0)) return false;
        plane.x() = (xh*yy-yh*xy)/determinant;
        plane.y() = (yh*xx-xh*xy)/determinant;
        plane.z() = h-plane.x()*x-plane.y()*y;
        return std::isfinite(plane.x()) && std::isfinite(plane.y()) && std::isfinite(plane.z());
    }
} }
