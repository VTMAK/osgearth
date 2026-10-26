# Terrain grading

The Prestige component is named `GradingLayer`, with factory tag `prestige:grading`.
It is intended for roads, building footprints, and inland water. The current algorithms include
terrain-derived road networks with detected junction planes, footprint-bounded road profiles, level pads,
and a water policy for channels, ponds, lakes, and river areas.
The older core `osgEarth::Contrib::FlatteningLayer` keeps its existing name and behavior.

`osgEarthPrestige::GradingLayer` generates terrain-derived road corridors and
level pads from GIS features. It is independent of the older
`osgEarth::Contrib::FlatteningLayer`. Enable the Prestige nodekit, include
`<osgEarthPrestige/GradingLayer>`, and link `osgEarth::osgEarthPrestige`.

For a whole-earth example, run `osgearth_imgui prestige_grading.earth` from
`tests` after setting up `osgearth_shell.bat` on Windows. The sample streams
ReadyMap imagery, elevation, and OSM roads, and opens at Stevens Canyon near
Mount Rainier. It infers road widths from OSM tags and uses twenty-meter shoulders,
five-meter profile samples, and thirty-meter smoothing. Bridges, tunnels, footpaths, and
highway areas are excluded. A Terrarium elevation alternative starts closed.

A **Prestige water grading** layer precedes the roads and grades OSM rivers,
streams, canals, ditches, ponds, lakes, and river areas from the same tiles.

In **Tools > Map**, toggle **Prestige road grading** or **Prestige water grading** to compare
terrain, or turn off **ReadyMap imagery** to inspect the shaded surface. **Tools > Viewpoints**
selects the road close-up, Rainier overview, equatorial roads, New Orleans, globe, Reflection
Lakes, and the Mississippi. The online OSM tiles contain clipped fragments. Fragments whose
clipped ends coincide are chained into one profile; overlapping (buffered) clips are not.

The implementation also corrects the shared tile-layer fallback conversion for
regional source profiles. A regional DEM's local maximum LOD must be converted
to the requesting map's LOD convention before selecting a parent tile. Without
this conversion, zooming beyond the DEM's native resolution can replace terrain
outside a flattened corridor with an excessively coarse source tile.

Load the nodekit and place the flattening layer **after** the elevation layers it
should modify:

```xml
<map>
    <libraries>osgEarthPrestige</libraries>

    <GDALElevation name="Terrain">
        <url>terrain.tif</url>
    </GDALElevation>

    <OGRFeatures name="Road centerlines">
        <url>roads.gpkg</url>
    </OGRFeatures>

    <prestige:grading name="Drivable terrain" features="Road centerlines">
        <line_width>10</line_width>
        <buffer_width>20</buffer_width>
        <sample_spacing>10</sample_spacing>
        <smoothing_distance>30</smoothing_distance>
        <elevation_resolution>5</elevation_resolution>
    </prestige:grading>
</map>
```

In C++, set the options and feature source before opening/adding the layer:

```cpp
osg::ref_ptr<osgEarthPrestige::GradingLayer> layer =
    new osgEarthPrestige::GradingLayer;
layer->setFeatureSource(roads);
layer->options().lineWidth() = osgEarth::Expression<double>(10.0);
layer->options().bufferWidth() = osgEarth::Expression<double>(20.0);
layer->options().smoothingDistance() = 30.0;
map->addLayer(layer);
```

## OSM road widths

The optional `osgEarthPrestige::OSMRoadWidthFilter` preprocesses OSM line features
in the feature source, before its tile cache. It writes `oe_road_width` in meters
and `oe_road_width_source` indicating the selected tag or fallback. Both the
flattening width expression and the road stroke read the same numeric attribute;
they do not run separate JavaScript width policies. The filter leaves original
OSM tags and geometry unchanged and does not annotate polygon features.

```xml
<filters>
    <!-- Place after any script that selects the surface roads to render. -->
    <prestige_osm_road_width attribute="oe_road_width" lane_width="3.5" max_width="40"/>
</filters>
```

Set `line_width` to `feature.properties.oe_road_width` and `max_line_width` to
at least the filter's `max_width`. Set the road style's `stroke-width` to that same
attribute expression and `stroke-width-units` to `meters`.

The first valid result wins:

1. `width`, then `width:carriageway`, then `est_width`, then a complete
   pipe-separated `width:lanes` list. Bare widths are meters; recognized length
   units and feet/inches such as `16'3"` are converted to meters.
2. `lanes` multiplied by `lane_width` (default 3.5 meters). If the total is
   absent or invalid, use `lanes:forward` plus `lanes:backward` and optional
   `lanes:both_ways`. Both directions must be supplied, except for an explicitly
   tagged one-way road. The total count is already for the whole way and is
   never doubled for two-way traffic.
3. Road-class estimates: one lane width for links; 5 meters for service roads;
   4.5 meters for living streets; 6 meters for residential/unclassified roads;
   and two lane widths for other road classes. These are demo assumptions,
   not measured widths or regional engineering standards.

Nonpositive, nonfinite, ambiguous, or oversized width tags fall through to the
next source. Lane counts must be whole numbers; missing lane widths are not
silently filled. Final fallback widths are bounded by `max_width`. The filter
ignores `maxwidth`, which limits vehicle width rather than measuring the road.
See the OSM definitions of [width](https://wiki.openstreetmap.org/wiki/Key:width)
and [lanes](https://wiki.openstreetmap.org/wiki/Key:lanes).

This preprocessing is specific to OSM roads and is independent of the generic
flattening layer. Building and airport polygons use their own footprints, and
the terrain transition buffer remains a separate setting. Inferring a road
width does not select an elevation policy for river polygons or other features.

## Surfaces

* **Lines and multilines:** lines sharing an end vertex where exactly two ends
  meet are chained into one path first, so way splits do not crease the profile.
  The layer samples the existing terrain along each path and smooths the longitudinal
  profile with a distance-weighted triangular filter (`smoothing_passes` times;
  two passes give gentler vertical curves). Each road is level across its full
  width: inside its half-width a road's own profile dominates, so a neighbor's
  shoulder cannot tilt it. Internal segment grades blend continuously at bends.
  A quintic shoulder transitions back to natural terrain with zero first and second
  derivatives at either edge of the transition.
* **Polygons, rings, and multipolygons:** flatten each component to a single
  terrain elevation sampled inside its footprint. The interior-point search
  handles concavity and holes. The shoulder also extends into holes; the part of
  a hole beyond the shoulder retains its original terrain.
* **Points and multipoints:** create round pads, using `line_width` as the pad
  diameter and sampling the terrain at each point.

Mixed multipart geometry is supported. Feature Z coordinates are ignored by default, and
source features are never modified. Distances are ground meters even when the
source or output uses geographic coordinates, Mercator, or projected feet.
Road stations, polygon anchors, and source sampling resolution are independent
of output tile boundaries and LOD. Centerline vertices are retained; intervals
are subdivided to at most `sample_spacing` meters.

## Detected junctions and anchors

In line mode (not paired footprints) a vertex shared by three or more road arms is a junction.
Only shared source vertices count, so roads that merely cross in plan (overpasses, underpasses)
never join. Each junction paves a disk of the widest member half-width plus `junction_margin`.
Junctions closer than their radii plus 8 meters are neighbors, for example a road crossing both
carriageways of a divided road; a junction's footprint covers its neighbors and the connections
between them. Its plane is a weighted least-squares fit to the members' smoothed profiles within
15 meters of that footprint. Member roads take the plane inside the footprint and return to their
own profiles over `junction_transition`, so intersections are free of undulations. Neighborhoods,
unlike transitive clusters, cannot grow without bound, so every paging job that can see a junction
fits the same plane. Paired footprints keep the compiler's own junctions instead; their line widths
are profile search widths, not pavement.

Point features become anchors when `anchor_height` is set: road profiles passing within one meter
of the point reach that ellipsoid height, ramping back to their own profile over `anchor_ramp`
meters with a quintic ease. Anchors prepare bridge, ramp, and overpass handoffs: a structure's end
height pins the ground road that meets it, and junction planes are fitted after anchoring.

## Inland water (water policy)

Set `policy=water` to grade OSM-style water features directly. Inside water the surface is set exactly
(pits fill, so water lies level). Outside it, banks within `buffer_width` are limited to `bank_slope`:
high ground is cut and low ground raised, so water never floats above its shore. Each nearby water
feature bounds the bank by its own slope cone; the deepest cut and the highest fill win.

* **Channels (lines):** the profile follows the digitized flow direction. Each sample takes the minimum
  over the preceding `flow_window` meters, removing bridge decks, culverts, and canopy shorter than the
  window, then is smoothed. The channel surface never sits more than `channel_cut_limit` below the
  smoothed terrain, so a line that strays up a valley side on coarse elevation cannot carve a canyon.
  The bed is the channel surface minus `depth`. Channels inside water bodies defer to the body, and
  a channel bed is bounded by nearby bodies' banks so streams descend into the lakes they enter.
* **Still bodies (polygons):** the union of all polygons is treated as one water shape, so source
  clipping does not matter. Terrain is sampled on a fixed global lattice (spacing
  `max(sample_spacing, water_level_radius/20)`); a node's level is the 0.2 quantile of wet samples
  within `water_level_radius`, minus `depth`. Bodies narrower than the radius are exactly level;
  larger ones follow their surface gently. Holes (islands) stay dry.
* **Flowing bodies (river and canal areas):** when `flowing` is true for a polygon, its surface takes
  the nearest channel line's surface minus `depth`: level across, falling with the river. Without a
  channel within the level radius the body is treated as still water.

Place water grading before road grading so roads sample the graded shoreline. The water policy
replaces neither the compiler-solved beds below nor bathymetry.

```xml
<prestige:grading name="Water grading" features="OSM inland water" min_level="14" max_data_level="14">
    <policy>water</policy>
    <line_width>({"river": 12, "canal": 8, "stream": 3})[feature.properties.waterway] || 4</line_width>
    <depth>({"river": 1.5, "canal": 1, "stream": 0.5})[feature.properties.waterway] || 0.6</depth>
    <flowing>feature.properties.water === "river" || feature.properties.waterway === "riverbank"</flowing>
    <buffer_width>50</buffer_width>
    <bank_slope>0.4</bank_slope>
    <max_line_width>20</max_line_width>
    <max_buffer_width>50</max_buffer_width>
</prestige:grading>
```

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `features` | required | Referenced or embedded feature source. |
| `road_footprints` | false | Require each feature to pair complete centerlines with pavement polygons. Polygons mask the sampled longitudinal profiles instead of becoming level pads. |
| `profile_elevation` | terrain | `terrain` samples existing elevation; `absolute` consumes solved centerline Z in road-footprint mode. |
| `profile_transition` | 0 | Absolute mode only: quintic fade into preceding elevation at both complete centerline endpoints, in meters; range 0–1000. Zero disables the fade. |
| `profile_offset` | 0 | Absolute mode only: numeric feature expression giving a receiver height offset in meters after datum conversion. Must be finite and within plus/minus 10,000 meters. |
| `footprint_inset` | 0 | Paired-footprint mode only: numeric expression for an inward bank transition, 0–1000 meters. Positive values require `buffer_width=0`. Shores and island boundaries stay unmodified. |
| `cut_only` | false | Sample preceding elevation even in footprint interiors and clamp the target to it before blending. Missing base elevation remains no-data. |
| `junction` | false | Boolean feature expression selecting common junction-plane fitting in road-footprint mode. |
| `profile_id_attribute` | empty | Attribute containing ordered, comma-separated centerline IDs. Required when junction fitting is selected; all road features must supply one ID per line component. |
| `junction_transition` | 30 | Distance in meters outside a junction footprint over which connected profiles blend back to their own grades; positive and at most 1,000. This does not expand terrain coverage. |
| `line_width` | 10 | Full road width or point-pad diameter, in meters. Supports numeric feature expressions. Zero disables line corridors. |
| `buffer_width` | 20 | Shoulder width on each side, in meters. Zero creates a hard boundary. Supports numeric feature expressions. |
| `max_line_width` | 100 | Upper bound on evaluated width; used to retrieve nearby features. |
| `max_buffer_width` | 100 | Upper bound on evaluated shoulder width; used to retrieve nearby features. |
| `sample_spacing` | 10 | Maximum spacing between road stations, in meters; range 1–100. |
| `smoothing_distance` | 30 | Physical radius of the longitudinal filter, in meters; range 0–1000. Zero disables profile smoothing. |
| `elevation_resolution` | 5 | Requested source-elevation resolution in meters; falls back to available coarser terrain. Must be positive. |
| `min_level` | 12 for the default global profile | First output LOD. Set this explicitly for coarse polygon work. Custom regional profiles retain the inherited LOD convention. |
| `policy` | road | `road` or `water`. Water is incompatible with `road_footprints`. |
| `smoothing_passes` | 2 | Repetitions of the profile filter, 1–8. |
| `junction_detection` | true | Line mode: detect junctions at vertices shared by three or more arms. |
| `junction_margin` | 1.5 | Detected junction radius beyond the widest member half-width, meters (0–100). |
| `anchor_height` | unset | Numeric expression; when set, point features are anchors at this ellipsoid height. |
| `anchor_ramp` | 60 | Numeric expression: anchor ramp length in meters, capped by `max_anchor_ramp`. |
| `max_anchor_ramp` | 200 | Upper bound of anchor ramps (0–2000 m); widens the query halo when anchors are enabled. |
| `depth` | 1 | Water: numeric expression for channel bed depth or body depth below its level, 0–100 m. |
| `flow_window` | 200 | Water channels: upstream minimum window in meters (0–5000). |
| `channel_cut_limit` | 4 | Water channels: deepest the channel surface may sit below the smoothed terrain, meters. |
| `water_level_radius` | 200 | Water bodies: level estimate radius in meters (10–2000). |
| `bank_slope` | 0.33 | Water: steepest bank rise per meter from the water's edge (0.01–10). |
| `flowing` | false | Water bodies: Boolean expression marking river or canal areas. |

## Compiler-owned road footprints

With `road_footprints=true`, each feature contains a GeometryCollection of one or more complete LineStrings
and one or more associated Polygons. Multipart lines/polygons are accepted in the default terrain mode.
In that mode profile heights come exclusively from terrain sampling along those centerlines;
footprint vertices and feature Z never supply elevations.
Polygons retain their holes, and their boundaries define the grading shoulder. A feature missing either
centerlines or polygons is rejected instead of silently flattening the polygon to one elevation.

In this mode `line_width` specifies a conservative **profile search diameter**, not pavement width.
It must cover the polygon and shoulder from its associated centerlines. `max_line_width` bounds it for
feature retrieval. Within each footprint the profile fully controls elevation; only the polygon shoulder
blends to the base. Without the optional junction policy, multiple associated approaches blend their
longitudinal profiles; that can form ridges where differently oriented roads impose different cross-slopes.
Disconnected or elevated roads must not be assigned to a ground footprint merely because they cross in plan.

For a common junction surface, set `junction` to a Boolean feature expression and `profile_id_attribute`
to the attribute containing ordered centerline IDs. Junction features and their approach-road features
must reference the same IDs and complete alignments. The layer fits one plane to terrain samples along
the junction's centerlines, within two meters of its pavement footprint. Station spacing is at most five
meters, with length weights to reduce digitization-density bias. A small regularization term stabilizes
collinear inputs without inventing an unconstrained transverse slope. Fewer than three valid samples
leave the original road-profile behavior in place.

Every identified member profile uses that same plane inside the junction. Outside it, a quintic transition
returns to the profile with zero blend derivative at both ends, matching height and slope at the handoffs.
Nearby junction influence bands blend continuously, with core interiors taking precedence over bands.
Disconnected roads are excluded by profile identity, not inferred from proximity. The correction changes
profile height only: existing pavement masks, holes, shoulders, and sparse no-data remain authoritative.
All controls are built before tile sampling and owned by that paging job. Fits use complete source
geometry, so changing output tiles does not change the fit.

These are terrain-derived planes, not surveyed grades. They can inherit implausible slopes from coarse
urban elevation. Closely spaced junction cores and complex interchanges still warrant explicit inspection;
this policy does not solve a regional engineering-grade optimization or superelevation design.

The road experiment's `road_grading.py` adapter pairs compiler source-way corridors and accepted junction
approaches with unchanged pavement coverage. Rounded fill retains shared network ownership; its nearest
ground-centerline height association is recorded explicitly as inferred.

For one native grading level, set `min_level` and `max_data_level` to the same value and leave `max_level`
unset. The three city scenes use global-geodetic LOD 16 and 257-square grids, approximately 1.2 meters
north-south sample spacing. Higher-detail terrain resamples these native tiles. No-data cells allow the
underlying elevation layer to participate. A finer tile entirely covered by parent grading can report no
new detail, allowing the terrain pager to retain its already graded parent. The terrain's 17-square render
grids and the elevation source's 257-square grids use different resolution mappings; native source LOD
does not directly equal visible terrain LOD.

```xml
<prestige:grading name="Road profiles" features="Compiled road profiles">
    <road_footprints>true</road_footprints>
    <junction>[grading_junction]</junction>
    <profile_id_attribute>source_way_ids</profile_id_attribute>
    <junction_transition>30</junction_transition>
    <line_width>[profile_search_width_m]</line_width>
    <max_line_width>500</max_line_width>
    <buffer_width>5</buffer_width>
    <max_buffer_width>5</max_buffer_width>
    <min_level>16</min_level>
    <max_data_level>16</max_data_level>
    <tile_size>257</tile_size>
</prestige:grading>
```

Evaluated widths are clamped to the declared maxima. Increase those maxima when
using wider roads or shoulders; each maximum is limited to 10,000 meters.
Negative or nonfinite feature widths skip that feature. The inherited elevation
tile size is respected (2–1024 samples per side, default 257).

## Solved road profiles and bridge approaches

Set `profile_elevation=absolute` with `road_footprints=true` when an offline compiler has already solved
the longitudinal profile. Each feature must contain one complete, continuous LineString and its paired
polygon footprint(s). Centerline Z supplies elevations in the feature SRS's vertical datum; without a
vertical datum, Z means ellipsoid meters. The layer converts to ellipsoid heights before applying the
optional `profile_offset`. Polygon Z still has no meaning.

Solved footprint interiors dominate neighboring grading shoulders. Their distance weighting uses a
one-centimeter regularizer, keeping exterior blends continuous while preventing a nearby higher road
from pulling up an otherwise constrained road interior. Overlapping interiors still blend equally;
terrain-sampled profiles and ordinary pads retain their existing weighting. Grid interpolation at a
deck handoff may still need a small compiler-owned support footprint around the boundary, clipped away
from lower receivers. This does not change the native elevation LOD or require a custom terrain mesh.

Absolute mode evaluates piecewise-linear height at the nearest centerline station. It does not resample,
smooth or refit the supplied Z values. Multiple line components, junction-plane fitting, invalid heights,
and coincident XY stations with different heights are rejected. Zero is a valid height. Complete source
geometry and source heights participate in duplicate detection. Supply a well-sampled local alignment;
individual segments longer than 100 km are rejected.

`profile_transition` fades the profile's influence at the complete alignment's endpoints, including
inside its footprint, so a bounded solved approach can rejoin terrain-derived roads. The polygon and
its shoulder still bound coverage. Outside the influence, cells remain no-data; shoulder and endpoint
transitions sample only the preceding elevation stack. This is an elevation overlay on the existing
regular grid, not a new terrain mesh.

The bridge experiment exports the same centerline XYZ for both the deck and approach grading. Only
the ground approaches and small landing support areas receive grading; the central span does not.
A separate `profile_offset` of -0.03 meters beneath the short deck overlap prevents coincident terrain
and deck surfaces. This receiver policy is recorded separately from the solved road height. It is not
a structural clearance or a blanket adjustment to the bridge geometry.

```xml
<prestige:grading name="Solved approaches" features="Paired absolute profiles">
    <road_footprints>true</road_footprints>
    <profile_elevation>absolute</profile_elevation>
    <profile_transition>10</profile_transition>
    <profile_offset>[profile_offset_m]</profile_offset>
    <line_width>80</line_width>
    <max_line_width>80</max_line_width>
    <buffer_width>5</buffer_width>
    <max_buffer_width>5</max_buffer_width>
    <min_level>16</min_level>
    <max_data_level>16</max_data_level>
</prestige:grading>
```

Keep stacked receivers separate. An elevated profile must not grade terrain along the entire span,
and a road beneath a bridge must retain its own profile and footprint. This mode consumes a solved
profile; it does not infer bridge clearances, solve elevated junctions, or assign water elevations.

## Inland water beds

Compiler-solved river and pond beds use the same paired-profile input. A river supplies a downstream
longitudinal profile; a pond supplies a constant reference. The compiler subtracts an explicit inferred
depth and exports bed XYZ separately from future water-reference XYZ. No water surface is rendered by
this elevation layer, and an OSM water footprint does not provide bathymetry.

Set `footprint_inset` to the desired bank width, `buffer_width=0`, and `cut_only=true`. Influence is zero
at the polygon boundary and rises smoothly to full bed depth inside the footprint. Interior rings receive
the same treatment, preserving islands. The quintic transition has zero first and second derivative at
both ends. Cut-only output is `base + influence * min(0, target-base)`, sampling the preceding elevation
stack throughout. No-data outside the footprint or in missing source terrain still falls through.

Place the water pass after road grading and exclude ground-road shoulders, landing supports, surface rail
corridors, and foundations in its masks. Bridge spans remain independent geometry. These masks are an
explicit compiler ownership decision; the layer does not infer transport structures from spatial overlap.
Complete connected river runs should be paired before export, so arbitrary OSM way splits do not become banks.
The prototype does not yet solve every confluence or differently attributed depth transition.

```xml
<prestige:grading name="Inland water beds" features="Compiled water beds">
    <road_footprints>true</road_footprints>
    <profile_elevation>absolute</profile_elevation>
    <footprint_inset>[bank_width_m]</footprint_inset>
    <cut_only>true</cut_only>
    <line_width>[profile_search_width_m]</line_width>
    <max_line_width>10000</max_line_width>
    <buffer_width>0</buffer_width>
    <max_buffer_width>0</max_buffer_width>
    <min_level>16</min_level>
    <max_data_level>16</max_data_level>
    <tile_size>257</tile_size>
</prestige:grading>
```

Both roads and beds retain one native LOD of regular elevation grids. A composed sample can differ slightly
outside the edited footprint when sparse parent fallback changes interpolation; cut-only constrains native
target values, not the error between two different composed sampling resolutions. Narrow channels may also
be unresolved in the visible terrain mesh. The Sion validation records these effects separately.

## Paging and accuracy

The layer queries a geographic halo large enough for everything that can change a
pixel of the tile: the maximum footprint, the smoothing window, detected junction
neighborhoods with their fit reach and transitions, anchor ramps, and, for water, two
flow windows and the level radius. Neighboring paging jobs therefore load identical
inputs and agree exactly along shared tile edges. The halo includes neighboring source
tiles and antimeridian splits. Road
geometry is indexed in double-precision Earth-centered coordinates. Raster
blocks without intersecting constraints are rejected before individual sample
queries. Only centerline/pad elevations and pixels in shoulders require source
elevation samples. Fully flattened pixels do not need a natural-terrain sample.

Each paging job owns its geometry, spatial index, and elevation working set. It
samples only open elevation layers preceding this layer in map order; an empty
source stack returns no tile. Cancellation or a processing failure publishes no
partial tile. Invalid profile samples are not interpolated through, and
unaffected pixels remain `NO_DATA_VALUE` so the underlying elevation stack is
preserved. The elevation pool's own fallback and missing-data resolution still
apply before this layer sees samples.

Output heights are above the ellipsoid; do not configure an output `vdatum`.
Source elevation layers may use vertical datums, which the elevation pool converts.

Persistent and layer-local output caching are disabled because the result
depends on other layers' elevation and on feature content. The elevation pool
still reuses source tiles. Configure the layer before opening it. After editing
source features, invalidate affected terrain using the application's normal
terrain-refresh mechanism; this layer does not watch external files for edits.

## Current limits

* Line fragments are chained only where their ends coincide within one centimeter.
  Distinct clipped fragments with a shared nonzero feature ID are retained; only
  matching XY geometry (including components and polygon holes) in an equivalent
  horizontal SRS is deduplicated. Water polygons act as their union, so clipped
  water bodies need no stitching.
* Water levels come from the terrain, not bathymetry or surveyed water surfaces. On
  coarse elevation a lake can sit several meters below its mapped shore; widen
  `buffer_width` so the bank slope can absorb it. Channels narrower than a few output
  samples alias on the native grid; grade water at a finer native LOD where that matters.
  The water lattice is not continuous across the antimeridian.
* Road segments follow short ellipsoid-projected chords. Individual input edges
  longer than 100 km are rejected; densify them upstream. Pad polygon vertices must
  lie within 100 km of their first vertex; water bodies of any size are accepted. Polygon edges use that component's
  local tangent plane and are intended for regional footprints such as buildings
  and airports.
* A tile request is limited to one million source samples. Transform failures,
  sampling failures, unsupported geometry sizes, and exceeded budgets return an
  error tile. Cancellation returns no tile.
* Terrain is a single-valued height field. Bridges, tunnels, and stacked roads
  require separate geometry. This layer smooths the road profile but does not
  enforce an engineering maximum grade, vertical-curve radius, crown, banking,
  drainage, or intersection topology.
* Terrain tessellation must resolve the desired road width. A narrow corridor
  cannot become drivable if the rendered terrain mesh spans it with a few large
  triangles. Source DEM resolution and elevation-pool quantization also limit
  achievable vertical accuracy.

## Validation

On Windows x64, `RelWithDebInfo`, the transient Google Benchmark harness measured
complete uncached 257-by-257 tile generation over the same warm, in-memory affine
terrain source. Each road contained 90 ten-meter segments in a 1024-meter square
tile. Both layers used ten-meter roads and twenty-meter shoulders. The new layer
used its default thirty-meter profile smoothing. Source elevation was limited to
LOD zero so both implementations sampled the same terrain data.

Median wall times from five repetitions, with at least 0.8 seconds per repetition:

| Roads in tile | Legacy | Prestige, smoothing enabled | Reduction |
| --- | ---: | ---: | ---: |
| 1 | 4.19 ms | 2.31 ms | 45% |
| 10 | 14.32 ms | 10.70 ms | 25% |
| 40 | 51.27 ms | 32.59 ms | 36% |

These are CPU tile-generation workloads, not network or frame-rate measurements.
Different feature density, source I/O, and terrain resolution will change the
result. The source and output caches were treated separately: source tiles were
warm, but every measured output tile was generated again.

The transient correctness harness passed 6,563 assertions across 13 cases:
affine-road agreement with the legacy layer; cross-road flatness; profile
smoothing; missing data; zero shoulders; cancellation before and during work;
empty source stacks; concurrent requests; shared tile edges; antimeridian edges;
point pads; polygon holes; continuous bends; random indexed-versus-exhaustive
queries; projected-foot coordinates; a buffered feature outside a tile at 75N;
nondefault tile size; and factory/configuration round trips. Tests and benchmarks
remain outside the repository, as required by the project's transient-test policy.

A live OSM regression reproduces the reported seam near Stevens Canyon at
46.76056 N, 121.65161 W. At output LOD 14, duplicate-ID filtering previously
discarded a road fragment in one parent tile, leaving 12 shared-edge samples
without flattening. Retaining distinct fragments makes both the coverage and
heights agree exactly at those samples; the corresponding LOD 18 comparison
also agrees exactly. A separate in-memory regression checks that distinct
same-ID fragments survive while repeated complete geometry is still removed.

A later transient harness for the network and water rewrite passed 408 assertions: flow-minimum
locality and bump removal; chaining through degree-two joints and directed water chaining; junction
detection (none for plan crossings) and bounded neighborhoods; water-shape union, holes, and banks;
exact pond levels; level road cross-sections on a cross slope (under 1 cm across, versus 0.64 m of
natural slope); a planar junction core (4 micrometers residual) straddling a tile edge with identical
seams; split ways matching unsplit ways; anchors; channel bumps removed without any downstream rise;
river areas following their channel; level, slope-limited ponds with identical results for whole and
clipped polygons and identical seams; and paired footprints whose profile footprint misses the tile.
Offline hillshades of Reflection Lakes, Stevens Canyon, and the NYC compiled roads were compared with
grading on and off.

The full existing suite could not complete: it contains a duplicate Chonk test
name, and, after temporarily distinguishing that name, the existing `glTF converts
channels and factors to shared PBR materials` test crashes with an access
violation. That crash also reproduces when the PBR test is run alone. The
temporary test-name change is not part of this implementation.
