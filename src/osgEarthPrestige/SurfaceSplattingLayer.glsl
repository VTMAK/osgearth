#pragma vp_name Prestige surface coordinates
#pragma vp_function oe_surface_vertex, vertex_view, 0.9

uniform mat4 oe_surface_viewToTexture;
out vec2 oe_surface_uv;
out vec3 oe_surface_viewPosition;

// Derive metric UVs from the rendered (already morphed) regular terrain vertices.
void oe_surface_vertex(inout vec4 vertex)
{
    oe_surface_uv = (oe_surface_viewToTexture * vertex).xy;
    oe_surface_viewPosition = vertex.xyz;
}

[break]
#pragma vp_name Prestige terrain surface material
#pragma vp_location fragment_coloring
#pragma import_defines(OE_USE_PBR, OE_IS_SHADOW_CAMERA, OE_IS_DEPTH_CAMERA, OE_IS_PICK_CAMERA)
#pragma include PBRMaterial.glsl

uniform sampler2D oe_surface_albedo;
uniform sampler2D oe_surface_normal;
uniform sampler2D oe_surface_pbr;
uniform sampler2D oe_surface_ao;
uniform sampler2D oe_surface_height;
uniform mat4 oe_surface_viewToTexture;
uniform float oe_surface_reliefDepth;
uniform vec4 oe_surface_factors;
uniform int oe_surface_flags;
uniform float oe_surface_normalStrength;
uniform vec2 oe_surface_slopeEdges;
in vec2 oe_surface_uv;
in vec3 oe_surface_viewPosition;
in vec3 oe_UpVectorView;
in vec3 vp_Normal;
#ifdef OE_USE_PBR
struct OE_PBR { float displacement, roughness, ao, metal; } oe_pbr;
#endif

// Parallax-occlusion relief: march the view ray down through the height map (1 = the terrain surface, 0 =
// oe_surface_reliefDepth meters below it) and return the texture coordinate where it first meets the relief.
// Texture-space steps come from the exact view-to-texture matrix, so the depth is metric at any texture period.
// The effect fades out between 5 and 20 mm pixels, where relief can no longer resolve.
vec2 oe_surface_relief(vec2 uv, vec2 dx, vec2 dy, vec3 N)
{
    // Pixel footprint in meters: the matrix scales view meters by 1/period.
    float footprint = max(length(dx), length(dy)) / length(oe_surface_viewToTexture[0].xyz);
    float weight = 1.0 - smoothstep(0.005, 0.02, footprint);
    if (weight <= 0.0) return uv;
    vec3 ray = normalize(oe_surface_viewPosition);
    // Grazing rays are limited, as usual for parallax mapping, to bound the search and its smearing.
    float descent = max(-dot(ray, N), 0.3);
    vec2 span = (mat3(oe_surface_viewToTexture) * ray).xy * (oe_surface_reliefDepth / descent);
    float steps = mix(8.0, 24.0, 1.0 - descent);
    vec2 delta = span / steps;
    float stepDepth = 1.0 / steps;
    vec2 current = uv;
    float depth = 0.0;
    float surface = 1.0 - textureGrad(oe_surface_height, current, dx, dy).r;
    for (int i = 0; i < 24 && depth < surface; ++i)
    {
        current += delta;
        depth += stepDepth;
        surface = 1.0 - textureGrad(oe_surface_height, current, dx, dy).r;
    }
    // Interpolate the crossing between the last step above the relief and the first below it.
    vec2 previous = current - delta;
    float after = surface - depth;
    float before = 1.0 - textureGrad(oe_surface_height, previous, dx, dy).r - (depth - stepDepth);
    float t = after / min(after - before, -1e-5);
    return mix(uv, mix(current, previous, clamp(t, 0.0, 1.0)), weight);
}

// Height-aware coverage: across the coverage ramp at a footprint edge, each texel of relief appears once the
// coverage passes a threshold set by its height, so unit tops (stones, pavers) stay whole and the joints give way
// first, instead of the whole surface fading. Full and empty coverage are unchanged. The layer beneath is not
// visible to this pass, so the other side of the edge blends as a flat surface.
float oe_surface_heightCoverage(float coverage, float height)
{
    const float width = 0.08;
    float threshold = mix(width, 1.0 - width, 1.0 - height);
    return smoothstep(threshold - width, threshold + width, coverage);
}

// Called by terrain after coverage morph/opacity and before first-layer or framebuffer blending. A composited layer
// (blended over the terrain color in-shader) weights its normal and oe_pbr by coverage; a GL-blended layer writes
// the material's own values and the framebuffer blend applies coverage once.
void oe_terrain_applyLayerMaterial(inout vec4 texel, bool composited)
{
#if !defined(OE_IS_SHADOW_CAMERA) && !defined(OE_IS_DEPTH_CAMERA) && !defined(OE_IS_PICK_CAMERA)
    // Steep terrain (a cut face the heightfield draws across a posting beside graded pavement) fades out of
    // coverage. The face normal comes from screen derivatives: the interpolated vertex normal smears across
    // the large triangles of a cut face.
    vec3 face = cross(dFdx(oe_surface_viewPosition), dFdy(oe_surface_viewPosition));
    float level = abs(dot(normalize(face), normalize(oe_UpVectorView)));
    texel.a *= smoothstep(oe_surface_slopeEdges.x, oe_surface_slopeEdges.y, level);
    // Derivatives of the undisplaced coordinate keep mip selection continuous across relief edges.
    vec2 dx = dFdx(oe_surface_uv), dy = dFdy(oe_surface_uv);
    vec2 uv = oe_surface_uv;
    float height = 0.5;
    if ((oe_surface_flags & 16) != 0)
    {
        uv = oe_surface_relief(uv, dx, dy, normalize(vp_Normal));
        height = textureGrad(oe_surface_height, uv, dx, dy).r;
        texel.a = oe_surface_heightCoverage(texel.a, height);
    }
    vec4 albedo = textureGrad(oe_surface_albedo, uv, dx, dy);
    if ((oe_surface_flags & 8) != 0)
        albedo.rgb = mix(1.055 * pow(max(albedo.rgb, 0.0), vec3(1.0/2.4)) - 0.055,
            12.92 * albedo.rgb, lessThanEqual(albedo.rgb, vec3(0.0031308)));
    texel.rgb = albedo.rgb;
    // Material alpha can only thin the compiler's footprint (cracks and wear in paint show the surface beneath);
    // it never enlarges it.
    texel.a *= albedo.a;
    float coverage = composited ? texel.a : 1.0;
    if ((oe_surface_flags & 1) != 0)
    {
        vec3 N = normalize(vp_Normal);
        vec3 p1 = dFdx(oe_surface_viewPosition), p2 = dFdy(oe_surface_viewPosition);
        vec2 t1 = dx, t2 = dy;
        vec3 q2 = cross(p2, N), q1 = cross(N, p1);
        vec3 T = q2*t1.x + q1*t2.x, B = q2*t1.y + q1*t2.y;
        float magnitude = max(dot(T,T), dot(B,B));
        float scale = inversesqrt(max(magnitude, 1e-20));
        vec3 bump = textureGrad(oe_surface_normal, uv, dx, dy).xyz*2.0-1.0;
        bump.xy *= oe_surface_normalStrength;
        vec3 perturbed = mat3(T*scale, B*scale, N)*bump;
        vp_Normal = normalize(mix(N, normalize(perturbed), coverage));
    }
#ifdef OE_USE_PBR
    vec4 packedValue = oe_pbr_default(oe_surface_factors.x);
    if ((oe_surface_flags & 2) != 0) packedValue = textureGrad(oe_surface_pbr, uv, dx, dy);
    float ao = (oe_surface_flags & 4) != 0 ? textureGrad(oe_surface_ao, uv, dx, dy).r : -1.0;
    vec4 material = oe_pbr_decode(packedValue, oe_surface_factors, ao);
    oe_pbr.roughness = mix(oe_pbr.roughness, material.g, coverage);
    oe_pbr.ao = mix(oe_pbr.ao, material.b, coverage);
    oe_pbr.metal = mix(oe_pbr.metal, material.a, coverage);
    oe_pbr.displacement = mix(oe_pbr.displacement, height, coverage);
#endif
#endif
}
