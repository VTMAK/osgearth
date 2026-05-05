
#pragma include FrustumGrid.h.glsl

#define MAX_DECALS_PER_TILE 16


struct DecalTile
{
    uint count;
    uint indices[MAX_DECALS_PER_TILE];
    uint padding[3];
};

struct Decal
{
    mat4 mvm;
    mat4 mvmInverse;
    float a; // halfX extent (ortho) or zmin (persp)
    float b; // halfY extent (ortho) or zmax (persp)
    float c; // halfZ extent (ortho) or culling radius (persp)
    int textureIndex;
    float opacity;
    float distance; // > 0 = persp
    float tanHalfFovY;
    float aspect;
};


layout(binding = OE_BINDING_DECALS, std430) readonly buffer Decals
{
    Decal oe_decals[];
};

// each MAX_DECALS_PER_TILE+1 ints holds the count in the first Index followed
// by the decal indices in each subsequent index.
layout(binding = OE_BINDING_DECAL_TILES, std430) buffer DecalTiles
{
    DecalTile oe_decalTiles[];
};
