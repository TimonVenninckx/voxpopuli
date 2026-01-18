#include "template.h"

//#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

inline float intersect_cube( Ray& ray )
{
	// branchless slab method by Tavian
	const float tx1 = -ray.O.x * ray.rD.x, tx2 = (1 - ray.O.x) * ray.rD.x;
	float ty, tz, tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
	const float ty1 = -ray.O.y * ray.rD.y, ty2 = (1 - ray.O.y) * ray.rD.y;
	ty = min( ty1, ty2 ), tmin = max( tmin, ty ), tmax = min( tmax, max( ty1, ty2 ) );
	const float tz1 = -ray.O.z * ray.rD.z, tz2 = (1 - ray.O.z) * ray.rD.z;
	tz = min( tz1, tz2 ), tmin = max( tmin, tz ), tmax = min( tmax, max( tz1, tz2 ) );
	if (tmin == tz) ray.axis = 2; else if (tmin == ty) ray.axis = 1;
	return tmax >= tmin ? tmin : 1e34f;
}

inline bool point_in_cube( const float3& pos )
{
	// test if pos is inside the cube
	return pos.x >= 0 && pos.y >= 0 && pos.z >= 0 &&
		pos.x <= 1 && pos.y <= 1 && pos.z <= 1;
}

Scene::Scene()
{
    grid = (uint*)MALLOC64(GRIDSIZE3 * sizeof(uint));
    if (true) {
        //memset(grid, 0, GRIDSIZE3 * sizeof(uint));
        gzFile f = gzopen("assets/viking.bin", "rb");
        int3 size;
        gzread(f, &size, sizeof(int3));
        gzread(f, grid, size.x * size.y * size.z * 4);
        gzclose(f);
    }
    else {
        memset(grid, 0, GRIDSIZE3 * sizeof(uint));
        for (int i = 0; i < 128 * 128 * 128; i++)
        {
            const int x = i & 127, y = (i >> 7) & 127, z = i >> 14;
            bool v = y < 10, p = z > 90 && z < 94;
            p &= (x & 31) > 10 && (x & 31) < 18 && y < 80;
            Set(x, y, z, p || v ? 0xffffff : 0);
        }
    }

    skyPixels = stbi_loadf("assets/kloofendal_48d_partly_cloudy_puresky_4k.hdr", &skyWidth, &skyHeight, &skyBpp,0 );
    for (int i{ 0 }; i < skyWidth * skyHeight * 3; i++)
        skyPixels[i] = sqrtf(skyPixels[i]);

    return;



	// allocate room for the world
	// initialize the scene using Perlin noise, parallel over z
    #pragma omp parallel for schedule(dynamic)
    for (int z = 0; z < WORLDSIZE; z++)
    {
        const float fz = (float)z / WORLDSIZE;
        for (int y = 0; y < WORLDSIZE; y++)
        {
            float fx = 0, fy = (float)y / WORLDSIZE;
            for (int x = 0; x < WORLDSIZE; x++, fx += 1.0f / WORLDSIZE)
            {
                const float n = noise3D(fx, fy, fz);
                Set(x, y, z, n > 0.09f ? 0x020101 * y : 0);
            }
        }
    }
    //for (int z = 0; z < 128; z++)
    //    for (int y = 0; y < 128; y++) for (int x = 0; x < 128; x++)
    //        if (x < 2 || x > 125 || z > 125 || y < 2 || y > 125)
    //            //Set(x, y, z, 0xeeeeee);
    //            Set(x, y, z, y == 1 ? 0x19999bb : 0xffffff);
    //        else if (y > 30 && y < 50 && z > 50 && z < 70 && x > 20)
    //            if (x < 40) Set(x, y, z, 0x3ff7777);
    //            else if (x > 55 && x < 75) Set(x, y, z, 0x2aaffaa);
    //            else if (x > 90 && x < 110) Set(x, y, z, 0x7777ff);



    /*for (int z = 0; z < WORLDSIZE; z++) {
        for (int x = 0; x < WORLDSIZE; x++) {
            Set(x, 0, z, 0);
            Set(x, WORLDSIZE - 1, z, 0);
        }
    }
    for (int y = 0; y < WORLDSIZE; y++) {
        for (int z = 0; z < WORLDSIZE; z++) {
            Set(0, y, z, 0);
            Set(WORLDSIZE - 1,y, z, 0);
        }

        for (int x = 0; x < WORLDSIZE; x++) {
            Set(x,y,0,0);
            Set(x,y,WORLDSIZE - 1, 0);
        }
    }*/
}

void Scene::Set( const uint x, const uint y, const uint z, const uint v )
{
	grid[x + y * GRIDSIZE + z * GRIDSIZE2] = v;
}

bool Scene::Setup3DDDA( Ray& ray, DDAState& state ) const
{
	// if ray is not inside the world: advance until it is
	state.t = 0;
	bool startedInGrid = point_in_cube( ray.O );
	if (!startedInGrid)
	{
		state.t = intersect_cube( ray );
		if (state.t > 1e33f) return false; // ray misses voxel data entirely
	}
	// setup amanatides & woo - assume world is 1x1x1, from (0,0,0) to (1,1,1)
	static const float cellSize = 1.0f / GRIDSIZE;
	state.step = make_int3( 1 - ray.Dsign * 2 );
	const float3 posInGrid = GRIDSIZE * (ray.O + (state.t + 0.00005f) * ray.D);
	const float3 gridPlanes = (ceilf( posInGrid ) - ray.Dsign) * cellSize;
	const int3 P = clamp( make_int3( posInGrid ), 0, GRIDSIZE - 1 );
	state.X = P.x, state.Y = P.y, state.Z = P.z;
	state.tdelta = cellSize * float3( state.step ) * ray.rD;
	state.tmax = (gridPlanes - ray.O) * ray.rD;
	// detect rays that start inside a voxel
	uint cell = grid[P.x + P.y * GRIDSIZE + P.z * GRIDSIZE2];
	ray.inside = cell != 0 && startedInGrid;
	// proceed with traversal
	return true;
}

void Scene::FindNearest( Ray& ray ) const
{
	// nudge origin
	ray.O += EPSILON * ray.D;
	// setup Amanatides & Woo grid traversal
	DDAState s;
	if (!Setup3DDDA( ray, s )) return;
	uint cell, lastCell = 0, axis = ray.axis;
	if (ray.inside)
	{
		// start stepping until we find an empty voxel
		while (1)
		{
			cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
			if (!cell) break;
			lastCell = cell;
			if (s.tmax.x < s.tmax.y)
			{
				if (s.tmax.x < s.tmax.z) { s.t = s.tmax.x, s.X += s.step.x, axis = 0; if (s.X >= GRIDSIZE) break; s.tmax.x += s.tdelta.x; }
				else { s.t = s.tmax.z, s.Z += s.step.z, axis = 2; if (s.Z >= GRIDSIZE) break; s.tmax.z += s.tdelta.z; }
			}
			else
			{
				if (s.tmax.y < s.tmax.z) { s.t = s.tmax.y, s.Y += s.step.y, axis = 1; if (s.Y >= GRIDSIZE) break; s.tmax.y += s.tdelta.y; }
				else { s.t = s.tmax.z, s.Z += s.step.z, axis = 2; if (s.Z >= GRIDSIZE) break; s.tmax.z += s.tdelta.z; }
			}
		}
		ray.voxel = lastCell; // we store the voxel we just left
	}
	else
	{
		// start stepping until we find a filled voxel
		while (1)
		{
			cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
			if (cell) break; else if (s.tmax.x < s.tmax.y)
			{
				if (s.tmax.x < s.tmax.z) { s.t = s.tmax.x, s.X += s.step.x, axis = 0; if (s.X >= GRIDSIZE) break; s.tmax.x += s.tdelta.x; }
				else { s.t = s.tmax.z, s.Z += s.step.z, axis = 2; if (s.Z >= GRIDSIZE) break; s.tmax.z += s.tdelta.z; }
			}
			else
			{
				if (s.tmax.y < s.tmax.z) { s.t = s.tmax.y, s.Y += s.step.y, axis = 1; if (s.Y >= GRIDSIZE) break; s.tmax.y += s.tdelta.y; }
				else { s.t = s.tmax.z, s.Z += s.step.z, axis = 2; if (s.Z >= GRIDSIZE) break; s.tmax.z += s.tdelta.z; }
			}
		}
		ray.voxel = cell;
	}
	ray.t = s.t;
	ray.axis = axis;
}

bool Scene::IsOccluded( Ray& ray ) const
{
	// nudge origin
	ray.O += EPSILON * ray.D;
	ray.t -= EPSILON * 2.0f;
	// setup Amanatides & Woo grid traversal
	DDAState s;
	if (!Setup3DDDA( ray, s )) return false;
	// start stepping
	while (s.t < ray.t)
	{
		const uint cell = grid[s.X + s.Y * GRIDSIZE + s.Z * GRIDSIZE2];
		if (cell && cell >> 24 != 2) /* we hit a solid voxel */ return s.t < ray.t;
		if (s.tmax.x < s.tmax.y)
		{
			if (s.tmax.x < s.tmax.z) { if ((s.X += s.step.x) >= GRIDSIZE) return false; s.t = s.tmax.x, s.tmax.x += s.tdelta.x; }
			else { if ((s.Z += s.step.z) >= GRIDSIZE) return false; s.t = s.tmax.z, s.tmax.z += s.tdelta.z; }
		}
		else
		{
			if (s.tmax.y < s.tmax.z) { if ((s.Y += s.step.y) >= GRIDSIZE) return false; s.t = s.tmax.y, s.tmax.y += s.tdelta.y; }
			else { if ((s.Z += s.step.z) >= GRIDSIZE) return false; s.t = s.tmax.z, s.tmax.z += s.tdelta.z; }
		}
	}
	return false;
}