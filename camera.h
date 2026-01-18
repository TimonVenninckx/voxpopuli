#pragma once

// default screen resolution
#define SCRWIDTH	640
#define SCRHEIGHT	400
// #define FULLSCREEN
#define DOUBLESIZE

typedef float4 Plane;
struct Frustum { Plane plane[4]; };
inline float distance(Plane& plane, float3& pos) {
    return dot(float3(plane), pos) - plane.w;
}

namespace Tmpl8 {

class Camera
{
public:
	Camera();
	~Camera();
	Ray GetPrimaryRay( const float x, const float y );
	bool HandleInput( const float t );
	float aspect = (float)SCRWIDTH / (float)SCRHEIGHT;
	float3 camPos, camTarget;
	float3 topLeft, topRight, bottomLeft;

    Frustum BuildFrustum();
};

}