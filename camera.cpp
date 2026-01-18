#include "template.h"

Camera::Camera()
{
	// try to load a camera
	FILE* f = fopen( "camera.bin", "rb" );
	if (f)
	{
		fread( this, 1, sizeof( Camera ), f );
		fclose( f );
	}
	else
	{
		// setup a basic view frustum
		camPos = float3( 0, .2, -2 );
		camTarget = float3( 0, .2,-1 );
		topLeft = float3( -aspect, 1, 0 );
		topRight = float3( aspect, 1, 0 );
		bottomLeft = float3( -aspect, -1, 0 );
	}
}

Camera::~Camera()
{
	// save current camera
	FILE* f = fopen( "camera.bin", "wb" );
	fwrite( this, 1, sizeof( Camera ), f );
	fclose( f );
}

Ray Camera::GetPrimaryRay( const float x, const float y )
{
	// calculate pixel position on virtual screen plane
	const float u = (float)x * (1.0f / SCRWIDTH);
	const float v = (float)y * (1.0f / SCRHEIGHT);
	const float3 P = topLeft + u * (topRight - topLeft) + v * (bottomLeft - topLeft);
	// return Ray( camPos, normalize( P - camPos ) );
	return Ray( camPos, P - camPos );
	// Note: no need to normalize primary rays in a pure voxel world
	// TODO: 
	// - if we have other primitives as well, we *do* need to normalize!
	// - there are far cooler camera models, e.g. try 'Panini projection'.
}

Frustum Tmpl8::Camera::BuildFrustum()
{
    Frustum f;
    f.plane[0] = cross(topLeft - bottomLeft, topLeft - camPos);
    f.plane[1] = cross(topRight - camPos, topLeft - bottomLeft);
    f.plane[2] = cross(topRight - topLeft, topLeft - camPos);
    f.plane[3] = cross(bottomLeft - camPos, topRight - topLeft);
    for (int i{ 0 }; i < 4; i++) {
        f.plane[i].w = distance(f.plane[i], camPos);
    }
    return f;
}

//bool Camera::HandleInput(const float t)
//{
//    static float r = 0, a;
//    r += t * 0.001f, a = 0.7f * sinf(r) + 1.6f;
//    camPos = 0.7f * float3(cosf(a) + 0.7, 0.7f, sinf(a) + 0.7f);
//    mat4 M = mat4::LookAt(camPos, float3(5, 3, 5) / 10);
//    float3 x(M[0], M[4], M[8]), y(M[1], M[5], M[9]), z(M[2], M[6], M[10]);
//    topLeft = camPos + 2 * z - aspect * x + y;
//    topRight = camPos + 2 * z + aspect * x + y;
//    bottomLeft = camPos + 2 * z - aspect * x - y;
//    return false;
//}



bool Camera::HandleInput( const float t )
{
	if (!WindowHasFocus()) return false;
	float speed = 0.0015f * t;
	float3 ahead = normalize( camTarget - camPos );
	float3 tmpUp( 0, 1, 0 );
	float3 right = normalize( cross( tmpUp, ahead ) );
	float3 up = normalize( cross( ahead, right ) );
	bool changed = false;
    if (true) {
        if (IsKeyDown(GLFW_KEY_UP)) camTarget += speed * up, changed = true;
        if (IsKeyDown(GLFW_KEY_DOWN)) camTarget -= speed * up, changed = true;
    }
    else {
        if (IsKeyDown(GLFW_KEY_UP)) camTarget -= speed * up, changed = true;
        if (IsKeyDown(GLFW_KEY_DOWN)) camTarget += speed * up, changed = true;
    }
	if (IsKeyDown( GLFW_KEY_LEFT )) camTarget -= speed * right, changed = true;
	if (IsKeyDown( GLFW_KEY_RIGHT )) camTarget += speed * right, changed = true;
	ahead = normalize( camTarget - camPos );
	right = normalize( cross( tmpUp, ahead ) );
	up = normalize( cross( ahead, right ) );
	if (IsKeyDown( GLFW_KEY_A )) camPos -= speed * right, changed = true;
	if (IsKeyDown( GLFW_KEY_D )) camPos += speed * right, changed = true;
	if (GetAsyncKeyState( 'W' )) camPos += speed * ahead, changed = true;
	if (IsKeyDown( GLFW_KEY_S )) camPos -= speed * ahead, changed = true;
	if (IsKeyDown( GLFW_KEY_R )) camPos += speed * up, changed = true;
	if (IsKeyDown( GLFW_KEY_F )) camPos -= speed * up, changed = true;
	camTarget = camPos + ahead;
	ahead = normalize( camTarget - camPos );
	up = normalize( cross( ahead, right ) );
	right = normalize( cross( up, ahead ) );
	topLeft = camPos + 2 * ahead - aspect * right + up;
	topRight = camPos + 2 * ahead + aspect * right + up;
	bottomLeft = camPos + 2 * ahead - aspect * right - up;
	if (!changed) return false;
	return true;
}