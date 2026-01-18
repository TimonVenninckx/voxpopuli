#include "template.h"

// -----------------------------------------------------------
// Calculate light transport via a ray
// -----------------------------------------------------------


struct Light {
    float3 pos;
    float3 color;
};


std::vector<Light> lights
{
    //Light { float3{  2,  2.5, 0.1}, float3{1.f, 0.f, 0.f} },
    //Light { float3{ 1.0, 0.0,-1.0}, float3{ 0.6f, .6f, 0.f} },
    Light { float3{ 0.0, 0.0, 0.0}, float3{ 1.f, 1.f,1.f} },
};


/*
refract()
according to glsl manual
k = 1 - eta * eta * (1 - dot(N, D) * dot(N, D));
if (k < 0)
    R = 0;
else
    R = eta * D - (eta * dot(N, D) + sqrtf(k)) * N;*/


float3 RandomPointOnLight() {
    return float3(RandomFloat() - 1, 3, RandomFloat() - 1);
}


float glassRefractiveIndex{ 1.46f };
float3 Renderer::Trace(Ray& ray, int depth, int, int)
{
    if (depth == 5)
        return 0;
    scene.FindNearest(ray);
    if (ray.voxel == 0) {
        // hit the sky
        uint u = scene.skyWidth * atan2f(ray.D.z, ray.D.x) * INV2PI - 0.5f;
        uint v = scene.skyHeight* acosf(ray.D.y) * INVPI - 0.5f;
        uint skyIdx = (u + v * scene.skyWidth) % (scene.skyWidth * scene.skyHeight);
        return 0.65f * float3(scene.skyPixels[skyIdx * 3], scene.skyPixels[skyIdx * 3 + 1], scene.skyPixels[skyIdx * 3 + 2]);
    }
    float3 I = ray.IntersectionPoint();
    float3 N = ray.GetNormal();
    float3 albedo = ray.GetAlbedo();
    if (ray.voxel >> 24 == 1) {
        float n1 = 1.f;  // air 
        float n2 = glassRefractiveIndex; // glass
        if (ray.inside)
            std::swap(n1, n2);

        float d = dot(N, -ray.D);
        float eta = n1 / n2;
        float k = 1 - eta * eta * (1 - pow2f(d));
        float3 reflected{};
        float3 refracted{};
        if (k > 0) {
            float3 R = eta * ray.D + (eta * d - sqrtf(k))* N;
            refracted = Trace(Ray(I, R), depth + 1);
        }
        float3 R = ray.D - 2 * N * dot(ray.D,N);
        reflected = Trace(Ray( I,R ), depth + 1);

        float R0 = pow2f((n1 - n2) / (n1 + n2));
        float fresnel = k > 0 ? R0 + (1 - R0) * pow5f(1 - d) : 1;

        if (RandomFloat() < fresnel) {
            return reflected;
        }
        else {
            return refracted * albedo;
        }
        /*return //lerp
            fresnel * reflected + 
            (1 - fresnel) * refracted * albedo;*/
    }
    float3 result{};
    for (auto& light : lights) {
        float3 L = (light.pos +  RandomPointOnLight()) - I;
        float distance = length(L);
        L = normalize(L);
        float cosa = max(0.0f, dot(N, L));
        Ray shadowRay(I, L, distance);
        if (scene.IsOccluded(shadowRay)){
            //result += float3(0.01f) * albedo; 
            continue;
        }

        result += 20 * albedo * light.color * cosa / pow2f(distance);
    }
    return result;
}

//float3 Renderer::Trace( Ray& ray, int depth, int, int /* we'll use these later */ )
//{
//    if (depth == 10)
//        return{ 0 };
//
//    const float3 skyColor{ .7f,.8f,1.f };
//
//	scene.FindNearest( ray );
//	if (ray.voxel == 0) return skyColor; // or a fancy sky color
//	float3 N = ray.GetNormal();
//	float3 I = ray.IntersectionPoint();
//	float3 albedo = ray.GetAlbedo();
//
//    if (ray.voxel >> 24 == 1) {
//        // reflective
//        float3 R = ray.D - 2 * N * dot(ray.D, N);
//        Ray reflection{ I,R };
//        return albedo * Trace(reflection, depth + 1);
//    }
//    else if (ray.voxel >> 24 == 2) {
//        //glass 
//        // refractive indexes
//        float n1 = 1.f;  // air 
//        float n2 = glassRefractiveIndex; // glass
//        if (ray.inside)
//            std::swap(n1, n2);
//        float d = dot(N, -ray.D);
//        float eta = n1 / n2;
//        float k = 1 - eta * eta * (1 - pow2f(d));
//
//        float3 reflected{ 0 };
//        float3 refracted{ 0 };
//        if (k > 0) {
//            float3 R = eta * ray.D + (eta * d - sqrtf(k)) * N;
//            refracted = Trace(Ray(I, R), depth + 1);
//        }
//        else {
//            float3 R = ray.D - 2 * N * d;
//            reflected = Trace(Ray(I, R), depth + 1);
//        }
//        float R0 = pow2f((n1 - n2) / (n1 + n2));
//        float fresnel = k > 0 ? R0 + (1 - R0) * pow5f(1 - d) : 1;
//        return albedo * (fresnel * reflected + (1 - fresnel) * refracted);
//    }
//    else if (ray.voxel >> 24 == 3) {
//        // reflective with shadow
//        float3 R = ray.D - 2 * N * dot(ray.D, N);
//        Ray reflection{ I,R };
//
//        return 0.9f * Trace(reflection, depth + 1) + 0.1f* CalculateDiffuse(N,I,albedo);
//    }
//    else {
//        // difuse
//        return CalculateDiffuse(N,I,albedo);
//    }
//}
                     // normal, intersection, albedo color
float3 Renderer::CalculateDiffuse(float3 N, float3 I, float3 albedo){
    float3 result{};
    for (auto& light : lights) {
        float3 L = light.pos - I;

        // diffuse material
        float distance = length(L);
        float cosa = max(0.0f, dot(N, normalize(L)));
        Ray colorRay{ I,L,distance };
        
        if (scene.IsOccluded(colorRay)) continue;

        result += albedo * light.color * (1.f / pow2f(distance)) * cosa;
    }
    if (result.x == 0.f && result.y == 0.f && result.z == 0.f)
        return albedo * 0.1f;
    return result;
}

// -----------------------------------------------------------
// Application initialization - Executed once, at app start
// -----------------------------------------------------------
void Renderer::Init()
{
    accumulator = new float3[SCRWIDTH * SCRHEIGHT];
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
}

// -----------------------------------------------------------
// Main application tick function - Executed every frame
// -----------------------------------------------------------
bool rotatingLight = false;
static int spp = 1;
void Renderer::Tick( float deltaTime )
{
    static Timer g;
    if (rotatingLight) {
        lights[0].pos = float3(
            0.1f,
            0.f + 3.f * cos(0.6f * g.elapsed()),
            0.f + 3.f * sin(0.6f * g.elapsed())
        );
    }
    //std::cout << "camera position:\t" << camera.camPos.x << ' ' << camera.camPos.y << ' ' << camera.camPos.z << "\n";

	// high-resolution timer, see template.h
	Timer t;
	// pixel loop: lines are executed as OpenMP parallel tasks (disabled in DEBUG)
    const float scale = 1.f / spp++;
#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCRHEIGHT; y++)
	{
		// trace a primary ray for each pixel on the line
		for (int x = 0; x < SCRWIDTH; x++)
		{
			Ray r = camera.GetPrimaryRay( (float)x + RandomFloat(), (float)y + RandomFloat());
			accumulator [x + y * SCRWIDTH] += Trace( r );
            float3 average = accumulator[x + y * SCRWIDTH] * scale;
            //float3 c = clamp(average, 0.f, 1.f);
            //c = pow2f(c, 1.0f / 2.2f);
			screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(average);
		}
	}
	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf( "%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000 );
	// handle user input
    if (camera.HandleInput(deltaTime)) {
        spp = 1;
        memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
    }
}

// -----------------------------------------------------------
// Update user interface (imgui)
// -----------------------------------------------------------
void Renderer::UI()
{
	// ray query on mouse
	Ray r = camera.GetPrimaryRay( (float)mousePos.x, (float)mousePos.y );
	scene.FindNearest( r );
	ImGui::Text( "voxel: %i", r.voxel );

    ImGui::SliderFloat("glass refraction index", &glassRefractiveIndex, 0.1f, 20.f);
    ImGui::Checkbox("rotating light", &rotatingLight);
}