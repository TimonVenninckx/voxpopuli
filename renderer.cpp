#include "template.h"

// -----------------------------------------------------------
// Calculate light transport via a ray
// -----------------------------------------------------------


struct Light {
    float3 pos;
    float3 color;
};

constexpr float FSTR = 500.f;
std::vector<Light> lights
{
    //Light { float3{ 3.0,  1.5, 0.1}, float3{ 0.f, FSTR, 0.f} },
    //Light { float3{ 0.0,  1.5,-3.0}, float3{ 0.f,  0.f,FSTR} },
    //Light { float3{ 0.0,  1.5, 3.0}, float3{ FSTR, 0.f, 0.f} },
    Light { float3{ 0.0,  4.5,16.0}, float3{ FSTR, FSTR, FSTR} },
};



// imgui stuff
bool rotatingLight = false;
bool toneMapping = false;
bool blueNoise = true;

Frustum previous;
float3 prevCamPos;


float3 RandomPointOnLight(float3 dirTowardsLight,float r0, float r1) {
    // create plane from normal
    // build orthonormal basis(u,v) perpendicular to normal
    float3 normal = normalize(dirTowardsLight);
    float3 u = normalize(abs(normal.x) > .1f ? cross(float3(0, 1, 0), normal) : cross(float3(1, 0, 0), normal));
    float3 v = cross(normal, u);
    
    float r = sqrt(r0);
    float theta = 2 * PI * r1;
    return (u * cos(theta) + v * sin(theta)) * r;
    //return u * r0 + v * r1;
}


//float3 Renderer::Trace(Ray& ray, int depth, int x, int y)
//{
//    scene.FindNearest(ray);
//    if (ray.voxel == 0) return float3(0.4f, 0.5f, 1.0f);
//    float3 I = ray.IntersectionPoint();
//    float3 P((RandomFloat() + 1) / 2, (RandomFloat() + 6) / 5, 3);
//    float3 L = normalize(P - I);
//    if (scene.IsOccluded(Ray(I, L, 10))) return 0.2f;
//    return ray.GetAlbedo() * 2 * max(0.2f, dot(ray.GetNormal(), L));
//}

float glassRefractiveIndex{ 1.46f };
float3 Renderer::Trace(Ray& ray, int depth, int x, int y)
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
            refracted = Trace(Ray(I, R), depth + 1,x,y);
        }
        float3 R = ray.D - 2 * N * dot(ray.D,N);
        reflected = Trace(Ray( I,R ), depth + 1, x, y);

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
        float r0 = RandomFloat();
        float r1 = RandomFloat();
        static Surface blue("assets/LDR_RG01_0.png");
        if (blueNoise) {
            int pixel = blue.pixels[(x & 63) + (y & 63) * 64];
            int red = pixel >> 16, green = (pixel >> 8) & 255;
            r0 = fracf(red / 256.0f + frame * 0.61803399f);
            r1 = fracf(green / 256.0f + frame * 0.61803399f);
        }
        float3 lightOffset = RandomPointOnLight(light.pos - I,r0,r1);

        float3 L = (light.pos + lightOffset) - I;
        float distance = length(L);
        L = normalize(L);
        float cosa = max(0.0f, dot(N, L));
        Ray shadowRay(I, L, distance);
        if (scene.IsOccluded(shadowRay)){
            //result += float3(0.01f) * albedo; // ambient test 
            continue;
        }

        result += albedo * light.color * cosa / pow2f(distance);
    }
    return result;
}
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
    history = new float3[SCRWIDTH * SCRHEIGHT];
    
    memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
    memset(history,     0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
}

void Renderer::Shutdown() {
    delete accumulator;
    delete history;
}




// -----------------------------------------------------------
// Main application tick function - Executed every frame
// -----------------------------------------------------------
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
    frame++;
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
			Ray r = camera.GetPrimaryRay( (float)x, (float)y);
            float3 historySample{};
            float  historyWeight{};
            float3 sample = Trace(r, 0, x, y);
            if (r.voxel > 0) {
                float3 P = r.IntersectionPoint();
                float dLeft = distance(previous.plane[0], P);
                float dRight= distance(previous.plane[1], P);
                float dUp   = distance(previous.plane[2], P);
                float dDown = distance(previous.plane[3], P);

                float prev_x = (SCRWIDTH * dLeft) / (dLeft + dRight);
                float prev_y = (SCRHEIGHT * dUp) / (dUp + dDown);
            
                if (prev_x >= 0 && prev_x < SCRWIDTH - 1 && prev_y >= 0 && prev_y < SCRHEIGHT - 2) {
                    Ray prevR(prevCamPos, P - prevCamPos);
                    scene.FindNearest(prevR);
                    float3 Q = prevCamPos + prevR.t * prevR.D;
                    if (sqrLength(P - Q) < 0.001f) {
                        historyWeight = 0.95f;
                        int a = (int)prev_x + (int)prev_y * SCRWIDTH;

                        float fx = fracf(prev_x);
                        float fy = fracf(prev_y);
                        float w0 = (1 - fx) * (1 - fy);
                        float w1 = fx * (1 - fy);
                        float w2 = (1 - fx) * fy;
                        float w3 = fx * fy;
                        float3 p0 = history[a];
                        float3 p1 = history[a + 1];
                        float3 p2 = history[a + SCRWIDTH];
                        float3 p3 = history[a + SCRWIDTH + 1];

                        historySample = p0 * w0 + p1 * w1 + p2 * w2 + p3 * w3;
                    }
                }
            }
            float3 pixel = historyWeight * historySample + (1.0 - historyWeight) * sample;
            accumulator[x + y * SCRWIDTH] = pixel;
            //screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(Trace(r, 0, x, y));
		}
	}

#pragma omp parallel for schedule(dynamic)
    for (int y = 1; y < SCRHEIGHT - 1; y++) {
        for (int x = 1; x < SCRWIDTH - 1; x++) {
            float3 pixel{};
            //pixel += accumulator[x + -1 + (y - 1) * SCRWIDTH];
            pixel += accumulator[x      + (y - 1) * SCRWIDTH] *-.5f;
            //pixel += accumulator[x +  1 + (y - 1) * SCRWIDTH];

            pixel += accumulator[x + -1 + (y    ) * SCRWIDTH] *-.5f;
            pixel += accumulator[x      + (y    ) * SCRWIDTH] * 3.f;
            pixel += accumulator[x +  1 + (y    ) * SCRWIDTH] *-.5f;

            //pixel += accumulator[x + -1 + (y + 1) * SCRWIDTH];
            pixel += accumulator[x      + (y + 1) * SCRWIDTH] *-.5f;
            //pixel += accumulator[x +  1 + (y + 1) * SCRWIDTH];
            pixel = clamp(pixel, 0.f, 1.f);
            screen->pixels[x + y * SCRWIDTH] = RGBF32_to_RGB8(pixel);
        }
    }


	// performance report - running average - ms, MRays/s
	static float avg = 10, alpha = 1;
	avg = (1 - alpha) * avg + alpha * t.elapsed() * 1000;
	if (alpha > 0.05f) alpha *= 0.5f;
	float fps = 1000.0f / avg, rps = (SCRWIDTH * SCRHEIGHT) / avg;
	printf( "%5.2fms (%.1ffps) - %.1fMrays/s\n", avg, fps, rps / 1000 );
	// handle user input
    previous = camera.BuildFrustum();
    prevCamPos = camera.camPos;
    if (camera.HandleInput(deltaTime)) {
        //spp = 1;
        //memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
    }
    std::swap(history, accumulator);
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

    bool setValue{};
    if(ImGui::SliderFloat("glass refraction index", &glassRefractiveIndex, 0.5f, 2.f)) setValue = true;
    ImGui::Checkbox("rotating light", &rotatingLight);
    ImGui::Checkbox("Tone mapping", &toneMapping);
    if (ImGui::Checkbox("Blue Noise", &blueNoise)) setValue = true;
    
    if(setValue){
        spp = 1;
        memset(accumulator, 0, SCRWIDTH * SCRHEIGHT * sizeof(float3));
    }
}