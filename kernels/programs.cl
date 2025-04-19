/* -*- mode: c++ -*- */

#ifndef M_PI
#define M_PI 3.14159265f  // Added 'f' for float and conditional define
#endif

// Random generator
uint stepLCG(uint *z, uint A, uint C){
  return (*z) = (A * (*z) + C);
}

uint stepCTG(uint *z, uint S1, uint S2, uint S3, uint M){
  uint b=((((*z)<<S1)^(*z))>>S2);
  return (*z) = ((((*z)&M)<<S3)^b);
}

float getRandom(uint* rng1, uint* rng2, uint* rng3, uint* rng4){
  return 2.3283064365387e-10f * (
                                stepCTG(rng1, 13, 19, 12, 4294967294UL) ^
                                stepCTG(rng2, 2, 25, 4, 4294967288UL) ^
                                stepCTG(rng3, 3, 11, 17, 4294967280UL) ^
				stepLCG(rng4,1664525,1013904223UL)
                                );
}

// Scattering simulation
struct photon{
  float4 origin;
  float4 direction;
  float3 energy;   // changed to float3 for RGB
};

__kernel
void resetSimulation(const int resolution, __global float* simulationBuffer){
  int4 id = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
  int channel = get_global_id(3);
  simulationBuffer[channel*resolution*resolution*resolution + 
                  id.x + id.y*resolution + id.z*resolution*resolution] = 0.0f;
}

float4 getRandomDirection(uint* rng1, uint* rng2, uint* rng3, uint* rng4){
  float x, y, z;
  bool inside = false;
  while(!inside){
    x = getRandom(rng1, rng2, rng3, rng4) * 2.0f - 1.0f;
    y = getRandom(rng1, rng2, rng3, rng4) * 2.0f - 1.0f;
    z = getRandom(rng1, rng2, rng3, rng4) * 2.0f - 1.0f;
    if( (x*x + y*y + z*z) <= 1.0f){
      inside = true;
    }
  }
  if( x*x + y*y + z*z == 0.0f){
    x = 0.0f;
    y = 1.0f;
    z = 0.0f;
  }
  float vlen = sqrt(x*x + y*y + z*z);
  return (float4)(x/vlen, y/vlen, z/vlen, 0);
}

__constant float eps = 0.000001f;

// Change albedo to float3
__constant float3 albedo = (float3)(0.9f, 0.8f, 0.7f); // RGB albedos

__constant float densityScale = 1.0f;

float3 getDensity(float4 p){

   float base = 0.5f * densityScale;
   
   // ===== NEW OBJECTS =====
    // Add sphere at center
    float3 center = (float3)(0.5f, 0.5f, 0.5f);
    float radius = 0.15f;
    if(length(p.xyz - center) < radius) return (float3) (100.0f);

    // Add cube
    float3 cubeMin = (float3)(0.7f, 0.1f, 0.1f);
    float3 cubeMax = (float3)(0.9f, 0.3f, 0.3f);
    if(all(p.xyz > cubeMin) && all(p.xyz < cubeMax)) return (float3)(2.0f);

  // oszlopok
  for(float ix = 0.3f; ix < 0.9f; ix += 0.4f){
    for(float iy = 0.3f; iy < 0.9f; iy += 0.4f){
      float px2 = (p.x - ix) * (p.x - ix);
      float py2 = (p.y - iy) * (p.y - iy);
      if( px2 + py2 < 0.001f){
        return (float3)(100.0f);
      }
    }
  }

  for(float ix = 0.3f; ix < 0.9f; ix += 0.4f){
    for(float iz = 0.3f; iz < 0.9f; iz += 0.4f){
      float px2 = (p.x - ix) * (p.x - ix);
      float pz2 = (p.z - iz) * (p.z - iz);
      if( px2 + pz2 < 0.001f){
        return (float3)(100.0f);
      }
    }
  }

  for(float iz = 0.3f; iz < 0.9f; iz += 0.4f){
    for(float iy = 0.3f; iy < 0.9f; iy += 0.4f){
      float pz2 = (p.z - iz) * (p.z - iz);
      float py2 = (p.y - iy) * (p.y - iy);
      if( pz2 + py2 < 0.001f){
        return (float3)(100.0f);
      }
    }
  }

  // teto lyukkal
  if(p.y > 0.78f && p.y < 0.83f &&
     ( (p.x - 0.5f) * (p.x - 0.5f) + (p.z - 0.5f) * (p.z - 0.5f) ) > 0.001f)
    return (float3)(100.0f);


  // falak
  if(p.x < 0.02f) return (float3)(100.0f);
  if(p.y < 0.02f) return (float3)(100.0f);
  if(p.z > 0.98f) return (float3)(100.0f);

  // alul besurusodik
  if(p.y < 0.2f) return (float3)((1.0f - p.y) * 5.0f);

  return (float3)(0.5f * densityScale);
}

void tracePhotonRM(__global struct photon* p, float rnd){

  float s = -log(rnd) / densityScale;

  float t = 0.0f;
  float dt = 1.0f / 256.0f;
  float sum = 0.0f;
  float3 sigmat = 0.0f;

  while(sum < s){
    float4 samplePos = p->origin + t * p->direction;
    if(samplePos.x < 0.0f || samplePos.x > 1.0f ||
       samplePos.y < 0.0f || samplePos.y > 1.0f ||
       samplePos.z < 0.0f || samplePos.z > 1.0f){
      p->energy = 0.0f;
      break;
    } else {
      sigmat = getDensity(samplePos);
      sum += (sigmat.x + sigmat.y + sigmat.z)/3.0f * dt;  // Average density
      t += dt;
    }
  }

  p->origin = p->origin + p->direction * t;
  p->direction = p->direction;
  p->energy = p->energy * albedo;
}

void storePhoton(__global struct photon* p, const int resolution, __global float* simulationBuffer){
   if(p->energy.x < 0.01f && p->energy.y < 0.01f && p->energy.z < 0.01f) return;
    
    int x = p->origin.x * resolution;
    int y = p->origin.y * resolution;
    int z = p->origin.z * resolution;
    
    if(x < 0 || x >= resolution || y < 0 || y >= resolution || z < 0 || z >= resolution) return;
    
    // Store RGB channels in separate buffers
    simulationBuffer[x + y*resolution + z*resolution*resolution] += p->energy.x; // Red
    simulationBuffer[resolution*resolution*resolution + x + y*resolution + z*resolution*resolution] += p->energy.y; // Green
    simulationBuffer[2*resolution*resolution*resolution + x + y*resolution + z*resolution*resolution] += p->energy.z; // Blue
}

__kernel
void simulation(const int iteration, __global uint4* seed, __global struct photon* photons,
                const int resolution, __global float* simulationBuffer,
				const float4 lightSourcePosition){
  int id = get_global_id(0);

  // random generator setup
  uint rng1 = seed[id].s0;
  uint rng2 = seed[id].s1;
  uint rng3 = seed[id].s2;
  uint rng4 = seed[id].s3;


  // scatter simulation
   if(0 == iteration || dot(photons[id].energy, 1.0f) < 0.2f){
    // ===== AREA LIGHT =====
    float rx = getRandom(&rng1, &rng2, &rng3, &rng4) * 0.1f - 0.05f;
    float rz = getRandom(&rng1, &rng2, &rng3, &rng4) * 0.1f - 0.05f;
    photons[id].origin = (float4)(0.5f + rx, 0.99f, 0.5f + rz, 0.0f); // From 1.0f
    
    // Hemisphere sampling
    photons[id].direction = getRandomDirection(&rng1, &rng2, &rng3, &rng4);
   photons[id].direction.y = -fabs(photons[id].direction.y); // Force downward
    
    photons[id].energy = 50.0f;
  } else {
    photons[id].direction = getRandomDirection(&rng1, &rng2, &rng3, &rng4);
  }

  tracePhotonRM(&photons[id], getRandom(&rng1, &rng2, &rng3, &rng4));
  storePhoton(&photons[id], resolution, simulationBuffer);

  // random state store
  seed[id].s0 = rng1;
  seed[id].s1 = rng2;
  seed[id].s2 = rng3;
  seed[id].s3 = rng4;
}

//**************
// Visualization
//**************

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

int intersectBox(float4 r_o, float4 r_d, float4 boxmin, float4 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float4 invR = (float4)(1.0f,1.0f,1.0f,1.0f) / r_d;
    float4 tbot = invR * (boxmin - r_o);
    float4 ttop = invR * (boxmax - r_o);

    // re-order intersections to find smallest and largest on each axis
    float4 tmin = min(ttop, tbot);
    float4 tmax = max(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = max(max(tmin.x, tmin.y), max(tmin.x, tmin.z));
    float smallest_tmax = min(min(tmax.x, tmax.y), min(tmax.x, tmax.z));

        *tnear = largest_tmin;
        *tfar = smallest_tmax;

        return smallest_tmax > largest_tmin;
}

float getIllumination(const float4 p, const int resolution, __global float* simulationBuffer){
  int x = p.x * resolution;
  int y = p.y * resolution;
  int z = p.z * resolution;

  if(x > resolution -1 || x < 0) return(0.0f);
  if(y > resolution -1 || y < 0) return(0.0f);
  if(z > resolution -1 || z < 0) return(0.0f);

  return simulationBuffer[x + y * resolution + z * resolution * resolution];
}

__kernel
void visualization(const int width, const int height, __global float4* visualizationBuffer,
                   const int resolution, __global float* simulationBuffer,
				   const int iteration, const float16 invViewMatrix){
  int2 id = (int2)(get_global_id(0), get_global_id(1));

  float2 uv = (float2)( (id.x / (float) width)*2.0f-1.0f, (id.y / (float) height)*2.0f-1.0f );

  float4 boxMin = (float4)(-1.0f, -1.0f, -1.0f,1.0f);
  float4 boxMax = (float4)(1.0f, 1.0f, 1.0f,1.0f);

  // calculate eye ray in world space
  struct photon eyeRay;

  eyeRay.origin = (float4)(invViewMatrix.sC, invViewMatrix.sD, invViewMatrix.sE, invViewMatrix.sF);

  float4 temp = normalize(((float4)(uv.x, uv.y, -2.0f, 0.0f)));
  eyeRay.direction.x = dot(temp, ((float4)(invViewMatrix.s0,invViewMatrix.s1,invViewMatrix.s2,invViewMatrix.s3)));
  eyeRay.direction.y = dot(temp, ((float4)(invViewMatrix.s4,invViewMatrix.s5,invViewMatrix.s6,invViewMatrix.s7)));
  eyeRay.direction.z = dot(temp, ((float4)(invViewMatrix.s8,invViewMatrix.s9,invViewMatrix.sA,invViewMatrix.sB)));
  eyeRay.direction.w = 0.0f;

  float3 sum = 0.0f;
  float transparency = 1.0f;

  float tnear, tfar;
  int hit = intersectBox(eyeRay.origin, eyeRay.direction, boxMin, boxMax, &tnear, &tfar);
if(hit){
    float maxStep = 256.0f;
    float step = (tfar - tnear)/maxStep;
    float t = tfar - 0.0001f;
    
    for(int i=0; i<maxStep; ++i){
      float4 pos = ((eyeRay.origin + t*eyeRay.direction) + 1.0f)/2.0f;
      
      // Getting RGB channels
      float r = getIllumination(pos, resolution, simulationBuffer);
      float g = getIllumination(pos, resolution, &simulationBuffer[resolution*resolution*resolution]);
      float b = getIllumination(pos, resolution, &simulationBuffer[2*resolution*resolution*resolution]);
      
float3 illumination = (float3)(r, g, b) / (max(iteration, 1) * 50.0f) / (4.0f * M_PI) / pow((float)resolution, 3.0f) * step;

      float3 density = getDensity(pos);

      sum = (1.0f - density*step) * sum + illumination;
      transparency *= (1.0f - (density.x + density.y + density.z)/3.0f * step);
      t -= step;
    }
  }

  if(id.x < width && id.y < height){
    
 float3 sky = (float3)(0.17f, 0.4f, 0.6f) * transparency;

	float3 tmpColor = sky + sum;

	float gamma = 0.8f;
	tmpColor.x = pow(tmpColor.x, gamma);
    tmpColor.y = pow(tmpColor.y, gamma);
    tmpColor.z = pow(tmpColor.z, gamma);

	visualizationBuffer[id.x + id.y * width] = (float4)(tmpColor, 1.0f);
  }
}

