// -*- mode: c++ -*-

//**************
// Visualization
//**************

struct ray{
  float4 origin;
  float4 direction;
};

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

float getDensityFromVolume(const float4 p, const int resolution, __global float* volumeData){
    // Convert [0,1] coordinates to [0, resolution-1] indices
    float x = p.x * (resolution - 1);
    float y = p.y * (resolution - 1);
    float z = p.z * (resolution - 1);

    // Clamp to valid range
    x = clamp(x, 0.0f, (float)(resolution - 1 - 1e-4));
    y = clamp(y, 0.0f, (float)(resolution - 1 - 1e-4));
    z = clamp(z, 0.0f, (float)(resolution - 1 - 1e-4));

    // Integer coordinates
    int x0 = (int)floor(x);
    int y0 = (int)floor(y);
    int z0 = (int)floor(z);
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;

    // Fractional parts
    float dx = x - x0;
    float dy = y - y0;
    float dz = z - z0;

    // Trilinear interpolation
    float c000 = volumeData[x0 + y0 * resolution + z0 * resolution * resolution];
    float c001 = volumeData[x0 + y0 * resolution + z1 * resolution * resolution];
    float c010 = volumeData[x0 + y1 * resolution + z0 * resolution * resolution];
    float c011 = volumeData[x0 + y1 * resolution + z1 * resolution * resolution];
    float c100 = volumeData[x1 + y0 * resolution + z0 * resolution * resolution];
    float c101 = volumeData[x1 + y0 * resolution + z1 * resolution * resolution];
    float c110 = volumeData[x1 + y1 * resolution + z0 * resolution * resolution];
    float c111 = volumeData[x1 + y1 * resolution + z1 * resolution * resolution];

    return 
        (1-dx)*(1-dy)*(1-dz)*c000 +
        dx*(1-dy)*(1-dz)*c100 +
        (1-dx)*dy*(1-dz)*c010 +
        dx*dy*(1-dz)*c110 +
        (1-dx)*(1-dy)*dz*c001 +
        dx*(1-dy)*dz*c101 +
        (1-dx)*dy*dz*c011 +
        dx*dy*dz*c111;
}
float4 getNormalFromVolume(const float4 p, const int resolution, __global float* volumeData){
  float4 normal;

  normal.x = getDensityFromVolume((float4)(p.x + 2.0f/resolution, p.y, p.z, 0.0f), resolution, volumeData) -
    getDensityFromVolume((float4)(p.x - 2.0f/resolution, p.y, p.z, 0.0f), resolution, volumeData);
  normal.y = getDensityFromVolume((float4)(p.x, p.y + 2.0f/resolution, p.z, 0.0f), resolution, volumeData) -
    getDensityFromVolume((float4)(p.x, p.y - 2.0f/resolution, p.z, 0.0f), resolution, volumeData);
  normal.z = getDensityFromVolume((float4)(p.x, p.y, p.z + 2.0f/resolution, 0.0f), resolution, volumeData) -
    getDensityFromVolume((float4)(p.x, p.y, p.z - 2.0f/resolution, 0.0f), resolution, volumeData);
  normal.w = 0.0f;

  if(dot(normal, normal) < 0.001f){
    normal = (float4)(0.0f, 0.0f, 1.0f, 0.0f);
  }

  return normalize(normal);
}

// Iso-sourface raycasting
__kernel
void isosurface(const int width, const int height, __global float4* visualizationBuffer,
                   const int resolution, __global float* volumeData,
		   const float isoValue, const float16 invViewMatrix){
  int2 id = (int2)(get_global_id(0), get_global_id(1));

  float2 uv = (float2)( (id.x / (float) width)*2.0f-1.0f, (id.y / (float) height)*2.0f-1.0f );

  float4 boxMin = (float4)(-1.0f, -1.0f, -1.0f,1.0f);
  float4 boxMax = (float4)(1.0f, 1.0f, 1.0f,1.0f);

  // calculate eye ray in world space
  struct ray eyeRay;

  eyeRay.origin = (float4)(invViewMatrix.sC, invViewMatrix.sD, invViewMatrix.sE, invViewMatrix.sF);

  float4 temp = normalize(((float4)(uv.x, uv.y, -2.0f, 0.0f)));
  eyeRay.direction.x = dot(temp, ((float4)(invViewMatrix.s0,invViewMatrix.s1,invViewMatrix.s2,invViewMatrix.s3)));
  eyeRay.direction.y = dot(temp, ((float4)(invViewMatrix.s4,invViewMatrix.s5,invViewMatrix.s6,invViewMatrix.s7)));
  eyeRay.direction.z = dot(temp, ((float4)(invViewMatrix.s8,invViewMatrix.s9,invViewMatrix.sA,invViewMatrix.sB)));
  eyeRay.direction.w = 0.0f;

  float4 color = (float4)(0.0f);

  float tnear, tfar;
  int hit = intersectBox(eyeRay.origin, eyeRay.direction, boxMin, boxMax, &tnear, &tfar);

  if(hit){
      if(tnear < 0.0f) tnear = 0.0f;
      float maxStep = 256.0f;
      float step = (tfar - tnear)/maxStep;
      float t = tnear;
      
      for(int i=0; i<maxStep; i++){
          float4 pos = (eyeRay.origin + t * eyeRay.direction + 1.0f)/2.0f; // Map to [0,1]^3
          float density = getDensityFromVolume(pos, resolution, volumeData);
          
          if(density > isoValue){
            
            float4 normal = getNormalFromVolume(pos, resolution, volumeData);
            float3 lightDir = normalize((float3)(0.3f, -1.0f, 0.2f)); // Example light direction
            float diffuse = clamp(dot(normal.xyz, lightDir), 0.0f, 1.0f);
            
            float hemisphere = 0.5f + 0.5f * dot(normal.xyz, lightDir);
            color = (float4)(hemisphere, hemisphere, hemisphere, 1.0f);
            break;
          }
          t += step;
      }
  }

  if(id.x < width && id.y < height){
    visualizationBuffer[id.x + id.y * width] = color;
  }
}

// Alpha blended
__kernel
void alphaBlended(const int width, const int height, __global float4* visualizationBuffer,
		  const int resolution, __global float* volumeData,
		  const float alphaExponent, const float alphaCenter,
		  const float16 invViewMatrix){
  int2 id = (int2)(get_global_id(0), get_global_id(1));

  float2 uv = (float2)( (id.x / (float) width)*2.0f-1.0f, (id.y / (float) height)*2.0f-1.0f );

  float4 boxMin = (float4)(-1.0f, -1.0f, -1.0f,1.0f);
  float4 boxMax = (float4)(1.0f, 1.0f, 1.0f,1.0f);

  // calculate eye ray in world space
  struct ray eyeRay;

  eyeRay.origin = (float4)(invViewMatrix.sC, invViewMatrix.sD, invViewMatrix.sE, invViewMatrix.sF);

  float4 temp = normalize(((float4)(uv.x, uv.y, -2.0f, 0.0f)));
  eyeRay.direction.x = dot(temp, ((float4)(invViewMatrix.s0,invViewMatrix.s1,invViewMatrix.s2,invViewMatrix.s3)));
  eyeRay.direction.y = dot(temp, ((float4)(invViewMatrix.s4,invViewMatrix.s5,invViewMatrix.s6,invViewMatrix.s7)));
  eyeRay.direction.z = dot(temp, ((float4)(invViewMatrix.s8,invViewMatrix.s9,invViewMatrix.sA,invViewMatrix.sB)));
  eyeRay.direction.w = 0.0f;
  

  float4 sum = (float4)(0.0f);

  float tnear, tfar;
  int hit = intersectBox(eyeRay.origin, eyeRay.direction, boxMin, boxMax, &tnear, &tfar);
if(hit){
    if(tnear < 0.0f) tnear = 0.0f;
    float maxStep = 256.0f;
    float step = (tfar - tnear) / maxStep;
    float t = tfar - 0.0001f;  // Start from back
    
    // ======== RAY MARCHING LOOP ========
    for(int i=0; i < maxStep; ++i){
        float4 pos = ((eyeRay.origin + t * eyeRay.direction) + 1.0f) / 2.0f;
        float density = getDensityFromVolume(pos, resolution, volumeData);
        
        // ********** TASK A: X-RAY **********
        // float alpha = pow(density, alphaExponent) * step;
        // float4 color = (float4)(1.0f);
        
        // ********** TASK B: SHADED **********
        float alpha = clamp(alphaExponent * (density - alphaCenter) + 0.5f, 0.0f, 1.0f);
        float4 normal = getNormalFromVolume(pos, resolution, volumeData);
        float3 lightDir = normalize((float3)(0.3f, -2.0f, 0.0f));
        float diffuse = 0.5f + 0.5f * dot(normal.xyz, lightDir);
        float4 color = (float4)(diffuse, diffuse, diffuse, 1.0f);
        
        // ********** TASK C: TRANSFER FUNCTION **********
        color *= (float4)(density, density*density, density*density*density, 1.0f) + 0.1f;

        sum = (1.0f - alpha) * sum + alpha * color;
        
        t -= step;
        if(t < tnear) break;
    }
}

  if(id.x < width && id.y < height){
    visualizationBuffer[id.x + id.y * width] = (float4)(sum);
  }

  
}

