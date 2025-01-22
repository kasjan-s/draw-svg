#include "texture.h"
#include "color.h"
#include "CS248/misc.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {


Sampler2D::~Sampler2D() {
    // Virtual destructor: no need to do anything for now.
}

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Advanced Task
  // Implement mipmap for trilinear filtering

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation

    if (level < 0 || level >= tex.mipmap.size()) {
      return Color(1,0,1,1);
    }
  
    // Get the appropriate mip level
    MipLevel& mip = tex.mipmap[level];

    // Ensure u and v are in the range [0, 1]
    u = clamp(u, 0.0f, 1.0f);
    v = clamp(v, 0.0f, 1.0f);

    // Calculate the coordinates for the texture
    int x = static_cast<int>(u * mip.width);
    int y = static_cast<int>(v * mip.height);


    // Calculate the texel index in the mip level's texels array
    int texel_index = y * mip.width + x;

    // Return the color at the nearest texel
    unsigned char r = mip.texels[texel_index * 4 + 0];  // Red channel
    unsigned char g = mip.texels[texel_index * 4 + 1];  // Green channel
    unsigned char b = mip.texels[texel_index * 4 + 2];  // Blue channel
    unsigned char a = mip.texels[texel_index * 4 + 3];  // Alpha channel

    // Return the color object (normalize the texel values to [0, 1] range)
    return Color(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 4: Implement bilinear filtering

// Ensure valid level

    if (level < 0 || level >= tex.mipmap.size()) {
      return Color(1,0,1,1);
    }
    // Fetch the mip level texture data
    const MipLevel& mip = tex.mipmap[level];

    // Clamp coordinates to [0, 1]
    u = clamp(u, 0.0f, 1.0f);
    v = clamp(v, 0.0f, 1.0f);

    // Scale u and v to the size of the texture
    float tex_u = u * (mip.width - 1);
    float tex_v = v * (mip.height - 1);

    // Compute integer coordinates of the four nearest texels
    int x0 = static_cast<int>(tex_u);
    int y0 = static_cast<int>(tex_v);
    int x1 = std::min(x0 + 1, static_cast<int>(mip.width - 1));
    int y1 = std::min(y0 + 1, static_cast<int>(mip.height - 1));

    // Compute the interpolation weights
    float dx = tex_u - x0;
    float dy = tex_v - y0;
    float wx0 = 1.0f - dx;
    float wx1 = dx;
    float wy0 = 1.0f - dy;
    float wy1 = dy;

    // Access the texels at the four corners of the bilinear interpolation
    const unsigned char* texel00 = &mip.texels[4 * (y0 * mip.width + x0)];
    const unsigned char* texel01 = &mip.texels[4 * (y0 * mip.width + x1)];
    const unsigned char* texel10 = &mip.texels[4 * (y1 * mip.width + x0)];
    const unsigned char* texel11 = &mip.texels[4 * (y1 * mip.width + x1)];

    // Extract RGBA values from texels
    Color c00(texel00[0] / 255.0f, texel00[1] / 255.0f, texel00[2] / 255.0f, texel00[3] / 255.0f);
    Color c01(texel01[0] / 255.0f, texel01[1] / 255.0f, texel01[2] / 255.0f, texel01[3] / 255.0f);
    Color c10(texel10[0] / 255.0f, texel10[1] / 255.0f, texel10[2] / 255.0f, texel10[3] / 255.0f);
    Color c11(texel11[0] / 255.0f, texel11[1] / 255.0f, texel11[2] / 255.0f, texel11[3] / 255.0f);

    // Interpolate in the x direction (horizontal)
    Color c0 = c00 * wx0 + c01 * wx1;  // Left side interpolation
    Color c1 = c10 * wx0 + c11 * wx1;  // Right side interpolation

    // Interpolate in the y direction (vertical)
    Color c = c0 * wy0 + c1 * wy1;

    return c;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Advanced Task
  // Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248
