#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "texture.h"

#include "triangulation.h"

using namespace std;

namespace CS248 {

namespace {
  int get_sample_coordinate(float coord, int sample_rate) {
    float remainder = coord - floor(coord);
    int extra = static_cast<int>(remainder * sample_rate);
    return sample_rate * coord + extra;
  }
}



// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function
  // check bounds
  if (sx < 0 || sx >= sample_rate * width) return;
  if (sy < 0 || sy >= sample_rate * height) return;

  Color sample_color;
  float inv255 = 1.0 / 255.0;
  size_t sample_width = width * sample_rate;

	sample_color.r = sample_buffer[4 * (sx + sy * sample_width)] * inv255;
	sample_color.g = sample_buffer[4 * (sx + sy * sample_width) + 1] * inv255;
	sample_color.b = sample_buffer[4 * (sx + sy * sample_width) + 2] * inv255;
	sample_color.a = sample_buffer[4 * (sx + sy * sample_width) + 3] * inv255;

	sample_color = alpha_blending(sample_color, color);

	sample_buffer[4 * (sx + sy * sample_width)] = (uint8_t)(sample_color.r * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 1] = (uint8_t)(sample_color.g * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 2] = (uint8_t)(sample_color.b * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 3] = (uint8_t)(sample_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// // check bounds
	if (x < 0 || x >= width) return;
	if (y < 0 || y >= height) return;

  float inv255 = 1.0 / 255.0;

  for (int i = 0; i < sample_rate; ++i) {
    for (int j = 0; j < sample_rate; ++j) {
      Color pixel_color;

      int sx = sample_rate * x + i;
      int sy = sample_rate * y + j;

      pixel_color.r = sample_buffer[4 * (sx + sy * sample_rate * width)] * inv255;
      pixel_color.g = sample_buffer[4 * (sx + sy * sample_rate * width) + 1] * inv255;
      pixel_color.b = sample_buffer[4 * (sx + sy * sample_rate * width) + 2] * inv255;
      pixel_color.a = sample_buffer[4 * (sx + sy * sample_rate * width) + 3] * inv255;

	    pixel_color = alpha_blending(pixel_color, color);

      sample_buffer[4 * (sx + sy * sample_rate * width)] = (uint8_t)(pixel_color.r * 255);
      sample_buffer[4 * (sx + sy * sample_rate * width) + 1] = (uint8_t)(pixel_color.g * 255);
      sample_buffer[4 * (sx + sy * sample_rate * width) + 2] = (uint8_t)(pixel_color.b * 255);
      sample_buffer[4 * (sx + sy * sample_rate * width) + 3] = (uint8_t)(pixel_color.a * 255);
    }
  }
}

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  update_sample_buffer();

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  update_sample_buffer();

}

void SoftwareRendererImp::set_pixel_buffer( unsigned char* pixel_buffer,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;
  update_sample_buffer();

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack
  Matrix3x3 previous_transformation = transformation;
  transformation = transformation * (element->transform);

	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}

  transformation = previous_transformation;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point(float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = static_cast<int>(round(x - 0.5));
  int sy = static_cast<int>(round(y - 0.5));

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  fill_pixel(sx, sy, color);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color,
                                          int width) {

  // Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  // ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);
  bresenham(x0, y0, x1, y1, color);

  // bresenham_width(x0, y0, x1, y1, color, width);
  
  // xiaolin(x0, y0, x1, y1, color);

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

void SoftwareRendererImp::bresenham(float x0, float y0, float x1, float y1, Color color, bool reflect) {
  // Switch the point order so we iterate from smaller x coordinate.
  if (x1 - x0 < 0) {
    return bresenham(x1, y1, x0, y0, color, reflect);
  }

  int float_sign = 1;
  if (x1 != x0) {
    // Slope value.
    float m = (y1-y0)/(x1-x0);

    if (m < 0) {
      float_sign = -1;
    }

    if (abs(m) > 1) {
      // Bresenham algorithm works for slopes with abs value less or equal to 1. 
      // For higher slopes we are reflecting the line around y=x slope by swapping x0 with y0, and x1 with y1.
      // This allows us to reuse Bresenham's algorithm for quadrants where slope is bigger than 1, as in this new
      // reflected coordinate system new slope value is m' = (x1-x0) / (y1-y0) = 1/m ==> it's between 0 and 1.
      return bresenham(y0, x0, y1, x1, color, true);
    }
  } else {
    // 1 point line.
    if (y1 == y0) {
      return rasterize_point(x0, y0, color);
    }

    // x0 == x1 is equivalent to slope m = infinity. If we reflect, we'll work equivalent line with m' = 0.
    return bresenham(y0, x0, y1, x1, color, true);
  }

  int sx0 = std::round(x0 - 0.5);
  int sx1 = std::round(x1 - 0.5);
  int sy0 = std::round(y0 - 0.5);
  int sy1 = std::round(y1 - 0.5);

  int dx = sx1 - sx0;
  int dy = sy1 - sy0;
  int y = sy0;
  int error = 0;
  for (int x = sx0; x <= sx1; ++x) {
      if (reflect) {
        fill_pixel(y, x, color);
      } else {
        fill_pixel(x, y, color);
    }

    error += float_sign * dy;
    if ((error << 1) >= float_sign * dx) {
      y += float_sign;
      error -= dx;
    }
  }
}

void SoftwareRendererImp::bresenham_width(float x0, float y0, float x1, float y1, Color color, int width, bool reflect) {
  // Switch the point order so we iterate from smaller x coordinate.
  if (x1 - x0 < 0) {
    return bresenham_width(x1, y1, x0, y0, color, width, reflect);
  }

  int float_sign = 1;
  if (x1 != x0) {
    // Slope value.
    float m = (y1-y0)/(x1-x0);

    if (m < 0) {
      float_sign = -1;
    }

    if (abs(m) > 1) {
      // Bresenham algorithm works for slopes with abs value less or equal to 1. 
      // For higher slopes we are reflecting the line around y=x slope by swapping x0 with y0, and x1 with y1.
      // This allows us to reuse Bresenham's algorithm for quadrants where slope is bigger than 1, as in this new
      // reflected coordinate system new slope value is m' = (x1-x0) / (y1-y0) = 1/m ==> it's between 0 and 1.
      return bresenham_width(y0, x0, y1, x1, color, width, true);
    }
  } else {
    // 1 point line.
    if (y1 == y0) {
      return rasterize_point(x0, y0, color);
    }

    // x0 == x1 is equivalent to slope m = infinity. If we reflect, we'll work equivalent line with m' = 0.
    return bresenham_width(y0, x0, y1, x1, color, width, true);
  }

  int sx0 = get_sample_coordinate(x0 - 0.5, sample_rate);
  int sx1 = get_sample_coordinate(x1 - 0.5, sample_rate);
  int sy0 = get_sample_coordinate(y0 - 0.5, sample_rate);
  int sy1 = get_sample_coordinate(y1 - 0.5, sample_rate);

  int dx = sx1 - sx0;
  int dy = sy1 - sy0;
  int y = sy0;
  int error = 0;
  
  for (int x = sx0; x <= sx1; ++x) {
    int sample_width = sample_rate * width;
    int sy = y - sample_width / 2;
    for (int w = 0; w < sample_width; ++w) {
      if (reflect) {
        fill_sample(sy + w, x, color);
      } else {
        fill_sample(x, sy + w, color);
      }
    }

    error += float_sign * dy;
    if ((error << 1) >= float_sign * dx) {
      y += float_sign;
      error -= dx;
    }
  }
}

void SoftwareRendererImp::xiaolin(float x0, float y0, float x1, float y1, Color color) {
  // First let's go from screen coordinates to pixel coordinates.
  // E.g. point (0.5,0.5) would correspond to the center of pixel (0,0).
  x0 -= 0.5;
  y0 -= 0.5;
  x1 -= 0.5;
  y1 -= 0.5;

  bool steep = abs(y1-y0) > abs(x1-x0);
  if (steep) {
    std::swap(x0, y0);
    std::swap(x1, y1);
  }

  if (x0 > x1) {
    std::swap(x0, x1);
    std::swap(y0, y1);
  }

  float slope = 1.0f;
  if (x1 != x0) {
    slope = (y1 - y0) / (x1 - x0);
  }

  Color pixel_color = color;

  // First endpoint
  int x_pixel_0 = static_cast<int>(round(x0));
  int y_pixel_0 = floor(y0 + slope * (x_pixel_0 - x0));
  float y_val_0 = y0 + slope * (x_pixel_0 - x0);
  float x0_gap = floor(x0 + 0.5) - (x0 - 0.5);

  float y_dist = y_val_0 - floor(y_val_0);
  if (steep) {
    fill_pixel(y_pixel_0, x_pixel_0, (1.0 - y_dist) * x0_gap * color);
    fill_pixel(y_pixel_0 + 1, x_pixel_0, y_dist * x0_gap * color);
  } else {
    fill_pixel(x_pixel_0, y_pixel_0, (1.0 - y_dist) * x0_gap * color);
    fill_pixel(x_pixel_0, y_pixel_0 + 1,  y_dist * x0_gap * color);
  }

  // Second endpoint
  int x_pixel_1 = static_cast<int>(round(x1));
  int y_pixel_1 = floor(y1 + slope * (x_pixel_1 - x1));
  float y_val_1 = y1 + slope * (x_pixel_1 - x1);
  float x1_gap = x1 + 0.5 - floor(x1 + 0.5);

  y_dist = y_val_1 - floor(y_val_1);
  if (steep) {
    fill_pixel(y_pixel_1, x_pixel_1, (1.0 - y_dist) * x1_gap * color);
    fill_pixel(y_pixel_1 + 1, x_pixel_1, y_dist * x1_gap * pixel_color);
  } else {
    fill_pixel(x_pixel_1, y_pixel_1, (1.0 - y_dist) * x1_gap * color);
    fill_pixel(x_pixel_1, y_pixel_1 + 1, y_dist * x1_gap * pixel_color);
  }

  float y = y_val_0 + slope;
  for (int x = x_pixel_0 + 1; x < x_pixel_1; ++x) {
    float y_dist = y - floor(y);

    if (steep) {
      fill_pixel(floor(y), x, (1.0 - y_dist) * pixel_color);
      fill_pixel(floor(y) + 1, x, y_dist * pixel_color);
    } else {
      fill_pixel(x, floor(y), (1.0 - y_dist) * pixel_color);
      fill_pixel(x, floor(y) + 1, y_dist * pixel_color);
    }

    y += slope;
  }
}

struct LineEquation {
  float A;
  float B;
  float C;

  explicit LineEquation(float x0, float y0, float x1, float y1) {
    A = y1 - y0;
    B = -(x1 - x0);
    C = y0 * (x1 - x0) - x0 * (y1 - y0);
  }

  float test(float x, float y) const {
    return A * x + B * y + C;
  }
};

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 1: 
  // Implement triangle rasterization
  float bound_x0 = min(x0, min(x1, x2));
  float bound_y0 = min(y0, min(y1, y2));
  float bound_x1 = max(x0, max(x1, x2));
  float bound_y1 = max(y0, max(y1, y2));

  vector<LineEquation> lines = {
    LineEquation(x0, y0, x1, y1),
    LineEquation(x1, y1, x2, y2),
    LineEquation(x2, y2, x0, y0)
  };

  float signed_area = x0*y1 - x1*y0 + x1*y2 - x2*y1 + x2*y0 - x0*y2;
  int orientation = signed_area >= 0 ? 1 : -1;

  for (int x = bound_x0; x <= bound_x1; ++x) {
    for (int y = bound_y0; y <= bound_y1; ++y) {
      for (int i = 0; i < sample_rate; ++i) {
        for (int j = 0; j < sample_rate; ++j) { 
          float sample_fraction = 1.0f / sample_rate;
          
          // +0.5 offset to convert from sample to screen coordinates.
          float sample_x = sample_fraction * i + x + 0.5;
          float sample_y = sample_fraction * j + y + 0.5;

          if (all_of(lines.begin(), lines.end(), [sample_x, sample_y, orientation](const LineEquation& line) {
              return orientation * line.test(sample_x, sample_y) <= 0;
            })) {
            fill_sample(x * sample_rate + i, y * sample_rate + j, color);
          }
        }
      }
    }
  }

  // Advanced Task
  // Implementing Triangle Edge Rules

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization

  // Create an instance of the sampler (adjust as needed to match your setup)
  Sampler2DImp sampler(SampleMethod::NEAREST);
  //Sampler2DImp sampler(SampleMethod::BILINEAR);

  // Image bounds
  int screen_x0 = static_cast<int>(std::floor(x0));
  int screen_y0 = static_cast<int>(std::floor(y0));
  int screen_x1 = static_cast<int>(std::ceil(x1));
  int screen_y1 = static_cast<int>(std::ceil(y1));

  float sample_fraction = 1.0f / sample_rate;

  // Iterate over the pixels in the bounding box
  for (int y = screen_y0; y <= screen_y1; ++y) {
      for (int x = screen_x0; x <= screen_x1; ++x) {
          for (int i = 0; i < sample_rate; ++i) {
            for (int j = 0; j < sample_rate; ++j) {
              float sample_x = sample_fraction * i + x + 0.5;
              float sample_y = sample_fraction * j + y + 0.5;

              // Map screen space (x, y) to normalized texture coordinates (u, v)
              float u = (sample_x - x0) / (x1 - x0);
              float v = (sample_y - y0) / (y1 - y0);

              // Sample the texture using bilinear interpolation via the sampler
              // Color tex_color = sampler.sample_nearest(tex, u, v, 0);
              Color tex_color = sampler.sample_bilinear(tex, u, v, 0);

              // Blend the sampled color into the framebuffer
              fill_sample(x * sample_rate + i, y * sample_rate + j, tex_color);
            }
          }
      }
  }

}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
  
  float sample_fraction = 1.0f / (sample_rate * sample_rate);

  for (int x = 0; x < width; ++x) {
    for (int y = 0; y < height; ++y) {
      Color pixel_color(0, 0, 0, 0);

      for (int i = 0; i < sample_rate; ++i) {
        for (int j = 0; j < sample_rate; ++j) {
          int sx = sample_rate * x + i;
          int sy = sample_rate * y + j;

          pixel_color.r += sample_fraction * sample_buffer[4 * (sx + sy * sample_rate * width)];
          pixel_color.g += sample_fraction * sample_buffer[4 * (sx + sy * sample_rate * width) + 1];
          pixel_color.b += sample_fraction * sample_buffer[4 * (sx + sy * sample_rate * width) + 2];
          pixel_color.a += sample_fraction * sample_buffer[4 * (sx + sy * sample_rate * width) + 3];
        }
      }

      pixel_buffer[4 * (x + y * width)] = pixel_color.r;
      pixel_buffer[4 * (x + y * width) + 1] = pixel_color.g;
      pixel_buffer[4 * (x + y * width) + 2] = pixel_color.b;
      pixel_buffer[4 * (x + y * width) + 3] = pixel_color.a;
    }
  }

  return;
}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  pixel_color.a = 1.0 - (1.0 - color.a) * (1.0 - pixel_color.a);
  pixel_color.r = (1.0 - color.a) * pixel_color.r + color.r;
  pixel_color.g = (1.0 - color.a) * pixel_color.g + color.g;
  pixel_color.b = (1.0 - color.a) * pixel_color.b + color.b;

  return pixel_color;
}


} // namespace CS248
