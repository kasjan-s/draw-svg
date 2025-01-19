#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

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

  Color pixel_color;
  float inv255 = 1.0 / 255.0;
  size_t sample_width = width * sample_rate;

	pixel_color.r = sample_buffer[4 * (sx + sy * sample_width)] * inv255;
	pixel_color.g = sample_buffer[4 * (sx + sy * sample_width) + 1] * inv255;
	pixel_color.b = sample_buffer[4 * (sx + sy * sample_width) + 2] * inv255;
	pixel_color.a = sample_buffer[4 * (sx + sy * sample_width) + 3] * inv255;

	pixel_color = ref->alpha_blending_helper(pixel_color, color);

	sample_buffer[4 * (sx + sy * sample_width)] = (uint8_t)(pixel_color.r * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 1] = (uint8_t)(pixel_color.g * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 2] = (uint8_t)(pixel_color.b * 255);
	sample_buffer[4 * (sx + sy * sample_width) + 3] = (uint8_t)(pixel_color.a * 255);
}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= width) return;
	if (y < 0 || y >= height) return;

	Color pixel_color;
	float inv255 = 1.0 / 255.0;
  
	pixel_color.r = pixel_buffer[4 * (x + y * width)] * inv255;
	pixel_color.g = pixel_buffer[4 * (x + y * width) + 1] * inv255;
	pixel_color.b = pixel_buffer[4 * (x + y * width) + 2] * inv255;
	pixel_color.a = pixel_buffer[4 * (x + y * width) + 3] * inv255;

	pixel_color = ref->alpha_blending_helper(pixel_color, color);

  for (int i = 0; i < sample_rate; ++i) {
    for (int j = 0; j < sample_rate; ++j) {
      int sx = sample_rate * x + i;
      int sy = sample_rate * y + j;

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

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  // pixel_buffer[4 * (sx + sy * width)] = (uint8_t)(color.r * 255);
  // pixel_buffer[4 * (sx + sy * width) + 1] = (uint8_t)(color.g * 255);
  // pixel_buffer[4 * (sx + sy * width) + 2] = (uint8_t)(color.b * 255);
  // pixel_buffer[4 * (sx + sy * width) + 3] = (uint8_t)(color.a * 255);

  fill_pixel(sx, sy, color);
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color,
                                          int width) {

  // Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  // ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);
  bresenham(x0, y0, x1, y1, color, width);

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

void SoftwareRendererImp::bresenham(float x0, float y0, float x1, float y1, Color color, int width, bool reflect) {
  // Switch the point order so we iterate from smaller x coordinate.
  if (x1 - x0 < 0) {
    return bresenham(x1, y1, x0, y0, color, width, reflect);
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
      return bresenham(y0, x0, y1, x1, color, width, true);
    }
  } else {
    // 1 point line.
    if (y1 == y0) {
      return rasterize_point(x0, y0, color);
    }

    // x0 == x1 is equivalent to slope m = infinity. If we reflect, we'll work equivalent line with m' = 0.
    return bresenham(y0, x0, y1, x1, color, width, true);
  }

  int sx0 = get_sample_coordinate(x0, sample_rate);
  int sx1 = get_sample_coordinate(x1, sample_rate);
  int sy0 = get_sample_coordinate(y0, sample_rate);
  int sy1 = get_sample_coordinate(y1, sample_rate);

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
          float sample_x = sample_fraction * i + x;
          float sample_y = sample_fraction * j + y;

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

}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
  
	float inv255 = 1.0 / 255.0;
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
  return pixel_color;
}


} // namespace CS248
