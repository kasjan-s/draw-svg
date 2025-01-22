#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvans coordinates.

  this->x = x;
  this->y = y;
  this->span = span; 

  // Compute the transformation matrix
  float scale = 1.0f / (2.0f * span); // Scaling factor
  float translate_x = -(x - span) * scale; // Translation in x
  float translate_y = -(y - span) * scale; // Translation in y

  double data[9] = {
        scale, 0,     translate_x, // Row 1
        0,     scale, translate_y, // Row 2
        0,     0,     1            // Row 3
  };

    // Initialize canvas_to_norm using the constructor
  svg_2_norm = Matrix3x3(data);

}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248
