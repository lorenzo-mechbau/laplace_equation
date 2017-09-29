gfx read node laplace_equation.part0.exnode;
gfx read elem laplace_equation.part0.exelem;

gfx def faces egroup laplace_equation;

gfx modify g_element laplace_equation general clear circle_discretization 6 default_coordinate Coordinate element_discretization "4*4*4" native_discretization none;
gfx modify g_element laplace_equation lines select_on invisible material default selected_material default_selected;
gfx modify g_element laplace_equation surfaces select_on material default data Phi spectrum default selected_material default_selected render_shaded;
gfx modify g_element laplace_equation node_points glyph sphere general size "0.025*0.025*0.025" centre 0,0,0 font default label cmiss_number select_on material blue selected_material default_selected;
gfx modify g_element laplace_equation cylinders constant_radius 0.005 select_on material gold selected_material default_selected render_shaded;

gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 567 height 653;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point -1.32423 -1.97615 2.23909 interest_point 0.5 0.5 0.5 up_vector 0.305326 0.386341 0.870354 view_angle 37.5028 near_clipping_plane 0.0353321 far_clipping_plane 12.6265 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;
