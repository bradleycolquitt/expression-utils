flip_geometry_y = function(x, max_val) {
  (x - c(0, max_val)) * c(1,-1)
}

flip_geometry_xy = function(x, max_val) {
  (x - c(max_val, max_val)) * c(-1,-1)
}

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
rotate_image = function(x, max_val, alpha) {
  centroid = c(max_val/2, max_val/2)
  ((x-centroid) * rot(alpha)) + centroid
}