#pragma once

#define PI 3.1415926536
#define TwoPI 6.2831853072

double   _r_input, _r_b, _r_s;
double   _r_input_squared, _r_b_squared, _r_s_squared;
double   _lx, _ly;
double   _sx, _sy;
double   _xo, _yo;
size_t   _nx, _ny, _n;
double*  x;
double*  y;
size_t*  tmp_active_i;
size_t*  tmp_active_j;
size_t   num_points;
size_t   num_active;
size_t*  cell_first_point;
size_t*  cell_second_point;
size_t*  active_i;
size_t*  active_j;
int      num_refinement_levels;
size_t** node_nodes;
size_t** _quad_nodes;
size_t** _delaunay;
size_t*  cell_point;
size_t** cell_point2;
size_t   _s1;
size_t   num_expected;
size_t   _num_fixed_quads;
size_t** _ghost_point;  // For a main points, it stores the ghost points
                        // associated with this main point For a ghost point, it
                        // stores the main point of that ghost point
size_t* _disk_color;  // 1D array tells the color of the disk (0 for blue, 1 for
                      // red)
double* _angles;      // for sorting
double  _min_angle, _max_angle, _min_edge_len, _max_edge_len, _min_ar, _max_ar;
size_t _s_d;  // number of points before inserting points at the circumcenter of
              // triangles with three vertices of the same color
size_t* _active;  // instead of deleting points, we mark it as active or
                  // inactive. (0 for inactive, 1 for active)
double   scale;
double   _alpha;
size_t** _quads;  // store the quad with its 4 vertices
size_t*  _list_of_removed_disks;
bool*    _circumcenter;
size_t   _iter_num, _store_point;
size_t   min_edge_ip1, min_edge_ip2;