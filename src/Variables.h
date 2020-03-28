#pragma once

#define PI 3.1415926536
#define TwoPI 6.2831853072

const double _r_input(0.02);
const double _r_s(_r_input);
const double _alpha(2.5);
const double _r_b(_r_s* _alpha);
const double _r_input_squared(_r_input*_r_input), _r_b_squared(_r_b*_r_b),
    _r_s_squared(_r_s*_r_s);
const double _lx(1.0), _ly(1.0);
const double _sx(_r_input / sqrt(2.0)), _sy(_r_input / sqrt(2.0));
const double _xo(-22.0 * _sx), _yo(-22.0 * _sy);
const size_t _nx(size_t(ceil(_lx / _sx) + 44.0)), _ny(_nx), _n(_nx*_ny);
const int    _num_refinement_levels(30);
size_t       _num_expected(10E6);

double* x;
double* y;
size_t* _tmp_active_i;
size_t* _tmp_active_j;
size_t  _num_points;
size_t  num_active;

size_t* _active_i;
size_t* _active_j;

size_t** _node_nodes;
size_t** _quad_nodes;
size_t** _delaunay;
size_t** _cell_point2;

size_t _num_fixed_quads;

size_t** _ghost_point;  // For a main points, it stores the ghost points
                        // associated with this main point For a ghost point, it
                        // stores the main point of that ghost point
size_t* _disk_color;  // 1D array tells the color of the disk (0 for blue, 1 for
                      // red)

double* _angles;  // for sorting

size_t _s_d;  // number of points before inserting points at the circumcenter of
              // triangles with three vertices of the same color

size_t* _active;  // instead of deleting points, we mark it as active or
                  // inactive. (0 for inactive, 1 for active)
double scale;

size_t** _quads;  // store the quad with its 4 vertices
size_t*  _list_of_removed_disks;
bool*    _circumcenter;


inline void init()
{
    num_active = 0;
    _num_points = 0;

    x = new double[_num_expected];
    y = new double[_num_expected];

    _angles = new double[20];

    _ghost_point = new size_t*[_num_expected];
    _disk_color = new size_t[_num_expected];
    _active = new size_t[_num_expected];

    for (size_t i = 0; i < _num_expected; i++) {
        _ghost_point[i] = new size_t[5];
        _ghost_point[i][0] = 0;
        _active[i] = 1;
    }

    _num_expected = _n;

    _active_i = new size_t[_num_expected];
    _active_j = new size_t[_num_expected];
    _tmp_active_i = new size_t[_num_expected];
    _tmp_active_j = new size_t[_num_expected];

    _cell_point2 = new size_t*[_n];

    for (size_t i = 0; i < _n; i++) {
        _cell_point2[i] = new size_t[9];
        _cell_point2[i][0] = 0;
    }

    for (size_t icell = 0; icell < _n; icell++) {
        size_t i = size_t(icell / _ny);
        size_t j = size_t(icell - i * _ny);

        if (i >= 10 && i < _ny - 10 && j >= 10 && j < _nx - 10) {
            _active_i[num_active] = i;
            _active_j[num_active] = j;
            num_active++;
            if (num_active > _num_expected - 5) {
                cout << "Error Do something (2).. We are losing the baby "
                     << endl;
                system("pause");
            }
        }
    }

    _list_of_removed_disks = new size_t[_num_expected];
    _list_of_removed_disks[0] = 0;

    _circumcenter = new bool[_num_expected];
    for (size_t i = 0; i < _num_expected; i++) {
        _circumcenter[i] = false;
    }
}

inline void cleanup()
{

    delete[] _active_i;
    delete[] _active_j;
    delete[] _tmp_active_i;
    delete[] _tmp_active_j;
    delete[] x;
    delete[] y;

    for (size_t i = 0; i < _num_expected; ++i) {
        delete[] _ghost_point[i];
    }
    delete[] _ghost_point;
    for (size_t i = 0; i < _n; ++i) {
        delete[] _cell_point2[i];
    }

    delete[] _cell_point2;
    delete[] _disk_color;
    delete[] _active;
    delete[] _list_of_removed_disks;
    delete[] _circumcenter;
    
    for (size_t i = 0; i < 4 * _num_points; i++) {
        delete[] _quads[i];
        delete[] _quad_nodes[i];
    }
    delete[] _quads;
    delete[] _quad_nodes;
}