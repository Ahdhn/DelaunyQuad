#include <iostream>
#include <string.h>
#include <time.h>

using namespace std;

// clang-format off
#include "Variables.h"
#include "Plotters.h"
#include "DartThrowing.h"
#include "DelaunayGen.h"
#include "ProcessDelaunay.h"
#include "QuadGen.h"
#include "Rand.h"
#include "PostProcess.h"
// clang-format on


int main()
{
    size_t i, j, d;
    size_t V, in_domain(0);

    _r_input = 0.02;
    _alpha = 2.5;

    _r_input_squared = _r_input * _r_input;

    _lx = 1.0;
    _ly = 1.0;

    _sx = _r_input / sqrt(2.0);
    _sy = _sx;

    _nx = size_t(ceil(_lx / _sx) + 44.0);
    _ny = _nx;
    _n = _nx * _ny;

    _xo = -22.0 * _sx;
    _yo = -22.0 * _sy;


    num_refinement_levels = 30;
    num_active = 0;
    num_points = 0;

    // srand ( time(NULL) ); // activate for different experiments
    InitiateRandNumGenerator(rand());

    num_expected = size_t(10E6);

    x = new double[num_expected];
    y = new double[num_expected];

    _angles = new double[20];

    _ghost_point = new size_t*[num_expected];
    _disk_color = new size_t[num_expected];
    _active = new size_t[num_expected];

    for (i = 0; i < num_expected; i++) {
        _ghost_point[i] = new size_t[5];
        _ghost_point[i][0] = 0;
        _active[i] = 1;
    }

    num_expected = _n;

    active_i = new size_t[num_expected];
    active_j = new size_t[num_expected];
    tmp_active_i = new size_t[num_expected];
    tmp_active_j = new size_t[num_expected];

    cell_point2 = new size_t*[_n];

    for (i = 0; i < _n; i++) {
        cell_point2[i] = new size_t[9];
        cell_point2[i][0] = 0;
    }

    for (size_t icell = 0; icell < _n; icell++) {
        i = size_t(icell / _ny);
        j = size_t(icell - i * _ny);

        if (i >= 10 && i < _ny - 10 && j >= 10 && j < _nx - 10) {
            active_i[num_active] = i;
            active_j[num_active] = j;
            num_active++;
            if (num_active > num_expected - 5) {
                cout << "Error Do something (2).. We are losing the baby "
                     << endl;
                system("pause");
            }
        }
    }

    _r_s = _r_input;  //=r_in which define the min seperation bewteen disks of
                      // different colors and the coverage
    _r_s_squared = _r_s * _r_s;

    _r_b = _r_s * _alpha;  //=r_out which define the min seperation between
                           // disks of the same color
    _r_b_squared = _r_b * _r_b;


    _list_of_removed_disks = new size_t[5000];
    _list_of_removed_disks[0] = 0;


    _circumcenter = new bool[num_expected];
    for (V = 0; V < num_expected; V++) {
        _circumcenter[0] = false;
    }

    DartThrowingWithColors(_r_s, _r_s_squared, _r_b_squared);

    if (true) {
        PLOT(_r_s, _r_b, 0, _sx, _sy);
    }

    delete[] active_i;
    delete[] active_j;
    delete[] tmp_active_i;
    delete[] tmp_active_j;

    _s_d = num_points;


    size_t mn = 3 * num_points;
    node_nodes = new size_t*[mn];  // initiation

    for (i = 0; i < mn; i++) {
        node_nodes[i] = new size_t[15 + 1];
        node_nodes[i][0] = 0;
    }

    DelaunayMeshGenerator();
    for (size_t AM = 0; AM < num_points; AM++) {
        SortDelaunay2(AM);
    }

    if (true) {
        PlotMesh();
    }


    CheckDelaunay2();


    if (true) {
        PlotMesh();
    }


    mn = 4 * num_points;
    _quads = new size_t*[mn];
    for (i = 0; i < mn; i++) {
        _quads[i] = new size_t[4];
    }

    _quad_nodes = new size_t*[mn];  // initiation
    for (i = 0; i < mn; i++) {
        _quad_nodes[i] = new size_t[10 + 1];
        _quad_nodes[i][0] = 0;
    }


    QuadMesh_Constructor2();
    for (V = 0; V < num_points; V++) {
        SortQuadMesh2(V);
    }
    if (true) {
        PlotQuadMesh();
    }

    PostProcess();

    if (true) {
        PlotQuadMesh("Quad Mesh postprocess.ps");
    }

    return 0;
}