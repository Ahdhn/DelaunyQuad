#include <string.h>
#include <time.h>

#include <iostream>

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
    init();

    // srand ( time(NULL) ); // activate for different experiments
    InitiateRandNumGenerator(rand());

    // 1)  Dart throwing
    DartThrowingWithColors(_r_s, _r_s_squared, _r_b_squared);
    // plot the samples
    PLOT(_r_s, _r_b, 0, _sx, _sy);


    _s_d = _num_points;
    const size_t nn_num = 3 * _num_points;
    _node_nodes = new size_t*[nn_num];
    for (size_t i = 0; i < nn_num; i++) {
        _node_nodes[i] = new size_t[15 + 1];
        _node_nodes[i][0] = 0;
    }

    // 2) Generate Delaunay mesh
    DelaunayMeshGenerator();
    // plot delaunay mesh
    PlotMesh();


    // 3) Process Delaunay mesh
    CheckDelaunay2();
    PlotMesh();


    const size_t num_quads = 4 * _num_points;
    _quads = new size_t*[num_quads];
    for (size_t i = 0; i < num_quads; i++) {
        _quads[i] = new size_t[4];
    }

    _quad_nodes = new size_t*[num_quads];
    for (size_t i = 0; i < num_quads; i++) {
        _quad_nodes[i] = new size_t[10 + 1];
        _quad_nodes[i][0] = 0;
    }

    // 4) Generate Quad mesh
    QuadMesh_Constructor2();
    // plot quad mesh
    PlotQuadMesh();


    // 5) Post process quad mesh
    PostProcess();
    PlotQuadMesh("Quad Mesh postprocess.ps");


    cleanup();
    return 0;
}