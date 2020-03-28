#pragma once
using namespace std;

#include <fstream>

inline void PlotQuadMesh(std::string file_name = "Quad Mesh.ps")
{
    fstream quad_mesh(file_name, ios::out);
    quad_mesh.precision(30);
    fstream CAD("quad.txt", ios::out);

    quad_mesh << "%!PS-Adobe-3.0" << endl;
    quad_mesh << "60 60 scale     % one unit = one Centimeter" << endl;

    double scale_x, scale_y, scale_mesh, shift_x, shift_y;
    scale_x = 11.0 / (_nx * _sx);

    scale_mesh = scale_x;
    shift_x = 0;
    shift_y = 0.05 * (28.0 - (_ny * _sy) * scale_mesh) + 1.5;

    scale = scale_mesh;
    quad_mesh << shift_x << " " << shift_y << " translate" << endl;

    quad_mesh << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    quad_mesh << "{newpath" << endl;
    quad_mesh << " moveto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " closepath" << endl;
    quad_mesh << " 0.03 setlinewidth" << endl;
    quad_mesh << " stroke" << endl;
    quad_mesh << "} def" << endl;

    quad_mesh << "/quad      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    quad_mesh << "{newpath" << endl;
    quad_mesh << " moveto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " lineto" << endl;
    quad_mesh << " closepath" << endl;
    quad_mesh << " 0.001 setlinewidth" << endl;
    quad_mesh << " stroke" << endl;
    quad_mesh << "} def" << endl;


    quad_mesh << "/Times-Roman findfont" << endl;
    quad_mesh << "0.1 scalefont" << endl;
    quad_mesh << "setfont" << endl;

    quad_mesh << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    quad_mesh << "{newpath" << endl;
    quad_mesh << "moveto" << endl;
    quad_mesh << "lineto" << endl;
    quad_mesh << "lineto" << endl;
    quad_mesh << "lineto" << endl;
    quad_mesh << "closepath" << endl;
    quad_mesh << "gsave" << endl;
    quad_mesh << "1.0 setgray fill" << endl;
    quad_mesh << " grestore" << endl;
    quad_mesh << "} def" << endl;

    quad_mesh << "/yellow_disk      % stack: x y r" << endl;
    quad_mesh << "{newpath" << endl;
    quad_mesh << " 0 360 arc closepath" << endl;
    quad_mesh << "  1.0 1.0 0.0 setrgbcolor" << endl;
    quad_mesh << " fill" << endl;
    quad_mesh << " 0.01 setlinewidth" << endl;
    quad_mesh << " stroke" << endl;
    quad_mesh << "} def" << endl;

    CAD << "l	" << 0 << "," << 0 << "," << 0 << "		" << 1 << "," << 0
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 1 << "," << 0 << "," << 0 << "		" << 1 << "," << 1
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 1 << "," << 1 << "," << 0 << "		" << 0 << "," << 1
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 0 << "," << 1 << "," << 0 << "		" << 0 << "," << 0
        << "," << 0 << endl
        << endl;

    /*
    quad_mesh << (0-_xo)* scale_mesh << "  " << (0-_yo)*scale_mesh << "  ";
    quad_mesh << (0-_xo)*  scale_mesh << "  " << (1-_yo) * scale_mesh << "  ";
    quad_mesh << (1-_xo)* scale_mesh << "  " << (1-_yo)* scale_mesh << "  ";
    quad_mesh << (1-_xo)* scale_mesh << "  " << (0-_yo)*scale_mesh << "  ";
    quad_mesh << "quad_bold"      << endl;
    */

    double xc, yc, xg, yg;


    for (size_t g = 0; g < num_points; g++) {
        if (_active[g] == 0) {
            continue;
        }

        xc = x[g];
        yc = y[g];

        /*quad_mesh << "/Times-Roman findfont" << endl;
        quad_mesh << "0.15 scalefont" << endl;
        quad_mesh << "0 0 0 setrgbcolor" << endl;
        quad_mesh << "setfont" << endl;
        quad_mesh << (xc-_xo) * scale_mesh <<" "<<(yc-_yo) * scale_mesh<<"
        moveto" << endl; quad_mesh << "("<<g; quad_mesh << ") " << "show" <<
        endl;*/

        quad_mesh << " 0 "
                  << " 0.0 "
                  << "0 setrgbcolor" << endl;
        quad_mesh << "0.006 setlinewidth" << endl;

        for (size_t k = 1; k <= _quad_nodes[g][0]; k++) {

            xg = x[_quad_nodes[g][k]];
            yg = y[_quad_nodes[g][k]];

            quad_mesh << (xc - _xo) * scale_mesh << " "
                      << (yc - _yo) * scale_mesh << " "
                      << "moveto" << endl;
            quad_mesh << (xg - _xo) * scale_mesh << " "
                      << (yg - _yo) * scale_mesh << " "
                      << "lineto" << endl;

            CAD << "l	" << xc << "," << yc << "," << 0 << "		" << xg
                << "," << yg << "," << 0 << endl
                << endl;
        }

        quad_mesh << "stroke" << endl;
    }

    for (size_t i = 0; i < num_points; i++) {
        if (_active[i] == 0) {
            continue;
        }
        if (_disk_color[i] == 1) {
            quad_mesh << "1 0 0 setrgbcolor" << endl;
            quad_mesh << (x[i] - _xo) * scale_mesh << " "
                      << (y[i] - _yo) * scale_mesh << " "
                      << 0.1 * _r_input * scale << " "
                      << "0 360 arc" << endl;  // point plotting
            // quad_mesh << ((x[i]-_xo)+(0.02*_r_input)) * scale_mesh << " " <<
            // ((y[i]-_yo)+(0.2*_r_input)) * scale_mesh << " moveto
            // ("<<i<<")show"<<endl;
            quad_mesh << "fill" << endl;
        } else {
            quad_mesh << "0 0 1 setrgbcolor" << endl;
            quad_mesh << (x[i] - _xo) * scale_mesh << " "
                      << (y[i] - _yo) * scale_mesh << " "
                      << 0.1 * _r_input * scale << " "
                      << "0 360 arc" << endl;  // point plotting
            // quad_mesh << ((x[i]-_xo)+(0.02*_r_input)) * scale_mesh << " " <<
            // ((y[i]-_yo)+(0.2*_r_input)) * scale_mesh << " moveto
            // ("<<i<<")show"<<endl;
            quad_mesh << "fill" << endl;
        }
    }


    if (false) {
        for (size_t V = 0; V < num_points; V++) {
            if (_active[V] == 0) {
                continue;
            }
            if (V >= _s_d) {
                quad_mesh << "1.0 0.5 0 setrgbcolor" << endl;
                quad_mesh << (x[V] - _xo) * scale_mesh << " "
                          << (y[V] - _yo) * scale_mesh << " " << 0.8 * _r_input
                          << " "
                          << "0 360 arc" << endl;  // point plotting
                // quad_mesh << ((x[i]-_xo)+(0.2*_r_input)) * scale_mesh << " "
                // << ((y[i]-_yo)+(0.2*_r_input)) * scale_mesh << " moveto
                // ("<<i<<")show"<<endl;
                quad_mesh << "fill" << endl;
            }
        }
    }


    quad_mesh << "0 0 0 setrgbcolor" << endl;
    quad_mesh << (0 - _xo) * scale_mesh << "  " << (0 - _yo) * scale_mesh
              << "  ";
    quad_mesh << (0 - _xo) * scale_mesh << "  " << (1 - _yo) * scale_mesh
              << "  ";
    quad_mesh << (1 - _xo) * scale_mesh << "  " << (1 - _yo) * scale_mesh
              << "  ";
    quad_mesh << (1 - _xo) * scale_mesh << "  " << (0 - _yo) * scale_mesh
              << "  ";
    quad_mesh << "quad_bold" << endl;

    if (false) {
        quad_mesh << (-100 - _xo) * scale << " " << (0 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (0 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (-100 - _yo) * scale << " "
                  << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
                  << " quad_white" << endl;
        quad_mesh << (0 - _xo) * scale << " " << (-100 - _yo) * scale << " "
                  << (0 - _xo) * scale << " " << (100 - _yo) * scale << " "
                  << (-100 - _xo) * scale << " " << (100 - _yo) * scale << " "
                  << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
                  << " quad_white" << endl;
        quad_mesh << (-100 - _xo) * scale << " " << (1 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (1 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
                  << (-100 - _xo) * scale << " " << (100 - _yo) * scale
                  << " quad_white" << endl;
        quad_mesh << (1 - _xo) * scale << " " << (-100 - _yo) * scale << " "
                  << (1 - _xo) * scale << " " << (100 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
                  << (100 - _xo) * scale << " " << (-100 - _yo) * scale
                  << " quad_white" << endl;
    }

    quad_mesh.close();
}
inline void PlotMesh()
{

    fstream mesh("mesh.ps", ios::out);
    mesh.precision(30);

    mesh << "%!PS-Adobe-3.0" << endl;
    mesh << "60 60 scale     % one unit = one Centimeter" << endl;

    fstream CAD("delaunay.txt", ios::out);

    double scale_x, scale_y, scale_mesh, shift_x, shift_y;
    scale_x = 11.0 / (_nx * _sx);

    scale_mesh = scale_x;
    shift_x = 0;

    shift_y = 0.05 * (28.0 - (_ny * _sy) * scale_mesh) + 1.5;
    scale = scale_mesh;
    mesh << shift_x << " " << shift_y << " translate" << endl;

    mesh << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    mesh << "{newpath" << endl;
    mesh << " moveto" << endl;
    mesh << " lineto" << endl;
    mesh << " lineto" << endl;
    mesh << " lineto" << endl;
    mesh << " closepath" << endl;
    mesh << " 0.03 setlinewidth" << endl;
    mesh << " stroke" << endl;
    mesh << "} def" << endl;

    mesh << "/quad      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    mesh << "{newpath" << endl;
    mesh << " moveto" << endl;
    mesh << " lineto" << endl;
    mesh << " lineto" << endl;
    mesh << " lineto" << endl;
    mesh << " closepath" << endl;
    mesh << " 0.001 setlinewidth" << endl;
    mesh << " stroke" << endl;
    mesh << "} def" << endl;

    mesh << "/Times-Roman findfont" << endl;
    mesh << "0.1 scalefont" << endl;
    mesh << "setfont" << endl;

    mesh << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    mesh << "{newpath" << endl;
    mesh << "moveto" << endl;
    mesh << "lineto" << endl;
    mesh << "lineto" << endl;
    mesh << "lineto" << endl;
    mesh << "closepath" << endl;
    mesh << "gsave" << endl;
    mesh << "1.0 setgray fill" << endl;
    mesh << " grestore" << endl;
    mesh << "} def" << endl;

    mesh << "/yellow_disk      % stack: x y r" << endl;
    mesh << "{newpath" << endl;
    mesh << " 0 360 arc closepath" << endl;
    mesh << "  1.0 1.0 0.0 setrgbcolor" << endl;
    mesh << " fill" << endl;
    mesh << " 0.01 setlinewidth" << endl;
    mesh << " stroke" << endl;
    mesh << "} def" << endl;

    double xc, yc, xg, yg;
    CAD << "l	" << 0 << "," << 0 << "," << 0 << "		" << 1 << "," << 0
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 1 << "," << 0 << "," << 0 << "		" << 1 << "," << 1
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 1 << "," << 1 << "," << 0 << "		" << 0 << "," << 1
        << "," << 0 << endl
        << endl;
    CAD << "l	" << 0 << "," << 1 << "," << 0 << "		" << 0 << "," << 0
        << "," << 0 << endl
        << endl;

    for (size_t g = 0; g < num_points; g++) {

        xc = x[g];
        yc = y[g];

        mesh << " 0 "
             << " 0 "
             << "0 setrgbcolor" << endl;
        // mesh<< "0.006 setlinewidth"<<endl;

        if (false) {

            mesh << " 0 "
                 << " 0 "
                 << "0 setrgbcolor" << endl;
            mesh << " 0 "
                 << " 0 "
                 << "0 setrgbcolor" << endl;
            mesh << (xc - _xo) * scale_mesh << " " << (yc - _yo) * scale_mesh
                 << " moveto (" << g << ")show" << endl;
        }


        mesh << "0.01 setlinewidth" << endl;

        for (size_t k = 1; k <= node_nodes[g][0]; k++) {

            xg = x[node_nodes[g][k]];
            yg = y[node_nodes[g][k]];

            mesh << (xc - _xo) * scale_mesh << " " << (yc - _yo) * scale_mesh
                 << " "
                 << "moveto" << endl;
            mesh << (xg - _xo) * scale_mesh << " " << (yg - _yo) * scale_mesh
                 << " "
                 << "lineto" << endl;

            CAD << "l	" << xc << "," << yc << "," << 0 << "		" << xg
                << "," << yg << "," << 0 << endl
                << endl;
        }

        mesh << "stroke" << endl;
    }

    for (size_t i = 0; i < num_points; i++) {


        if (_active[i] == 0) {
            continue;
        }
        if (_disk_color[i] == 1) {
            mesh << "1 0 0 setrgbcolor" << endl;
            mesh << (x[i] - _xo) * scale_mesh << " "
                 << (y[i] - _yo) * scale_mesh << " "
                 << 0.08 * _r_input * scale_mesh << " "
                 << "0 360 arc" << endl;  // point plotting
            mesh << "fill" << endl;
            mesh << "stroke" << endl;
            if (false) {
                mesh << "0.1 setlinewidth" << endl;
                mesh << ((x[i] - _xo)) * scale_mesh << " "
                     << ((y[i] - _yo)) * scale_mesh << " moveto (" << i
                     << ")show" << endl;
                mesh << "fill" << endl;
            }


        } else {
            mesh << "0 0 1 setrgbcolor" << endl;
            mesh << (x[i] - _xo) * scale_mesh << " "
                 << (y[i] - _yo) * scale_mesh << " "
                 << 0.08 * _r_input * scale_mesh << " "
                 << "0 360 arc" << endl;  // point plotting
            mesh << "fill" << endl;
            mesh << "stroke" << endl;
            if (false) {
                mesh << "0.1 setlinewidth" << endl;
                mesh << ((x[i] - _xo)) * scale_mesh << " "
                     << ((y[i] - _yo)) * scale_mesh << " moveto (" << i
                     << ")show" << endl;
                mesh << "fill" << endl;
            }
        }
    }

    mesh << "0 0 0 setrgbcolor" << endl;
    mesh << (0 - _xo) * scale_mesh << "  " << (0 - _yo) * scale_mesh << "  ";
    mesh << (0 - _xo) * scale_mesh << "  " << (1 - _yo) * scale_mesh << "  ";
    mesh << (1 - _xo) * scale_mesh << "  " << (1 - _yo) * scale_mesh << "  ";
    mesh << (1 - _xo) * scale_mesh << "  " << (0 - _yo) * scale_mesh << "  ";
    mesh << "quad_bold" << endl;

    if (false) {
        mesh << (-100 - _xo) * scale << " " << (0 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (0 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
        mesh << (0 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (0 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
        mesh << (-100 - _xo) * scale << " " << (1 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (1 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (100 - _yo) * scale
             << " quad_white" << endl;
        mesh << (1 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (1 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
    }

    mesh.close();
}
inline void PLOT(double r_in, double r_out, size_t iref, double ssx, double ssy)
{
    fstream file("MPS1.ps", ios::out);
    file.precision(30);
    file << "%!PS-Adobe-3.0" << endl;
    file << "60 60 scale     % one unit = one Centimeter" << endl;

    double xmin(0.0), ymin(0.0);
    double scale_x, scale_y;
    double shift_x, shift_y;

    scale_x = 11.0 / (_nx * _sx);
    scale_y = 10.0 / (_ny * _sy);

    scale = scale_x;
    shift_x = 0;
    shift_y = 0.05 * (28.0 - (_ny * _sy) * scale) + 1.5;


    file << shift_x << " " << shift_y << " translate" << std::endl;

    file << "/Times-Roman findfont" << endl;
    file << "0.2 scalefont" << endl;
    file << "setfont" << endl;


#pragma region Definitions of Shapes:
    file << "/seg      % stack: x1 y1 x2 y2" << std::endl;
    file << "{newpath" << std::endl;
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.05 setlinewidth" << std::endl;
    file << " 1 0 0 setrgbcolor" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;


    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl;
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.0 0 0.0 setrgbcolor" << endl;
    file << " 0.009 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;


    file << "/cell  % stack: x1 y1 x2 y2 x3 y3 x4 y4 " << endl;
    file << "{" << endl;
    file << " newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << "closepath" << endl;
    file << " 0.5 0.5 0.5 setrgbcolor" << endl;
    file << "fill" << endl;
    file << "} def" << endl;


    file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.03 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;


    file << "/pink_disk      % stack: x y r" << endl;
    file << "{newpath" << endl;
    file << " 0 360 arc closepath" << endl;
    file << "  0 1.0 0.17 setrgbcolor" << endl;
    file << " fill" << endl;
    file << " 0.008 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/yellow_disk      % stack: x y r" << endl;
    file << "{newpath" << endl;
    file << " 0 360 arc closepath" << endl;
    file << "  1.0 1.0 0.0 setrgbcolor" << endl;
    file << " fill" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;


    file << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << "moveto" << endl;
    file << "lineto" << endl;
    file << "lineto" << endl;
    file << "lineto" << endl;
    file << "closepath" << endl;
    file << "gsave" << endl;
    file << "1.0 setgray fill" << endl;
    file << " grestore" << endl;
    file << "} def" << endl;


#pragma endregion

    // plot domain boundaries in bold


    /*file <<(0)*scale<< "  " << (0)*scale << "  ";
    file <<(0)*scale<< "  " << (1) * scale << "  ";
    file <<(1)*scale<< "  " << (1)* scale << "  ";
    file <<(1)*scale<< "  " << (0)*scale << "  ";
    file << "quad_bold"      << endl;*/

    /*
   for(size_t i = 0 ;i < num_points ; i++){

       file << "/Times-Roman findfont" << endl;
       file << "0.15 scalefont" << endl;
       file << "1.0 0 0 setrgbcolor" << endl;
       file << "setfont" << endl;
       file << (x[i]-_xo) * scale <<" "<<(y[i]-_yo) * scale<<" moveto" << endl;
       file << "("<<i;
       file << ") " << "show" << endl;


       file << "0 0 0 setrgbcolor" << endl;
       file << "0.05 setlinewidth"<< endl;
       file << (x[i]-_xo) * scale<<" "<< (y[i]-_yo) * scale<<" "<< "0.01 0 360
   arc"<< endl; //point plotting
       //file << X[i] * scale << " " << Y[i] * scale << " moveto
   ("<<i<<")show"<<endl; file << "fill" << endl; file << "0.01
   setlinewidth"<<endl; file << (x[i]-_xo) * scale <<" "<< (y[i]-_yo) * scale
   <<" "<< _r_input * scale<< " 0 360 arc closepath"<< endl; // disk plotting
       file << "stroke"<<endl;

   }
   */


    // plot the grid


    if (iref != 0) {
        iref -= 1;
        size_t iactive, ii, jj;
        double x1, y1, x2, y2;

        for (iactive = 0; iactive < num_active; iactive++) {
            ii = active_i[iactive];
            jj = active_j[iactive];
            x1 = _xo + ii * ssx;
            y1 = _yo + jj * ssy;
            x2 = x1 + ssx;
            y2 = y1 + ssy;
            file << (x1 - _xo) * scale << "  " << (y1 - _yo) * scale << "  ";
            file << (x1 - _xo) * scale << "  " << (y2 - _yo) * scale << "  ";
            file << (x2 - _xo) * scale << "  " << (y2 - _yo) * scale << "  ";
            file << (x2 - _xo) * scale << "  " << (y1 - _yo) * scale << "  ";
            file << "cell" << endl;
        }

        for (iactive = 0; iactive < num_active; iactive++) {
            ii = active_i[iactive];
            jj = active_j[iactive];
            x1 = _xo + ii * ssx;
            y1 = _yo + jj * ssy;
            x2 = x1 + ssx;
            y2 = y1 + ssy;
            file << (x1 - _xo) * scale << "  " << (y1 - _yo) * scale << "  ";
            file << (x1 - _xo) * scale << "  " << (y2 - _yo) * scale << "  ";
            file << (x2 - _xo) * scale << "  " << (y2 - _yo) * scale << "  ";
            file << (x2 - _xo) * scale << "  " << (y1 - _yo) * scale << "  ";
            file << "quad" << endl;
        }
    }


    for (size_t i = 0; i < num_points; i++) {
        if (_active[i] == 0) {
            continue;
        }
        // file<< (x[i]-_xo)*scale<<" "<< (y[i]-_yo)*scale<<" "<<r_in*scale<<"
        // pink_disk"<<endl;
    }

    for (size_t i = 0; i < num_points; i++) {
        if (_active[i] == 0) {
            continue;
        }
        // file << (x[i]-_xo) * scale <<" "<< (y[i]-_yo) * scale <<" "<< r *
        // scale<< " 0 360 arc"<< endl; // disk plotting file << "fill" << endl;
        // file << "stroke"<<endl;
        // file << "0 0 0 setrgbcolor" << endl;
        if (_disk_color[i] == 0) {
            file << "0 0 1 setrgbcolor" << endl;
        }  // if red disks
        else {
            file << "1 0 0 setrgbcolor" << endl;
        }  // if blue disk
        file << "0.009 setlinewidth" << endl;
        file << (x[i] - _xo) * scale << " " << (y[i] - _yo) * scale << " "
             << r_in * scale << " 0 360 arc closepath"
             << endl;  // disk plotting
        file << "stroke" << endl;

        if (r_in != r_out) {
            file << (x[i] - _xo) * scale << " " << (y[i] - _yo) * scale << " "
                 << r_out * scale << " 0 360 arc closepath"
                 << endl;  // disk plotting
            file << "stroke" << endl;
        }
        // file << "0 0 0 setrgbcolor" << endl;
        // file << "0.1 setlinewidth"<< endl;
        file << (x[i] - _xo) * scale << " " << (y[i] - _yo) * scale << " "
             << r_in * 0.1 * scale << " "
             << "0 360 arc" << endl;  // point plotting
        file << "fill" << endl;
        file << "stroke" << endl;

        if (false) {
            file << (x[i] - _xo) * scale << " " << (y[i] - _yo) * scale
                 << " moveto (" << i << ")show" << endl;
            file << "fill" << endl;
            file << "stroke" << endl;
        }
    }

    file << "0 0 0 setrgbcolor" << endl;
    file << (0 - _xo) * scale << "  " << (0 - _yo) * scale << "  ";
    file << (0 - _xo) * scale << "  " << (1.0 - _yo) * scale << "  ";
    file << (1.0 - _xo) * scale << "  " << (1.0 - _yo) * scale << "  ";
    file << (1.0 - _xo) * scale << "  " << (0.0 - _yo) * scale << "  ";
    file << "quad_bold" << endl;


    if (false) {
        for (size_t i = 0; i < _nx; i++) {
            double x = _xo + i * _sx;
            for (size_t j = 0; j < _ny; j++) {
                double y = _yo + j * _sy;
                file << (x - _xo) * scale << "  " << (y - _yo) * scale << "  ";
                file << (x + _sx - _xo) * scale << "  " << (y - _yo) * scale
                     << "  ";
                file << (x + _sx - _xo) * scale << "  "
                     << (y + _sy - _yo) * scale << "  ";
                file << (x - _xo) * scale << "  " << (y + _sy - _yo) * scale
                     << "  ";
                file << "quad" << endl;
            }
        }
    }

    if (false) {
        file << (-100 - _xo) * scale << " " << (0 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (0 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
        file << (0 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (0 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
        file << (-100 - _xo) * scale << " " << (1 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (1 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (-100 - _xo) * scale << " " << (100 - _yo) * scale
             << " quad_white" << endl;
        file << (1 - _xo) * scale << " " << (-100 - _yo) * scale << " "
             << (1 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (100 - _yo) * scale << " "
             << (100 - _xo) * scale << " " << (-100 - _yo) * scale
             << " quad_white" << endl;
    }


    file.close();
}

inline void PlotSinglePoint(double xx, double yy)
{
    fstream file("Quad Mesh.ps", ios::app);
    file << "1.0 1.0 0.0 setrgbcolor" << endl;
    file << "0.01 setlinewidth" << endl;

    fstream file1("mesh.ps", ios::app);
    file1 << "1.0 1.0 0.0 setrgbcolor" << endl;
    file1 << "0.01 setlinewidth" << endl;

    fstream file2("MPS1.ps", ios::app);
    file2 << "1.0 1.0 0.0 setrgbcolor" << endl;
    file2 << "0.01 setlinewidth" << endl;

    file << (xx - _xo) * scale << " " << (yy - _yo) * scale << " "
         << 0.05 * _r_input * scale << " yellow_disk" << endl;
    file1 << (xx - _xo) * scale << " " << (yy - _yo) * scale << " "
          << 0.05 * _r_input * scale << " yellow_disk" << endl;
    file2 << (xx - _xo) * scale << " " << (yy - _yo) * scale << " "
          << 0.05 * _r_input * scale << " yellow_disk" << endl;
}