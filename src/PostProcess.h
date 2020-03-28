#pragma once
#include "QuadGen.h"

inline double GetAngle(double xa,
                       double ya,
                       double xb,
                       double yb,
                       double xc,
                       double yc)
{
    // return angle in degrees between b.a.c

    double dxba, dyba, dxca, dyca, dot, pcross, theta;

    dxba = xb - xa;
    dyba = yb - ya;

    dxca = xc - xa;
    dyca = yc - ya;

    dot = dxba * dxca + dyba * dyca;
    pcross = dxba * dyca - dyba * dxca;
    theta = atan2(pcross, dot);

    if (theta < 0) {
        theta += TwoPI;
    }

    return theta * 180 / PI;

}

inline void Move(double  Ax,
                 double  Ay,
                 double  Bx,
                 double  By,
                 double& Cx,
                 double& Cy)
{
    double lenAB = sqrt((Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By));
    double dist = 0.1 * _r_s;
    Cx = Ax - (Ax - Bx) / lenAB * dist;
    Cy = Ay - (Ay - By) / lenAB * dist;
}

inline void MoveThePoint(size_t ip1,
                         size_t ip,
                         size_t ip2)  // move the point for random mesh
{
    double xn, yn, dx, dy;
    size_t ig;
    for (size_t node = 1; node <= _quad_nodes[ip][0]; node++) {
        if (_quad_nodes[ip][node] == ip1) {
            continue;
        }
        if (_quad_nodes[ip][node] == ip2) {
            continue;
        }

        Move(x[ip], y[ip], x[_quad_nodes[ip][node]], y[_quad_nodes[ip][node]],
             xn, yn);
        dx = xn - x[ip];
        dy = yn - y[ip];

        x[ip] = xn;
        y[ip] = yn;

        if (x[ip] >= 0.0 && x[ip] < 1.0 && y[ip] >= 0.0 && y[ip] < 1.0) {
            if (_ghost_point[ip][0] == 0) {
                continue;
            }  // if no ghost
            for (size_t q = 1; q <= _ghost_point[ip][0]; q++) {
                ig = _ghost_point[ip][q];
                x[ig] = dx + x[ig];
                y[ig] = dy + y[ig];
            }
        }


        // PlotProblem(ip);

        if (GetAngle(x[ip], y[ip], x[ip1], y[ip1], x[ip2], y[ip2]) > 170) {
            MoveThePoint(ip1, ip, ip2);
        }

        break;
    }
}

inline void FixingAngle()
{
    size_t V, ip2, Vnext, i, ip1;

    for (i = 0; i < _num_points; i++) {
        if (!(x[i] >= 0.0 && x[i] < 1.0 && y[i] >= 0.0 && y[i] < 1.0)) {
            continue;
        }

        if (_quad_nodes[i][0] == 2) {
            continue;
        }
        SortQuadMesh2(i);

        for (V = 1; V <= _quad_nodes[i][0]; V++) {
            Vnext = V + 1;
            if (Vnext == _quad_nodes[i][0] + 1) {
                Vnext = 1;
            }

            ip1 = _quad_nodes[i][V];
            ip2 = _quad_nodes[i][Vnext];

            // if angle geater than 170
            if (GetAngle(x[i], y[i], x[ip1], y[ip1], x[ip2], y[ip2]) > 170) {
                if (false) {
                    PlotSinglePoint(x[i], y[i]);
                }

                _num_fixed_quads++;
                MoveThePoint(ip1, i, ip2);  // for trivalent node
                if (false) {
                    PlotQuadMesh();
                }
            }
        }
    }
}

inline void RemoveQuadConnectivity(size_t ip)
{
    // remove (ip) from all points connected to it
    // remove all points connected to (ip)
    size_t V, d, ip1;
    bool   found;
    for (V = 1; V <= _quad_nodes[ip][0]; V++) {
        ip1 = _quad_nodes[ip][V];

        // remove ip from ip1
        found = false;
        for (d = 1; d < _quad_nodes[ip1][0]; d++) {
            if (_quad_nodes[ip1][d] == ip) {
                found = true;
            }  // first find ip
            if (found) {
                _quad_nodes[ip1][d] =
                    _quad_nodes[ip1][d + 1];  // move the point on step to
                                              // preserve the sorted order
            }
        }
        if (!found && _quad_nodes[ip1][d] != ip) {
            cout << "Error at RemoveQuadConnectivity().." << endl;
            system("pause");
        }
        _quad_nodes[ip1][0]--;
    }
    // if(ip == 119) system("pause");

    _quad_nodes[ip][0] = 0;  // remove all points connected to (ip)
}

inline void RemovePointFromCell_Point(size_t ip)
{
    size_t icell, ii, jj;


    ii = size_t((x[ip] - _xo) / _sx);
    jj = size_t((y[ip] - _yo) / _sy);
    icell = ii * _ny + jj;
    if (_cell_point2[icell][0] == 1) {
        _cell_point2[icell][0] = 0;
    } else {
        for (size_t V = 1; V < _cell_point2[icell][0]; V++) {
            if (_cell_point2[icell][V] == ip) {
                _cell_point2[icell][V] =
                    _cell_point2[icell][_cell_point2[icell][0]];
                break;
            }
        }
        _cell_point2[icell][0]--;
    }
}
inline void RemoveNodesConnectedWithTwoDelaunayEdges()
{  
    for (size_t i = 0; i < _num_points; i++) {
        if (!(x[i] >= 0.0 && x[i] < 1.0 && y[i] >= 0.0 && y[i] < 1.0)) {
            continue;
        }
      
        if (_quad_nodes[i][0] == 2) {
            if (false) {
                PlotQuadMesh();
                PlotSinglePoint(x[i], y[i]);
            }

            RemoveQuadConnectivity(i);
            RemovePointFromCell_Point(i);
            _num_fixed_quads++;

            if (x[i] >= 0.0 && x[i] < 1.0 && y[i] >= 0.0 && y[i] < 1.0) {
                if (_ghost_point[i][0] == 0) {
                    continue;
                }  // if no ghost
                for (size_t q = 1; q <= _ghost_point[i][0]; q++) {
                    RemoveConnectivity(_ghost_point[i][q]);
                    RemovePointFromCell_Point(_ghost_point[i][q]);
                    _num_fixed_quads++;
                }
            }

            if (false) {
                PlotQuadMesh();
                PlotSinglePoint(x[i], y[i]);
            }
        }
    }
}

inline void PostProcess()
{
    _list_of_removed_disks = new size_t[1000];

    while (true) {

        _num_fixed_quads = 0;
        FixingAngle();
        RemoveNodesConnectedWithTwoDelaunayEdges();
        if (_num_fixed_quads == 0) {
            break;
        }
    }
}