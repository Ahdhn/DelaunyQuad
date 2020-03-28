#pragma once

inline void FindMainAddGhosts(double xg1, double yg1, size_t ig1)
{
    // 1-find the main point of (xg1,yg1)
    // 2- add ghost points for this main point (except for ig1)

    // 2 |   4    | 7
    //___|________|____
    //   |  main  |
    // 1 | domain | 6
    //___|________|____
    // 0 |   3    | 5
    //   |        |
    // ig1 could be at one of the above 8 parts outside the domain
    size_t V, d, s, icell;
    size_t color = _disk_color[ig1];
    double _r(_r_input * 3.0), x_main, y_main;
    double r_squared = _r * _r;
    double limit(1.0 - _r);
    double x_p, y_p;
    double dx, dy;
    // 0
    if (xg1 < 0 && yg1 < 0) {
        x_main = xg1 + 1.0;
        y_main = yg1 + 1.0;
    }
    // 1
    if (xg1 < 0 && yg1 > 0 && yg1 < 1.0) {
        x_main = xg1 + 1.0;
        y_main = yg1;
    }
    // 2
    if (xg1 < 0 && yg1 > 1.0) {
        x_main = xg1 + 1.0;
        y_main = yg1 - 1.0;
    }
    // 3
    if (xg1 > 0 && xg1 < 1.0 && yg1 < 0.0) {
        x_main = xg1;
        y_main = yg1 + 1.0;
    }
    // 4
    if (xg1 > 0 && xg1 < 1.0 && yg1 > 1.0) {
        x_main = xg1;
        y_main = yg1 - 1.0;
    }
    // 5
    if (xg1 > 1.0 && yg1 < 0.0) {
        x_main = xg1 - 1.0;
        y_main = yg1 + 1.0;
    }
    // 6
    if (xg1 > 1.0 && yg1 > 0.0 && yg1 < 1.0) {
        x_main = xg1 - 1.0;
        y_main = yg1;
    }
    // 7
    if (xg1 > 1.0 && yg1 > 1.0) {
        x_main = xg1 - 1.0;
        y_main = yg1 - 1.0;
    }

    if (x_main < 0.0 || y_main < 0.0 || x_main > 1.0 || y_main > 1.0) {
        cout << "Error (0) at FindMainAddGhosts()... main point outside domain"
             << endl;
        system("pause");
    }

    size_t ii = size_t((x_main - _xo) / _sx);
    size_t jj = size_t((y_main - _yo) / _sy);
    if (false) {
        // check if the main point is covered
        size_t i_left = ii - 3;
        size_t i_right = ii + 3;
        size_t j_down = jj - 3;
        size_t j_up = jj + 3;
        if (i_left < 0) {
            i_left = 0;
        }
        if (i_right >= _nx) {
            i_right = _nx - 1;
        }
        if (j_down < 0) {
            j_down = 0;
        }
        if (j_up >= _ny) {
            j_up = _ny - 1;
        }

        size_t id;

        double r_out_sq(_r_b * _r_b), r_in_sq(_r_s * _r_s), r_com_sq;


        for (V = i_left; V <= i_right; V++) {
            for (d = j_down; d <= j_up; d++) {
                icell = V * _ny + d;
                if (_cell_point2[icell][0] == 0) {
                    continue;
                }
                for (s = 1; s <= _cell_point2[icell][0]; s++) {
                    id = _cell_point2[icell][s];
                    dx = x_main - x[id];
                    dy = y_main - y[id];
                    dx *= dx;
                    dy *= dy;
                    if (_disk_color[id] == color) {
                        r_com_sq = r_out_sq;
                    } else {
                        r_com_sq = r_in_sq;
                    }


                    if (dx + dy < r_com_sq) {
                        cout << "Error (1) at FindMainAddGhosts().. main point "
                                "covered"
                             << endl;
                        system("pause");
                    }
                }
            }
        }
    }
    // add the main point
    x[_num_points] = x_main;
    y[_num_points] = y_main;
    _disk_color[_num_points] = color;
    _num_points++;
    icell = ii * _ny + jj;
    _cell_point2[icell][0]++;
    if (_cell_point2[icell][0] > 4) {
        cout << "Error (2) at FindMainAddGhosts .....more than 4 points in a "
                "cell.."
             << endl;
        system("pause");
    }
    size_t main_point = _num_points - 1;
    _cell_point2[icell][_cell_point2[icell][0]] = main_point;


    _ghost_point[ig1][0] = main_point;  // add the main point to the ghost point
    _ghost_point[main_point][0] = 1;    // add ghost point to the main point
    _ghost_point[main_point][_ghost_point[main_point][0]] = ig1;


    // find the other ghost points of this main point

    if (!(x_main < _r || x_main > limit || y_main < _r || y_main > limit)) {
        cout << "Error (3) at FindMainAddGhosts .....something is wrong"
             << endl;
        system("pause");
    }

    bool there;

    if (x_main < _r) {
        y_p = y_main;
        x_p = 1.0 + x_main;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {

            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;
            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;
            _disk_color[_num_points] = color;
            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (0) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    if (x_main > limit) {
        y_p = y_main;
        x_p = x_main - 1.0;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (1) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    if (y_main < _r) {
        x_p = x_main;
        y_p = 1.0 + y_main;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;
            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (2) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    if (y_main > limit) {
        x_p = x_main;
        y_p = y_main - 1.0;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (3) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    double t1(x_main - 1.0), t2(y_main - 1.0), yy_sq(y_main * y_main),
        xx_sq(x_main * x_main);
    t1 *= t1;
    t2 *= t2;

    if ((yy_sq + xx_sq) < r_squared) {
        x_p = 1.0 + x_main;
        y_p = 1.0 + y_main;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (4) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }


    if ((t1 + yy_sq) < r_squared) {
        x_p = x_main - 1.0;
        y_p = 1.0 + y_main;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);

            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (5) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    if ((xx_sq + t2) < r_squared) {
        x_p = 1.0 + x_main;
        y_p = y_main - 1.0;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {

            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (6) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }

    if ((t1 + t2) < r_squared) {
        x_p = x_main - 1.0;
        y_p = y_main - 1.0;

        there = false;
        dx = x_p - xg1;
        dy = y_p - yg1;
        dx *= dx;
        dy *= dy;
        if (dx + dy < 10E-5) {
            there = true;
        }

        if (!there) {
            ii = size_t(abs(_xo - x_p) / _sx);
            jj = size_t(abs(_yo - y_p) / _sy);
            icell = ii * _ny + jj;
            x[_num_points] = x_p;
            y[_num_points] = y_p;

            _ghost_point[_num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = _num_points;

            _disk_color[_num_points] = color;

            _cell_point2[icell][0]++;
            _cell_point2[icell][_cell_point2[icell][0]] = _num_points;
            _num_points++;
            if (_cell_point2[icell][0] > 10) {
                cout << "Error (8) at FindMainAddGhosts().." << endl;
                system("pause");
            }
        }
    }
}
inline void ReconstructDelaunay(size_t ipoint)
{
    size_t i, j, icell, ip;
    double xp, yp;
    int    c, r;
    double node_x, node_y;

    node_x = x[ipoint];
    node_y = y[ipoint];

    i = size_t((node_x - _xo) / _sx);
    j = size_t((node_y - _yo) / _sy);

    for (c = int(i - 10); c <= int(i + 10); c++) {
        for (r = int(j - 10); r <= int(j + 10); r++) {
            if (c >= _nx || r >= _ny || c < 0 || r < 0) {
                continue;
            }
            // if(c==i&&r==j){continue;}

            icell = c * _ny + r;

            if (_cell_point2[icell][0] == 0) {
                continue;
            }

            for (size_t AM = 1; AM <= _cell_point2[icell][0]; AM++) {
                ip = _cell_point2[icell][AM];
                if (ip == ipoint) {
                    continue;
                }
                xp = x[ip];
                yp = y[ip];

                if (_active[ip] == 0) {
                    cout << "Error at ReconstructDelaunay().. inactive point "
                            "in a cell "
                         << endl;
                    system("pause");
                }

                if (CheckConnectivity(ipoint, ip)) {
                    continue;
                }
                if (DelaunayEdgeExist(ipoint, ip, i, j, c, r, node_x, node_y,
                                      xp, yp)) {

                    _node_nodes[ipoint][0]++;
                    _node_nodes[ipoint][_node_nodes[ipoint][0]] = ip;
                    _node_nodes[ip][0]++;
                    _node_nodes[ip][_node_nodes[ip][0]] = ipoint;
                    if (_node_nodes[ipoint][0] > 19 || _node_nodes[ip][0] > 19) {
                        cout << "Error.. way to many points connected to point "
                                "in delaunay "
                             << endl;
                        system("pause");
                        PlotMesh();
                        PlotSinglePoint(x[ipoint], y[ipoint]);
                    }

                    // PlotMesh();
                }
            }
        }
    }
}
inline void RemoveConnectivity(size_t ip)
{
    // remove (ip) from all points connected to it
    // remove all points connected to (ip)
    size_t V, d, ip1;
    bool   found;
    for (V = 1; V <= _node_nodes[ip][0]; V++) {
        ip1 = _node_nodes[ip][V];
        // remove ip from ip1
        found = false;
        for (d = 1; d < _node_nodes[ip1][0]; d++) {
            if (_node_nodes[ip1][d] == ip) {
                found = true;
            }  // first find ip
            if (found) {
                _node_nodes[ip1][d] =
                    _node_nodes[ip1][d + 1];  // move the point on step to
                                             // preserve the sorted order
            }
        }
        if (!found && _node_nodes[ip1][d] != ip) {
            cout << "Error at RemoveConnectivity().." << endl;
            system("pause");
        }
        _node_nodes[ip1][0]--;
    }


    _node_nodes[ip][0] = 0;  // remove all points connected to (ip)
}
inline double GetCircumcentre(size_t  ip1,
                              size_t  ip2,
                              size_t  ip3,
                              double& xi,
                              double& yi)
{
    double x1 = 0.5 * (x[ip1] + x[ip2]);
    double y1 = 0.5 * (y[ip1] + y[ip2]);

    double dx1 = x[ip1] - x[ip2];
    double dy1 = y[ip1] - y[ip2];

    double x2 = 0.5 * (x[ip1] + x[ip3]);
    double y2 = 0.5 * (y[ip1] + y[ip3]);

    double dx2 = x[ip1] - x[ip3];
    double dy2 = y[ip1] - y[ip3];

    double dx, dy;

    if (dy1 == 0 && dy2 != 0) {
        xi = x1;
        yi = (x2 - xi) * dx2 / dy2 + y2;

        dx = xi - x[ip1];
        dy = yi - y[ip1];

        return (dx * dx + dy * dy);
    } else if (dy2 == 0 && dy1 != 0) {
        xi = x2;
        yi = (x1 - xi) * dx1 / dy1 + y1;

        dx = xi - x[ip1];
        dy = yi - y[ip1];

        return (dx * dx + dy * dy);
    }

    else if (dy1 == 0 && dy2 == 0) {
        cout << " Error ..Two Horizontal line " << endl;
        system("pause");
    }


    double m1 = -1 * dx1 / dy1;
    double m2 = -1 * dx2 / dy2;

    xi = (y2 - y1 + m1 * x1 - m2 * x2) / (m1 - m2);
    yi = (x1 - xi) * dx1 / dy1 + y1;

    dx = xi - x[ip1];
    dy = yi - y[ip1];

    return (dx * dx + dy * dy);
}

inline void FixDelaunay()
{
    double xc, yc, x_star, y_star, rr;
    double rr_sqrt;  // extremities of the grid
   
    size_t i, j, ip1, ip2, jnext, color, icell, ii, jj, s, id, q, w, ig;
      
    size_t V, d;

    size_t* r_b_covered =
        new size_t[2];  // array tells if the cell is covered by r_out of only
                        // red,only blue or both red and blue disks. if
                        // r_b_covered[0]=1 -->means covered by a blue disk
                        // other wise not covered by any blue if
                        // r_b_covered[1]=1 -->means covered by a red disk other
                        // wise not covered by any red

    double r_out_sq(_r_b * _r_b), r_in_sq(_r_s * _r_s);
    int    i_left, i_right, j_down, j_up;
    bool   again = true;

    while (again) {
        again = false;

        for (i = 0; i < _num_points; i++) {

            x_star = x[i];
            y_star = y[i];

            if (x_star < 0.0 || x_star > 1.0 || y_star < 0.0 || y_star > 1.0) {
                continue;
            }

            for (j = 1; j <= _node_nodes[i][0]; j++) {
                jnext = j + 1;
                if (j == _node_nodes[i][0]) {
                    jnext = 1;
                }

                ip1 = _node_nodes[i][j];
                ip2 = _node_nodes[i][jnext];

                if (_disk_color[ip1] == _disk_color[ip2] &&
                    _disk_color[i] == _disk_color[ip2]) {  // if three vertices
                                                           // with same color


                    if (_active[i] == 0 || _active[ip1] == 0 ||
                        _active[ip2] == 0) {
                        cout << "Error.. not active point in delaunay " << endl;
                        system("pause");
                    }


                    color = _disk_color[i];  // the color

                    rr = GetCircumcentre(i, ip1, ip2, xc, yc);
                    rr_sqrt = rr;
                    rr = sqrt(rr);


                    // if the circumcenter is outside, It must have a refelction
                    // that is inside the domain
                    // when looping we will reach it. (so that when added, we
                    // can add it mirror points correctly)

                    if (xc < 0.0 || xc >= 1.0 || yc < 0.0 || yc >= 1.0) {
                        continue;
                    }
                    again = true;

                    if (i >= _s_d && ip1 >= _s_d && ip2 >= _s_d) {
                        cout << "Three new inserted same color" << endl;
                        system("pause");
                    }

                    // cout<<"i= "<<i<<endl;

                    if (false) {
                        PlotMesh();
                        PlotSinglePoint(xc, yc);
                    }


                    size_t KO;
                    if (_list_of_removed_disks[0] >= 1) {
                        size_t baba =
                            _list_of_removed_disks[_list_of_removed_disks[0]];
                        // ins++;
                        x[baba] = xc;
                        y[baba] = yc;
                        if (color == 1) {
                            _disk_color[baba] = 0;
                        } else {
                            _disk_color[baba] = 1;
                        }
                        id = baba;  // lets just name the new point id

                        _circumcenter[_list_of_removed_disks
                                          [_list_of_removed_disks[0]]] = true;
                        KO = baba;
                        _num_points++;

                        ii = size_t((xc - _xo) / _sx);
                        jj = size_t((yc - _yo) / _sy);
                        icell = ii * _ny + jj;
                        _cell_point2[icell][0]++;
                        _cell_point2[icell][_cell_point2[icell][0]] = baba;

                        if (xc >= 0.0 && xc < 1.0 && yc >= 0.0 && yc < 1.0) {
                            MirrorPoint(baba, xc, yc, _r_input,
                                        _disk_color[baba]);
                        } else {
                            FindMainAddGhosts(xc, yc, _num_points - 1);
                        }
                        _list_of_removed_disks[0]--;
                    }

                    else {
                        x[_num_points] = xc;
                        y[_num_points] = yc;
                        if (color == 1) {
                            _disk_color[_num_points] = 0;
                        } else {
                            _disk_color[_num_points] = 1;
                        }
                        id = _num_points;  // lets just name the new point id

                        _circumcenter[_num_points] = true;
                        KO = _num_points;
                        _num_points++;

                        ii = size_t((xc - _xo) / _sx);
                        jj = size_t((yc - _yo) / _sy);
                        icell = ii * _ny + jj;
                        _cell_point2[icell][0]++;
                        _cell_point2[icell][_cell_point2[icell][0]] =
                            _num_points - 1;
                        if (xc >= 0.0 && xc < 1.0 && yc >= 0.0 && yc < 1.0) {
                            MirrorPoint(
                                _num_points - 1, xc, yc, _r_input,
                                _disk_color[_num_points -
                                            1]);  // add new mirror points
                        } else {
                            FindMainAddGhosts(xc, yc, _num_points - 1);
                        }
                    }

                    if (_ghost_point[KO][0] > 0 && xc >= 0.0 && xc < 1.0 &&
                        yc >= 0.0 && yc < 1.0) {
                        for (V = 1; V <= _ghost_point[KO][0]; V++) {
                            _circumcenter[_ghost_point[KO][V]] = true;
                        }
                    }


                    // DestroyDelaunayAndRebuild
                    // update the 6 layers around the new point
                    ii = size_t((xc - _xo) / _sx);
                    jj = size_t((yc - _yo) / _sy);

                    // update the 6 layers around the new point
                    i_left = ii - 10;
                    i_right = ii + 10;
                    j_down = jj - 10;
                    j_up = jj + 10;

                    if (i_left < 0) {
                        i_left = 0;
                    }
                    if (i_right >= _nx) {
                        i_right = _nx - 1;
                    }
                    if (j_down < 0) {
                        j_down = 0;
                    }
                    if (j_up >= _ny) {
                        j_up = _ny - 1;
                    }


                    // remove delaunay connectivity of the points in the 6
                    // layers
                    for (V = i_left; V <= i_right; V++) {
                        for (d = j_down; d <= j_up; d++) {
                            icell = V * _ny + d;
                            if (_cell_point2[icell][0] == 0) {
                                continue;
                            }
                            for (s = 1; s <= _cell_point2[icell][0]; s++) {
                                id = _cell_point2[icell][s];

                                if (_active[id] == 0) {
                                    cout << "Error(1).. not active point in "
                                            "delaunay "
                                         << endl;
                                    system("pause");
                                }

                                RemoveConnectivity(
                                    id);  // remove the connectivity of the main
                                          // point


                                // if (id) is inside the domain
                                if (x[id] >= 0.0 && x[id] < 1.0 &&
                                    y[id] >= 0.0 && y[id] < 1.0) {
                                    if (_ghost_point[id][0] == 0) {
                                        continue;
                                    }  // if no ghost points
                                    for (q = 1; q <= _ghost_point[id][0];
                                         q++) {  // remove the connetivity of
                                                 // the ghost point associated
                                                 // with this main point
                                        if (_active[_ghost_point[id][q]] == 0) {
                                            cout << "Error(2).. not active "
                                                    "point in delaunay "
                                                 << endl;
                                            system("pause");
                                        }
                                        RemoveConnectivity(_ghost_point[id][q]);
                                    }
                                } else {  // if it was a ghost point
                                    // then remove the connectivity of its main
                                    // point
                                    if (_active[_ghost_point[id][0]] == 0) {
                                        cout << "Error(3).. not active point "
                                                "in delaunay "
                                             << endl;
                                        system("pause");
                                    }
                                    RemoveConnectivity(_ghost_point[id][0]);
                                }
                            }
                        }
                    }

                    // PlotMesh();

                    // reconstruct delaunay of the points in the 6 layers again
                    for (V = i_left; V <= i_right; V++) {
                        for (d = j_down; d <= j_up; d++) {
                            icell = V * _ny + d;
                            if (_cell_point2[icell][0] == 0) {
                                continue;
                            }
                            for (s = 1; s <= _cell_point2[icell][0]; s++) {
                                id = _cell_point2[icell][s];
                                if (_active[id] == 0) {
                                    cout << "Error(0).. not active point in a "
                                            "cell "
                                         << endl;
                                    system("pause");
                                }

                                ReconstructDelaunay(
                                    id);  // reconstruct delaunay of this main
                                          // point

                                // if (id) is inside the domain
                                if (x[id] >= 0.0 && x[id] < 1.0 &&
                                    y[id] >= 0.0 && y[id] < 1.0) {
                                    if (_ghost_point[id][0] == 0) {
                                        continue;
                                    }  // if no ghost points
                                    for (q = 1; q <= _ghost_point[id][0];
                                         q++) {  // reconstruct delaunay of the
                                                 // ghost point associated with
                                                 // this main point
                                        if (_active[_ghost_point[id][q]] == 0) {
                                            cout << "Error(1).. not active "
                                                    "point in a cell "
                                                 << endl;
                                            system("pause");
                                        }

                                        ReconstructDelaunay(
                                            _ghost_point[id][q]);
                                    }
                                } else {  // if it was a ghost point
                                    // then reconstruct delaunay of its main
                                    // point
                                    if (_active[_ghost_point[id][0]] == 0) {
                                        cout << "Error(2).. not active point "
                                                "in a cell "
                                             << endl;
                                        system("pause");
                                    }

                                    ReconstructDelaunay(_ghost_point[id][0]);
                                }
                            }
                        }
                    }


                    // resort delaunay connectivity for the points in the 6
                    // layers
                    for (V = i_left; V <= i_right; V++) {
                        for (d = j_down; d <= j_up; d++) {
                            icell = V * _ny + d;
                            if (_cell_point2[icell][0] == 0) {
                                continue;
                            }
                            for (s = 1; s <= _cell_point2[icell][0]; s++) {
                                SortDelaunay2(
                                    _cell_point2[icell]
                                               [s]);  // sorting delaunay edges
                                                      // of this main point

                                id = _cell_point2[icell]
                                                [s];  // sort the delaunay edges
                                                      // connected to this main
                                                      // point
                                for (q = 1; q <= _node_nodes[id][0]; q++) {
                                    SortDelaunay2(_node_nodes[id][q]);
                                }

                                // if it's insdie the domain
                                if (x[id] >= 0.0 && x[id] < 1.0 &&
                                    y[id] >= 0.0 && y[id] < 1.0) {
                                    if (_ghost_point[id][0] == 0) {
                                        continue;
                                    }  // if no ghost points
                                    for (q = 1; q <= _ghost_point[id][0];
                                         q++) {  //  sort delaunay edges of this
                                                 //  ghost point (associated
                                                 //  with this main point)
                                        ig = _ghost_point[id][q];
                                        SortDelaunay2(ig);
                                        for (w = 1; w <= _node_nodes[ig][0];
                                             w++) {
                                            SortDelaunay2(
                                                _node_nodes[ig]
                                                          [w]);  // sort edge
                                                                 // connected to
                                                                 // this ghost
                                                                 // point
                                        }
                                    }
                                } else {  // if it's a ghost
                                    ig = _ghost_point[id][0];
                                    SortDelaunay2(
                                        ig);  // sort edges of the main point of
                                              // this ghost point

                                    for (q = 1; q <= _node_nodes[ig][0]; q++) {
                                        SortDelaunay2(
                                            _node_nodes[ig]
                                                      [q]);  // then sort the
                                                             // edges connect to
                                                             // the main
                                    }
                                }
                            }
                        }
                    }

                    if (false) {
                        PlotMesh();
                        PLOT(_r_s, _r_b, 0, 0, 0);
                    }

                    break;
                }
            }
        }
    }
}