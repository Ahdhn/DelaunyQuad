#pragma once

inline void FarPoint(double  xc,
                     double  yc,
                     double  xmin,
                     double  xmax,
                     double  ymin,
                     double  ymax,
                     double& xfar,
                     double& yfar)
{
    if (abs(xc - xmin) > abs(xc - xmax)) {
        xfar = xmin;
    } else {
        xfar = xmax;
    }
    if (abs(yc - ymin) > abs(yc - ymax)) {
        yfar = ymin;
    } else {
        yfar = ymax;
    }
}

inline void MirrorPoint(size_t main_point,
                        double xx,
                        double yy,
                        double _r,
                        size_t disk_color)
{

    _r *= 5.0;
    double r_squared = _r * _r;
    double x_p, y_p;
    double limit = 1.0 - _r;

    if (xx < _r || xx > limit || yy < _r || yy > limit) {

        size_t i, j, icell;
        if (xx < _r) {
            y_p = yy;
            x_p = 1.0 + xx;
        }

        if (xx < _r) {
            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;

            _disk_color[num_points] = disk_color;

            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (0) at MirrorPoint().." << endl;
                if (false) {
                    PlotSinglePoint(x_p, y_p);
                }
                system("pause");
            }
        }
        if (xx > limit) {
            y_p = yy;
            x_p = xx - 1.0;
        }

        if (xx > limit) {
            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (1) at MirrorPoint().." << endl;
                system("pause");
            }
        }
        if (yy < _r) {
            x_p = xx;
            y_p = 1.0 + yy;
        }

        if (yy < _r) {
            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (2) at MirrorPoint().." << endl;
                system("pause");
            }
        }
        if (yy > limit) {
            x_p = xx;
            y_p = yy - 1.0;
        }
        if (yy > limit) {
            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            if (false) {
                PlotSinglePoint(x_p, y_p);
            }

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (3) at MirrorPoint().." << endl;
                system("pause");
            }
        }
        size_t in(0);
        double t1(xx - 1.0), t2(yy - 1.0), yy_sq(yy * yy), xx_sq(xx * xx);
        t1 *= t1;
        t2 *= t2;
        if ((yy_sq + xx_sq) < r_squared) {
            x_p = 1.0 + xx;
            y_p = 1.0 + yy;
            in++;
        }

        if (in > 0) {
            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (4) at MirrorPoint().." << endl;
                system("pause");
            }
        }

        in = 0;

        if ((t1 + yy_sq) < r_squared) {
            x_p = xx - 1.0;
            y_p = 1.0 + yy;
            in++;
        }
        if (in > 0) {

            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;

            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (5) at MirrorPoint().." << endl;
                system("pause");
            }
        }

        in = 0;
        if ((xx_sq + t2) < r_squared) {
            x_p = 1.0 + xx;
            y_p = yy - 1.0;
            in++;
        }

        if (in > 0) {

            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (6) at MirrorPoint().." << endl;
                system("pause");
            }
        }

        in = 0;

        if ((t1 + t2) < r_squared) {
            x_p = xx - 1.0;
            y_p = yy - 1.0;
            in++;
        }

        if (in > 0) {

            i = size_t(abs(_xo - x_p) / _sx);
            j = size_t(abs(_yo - y_p) / _sy);
            icell = i * _ny + j;
            x[num_points] = x_p;
            y[num_points] = y_p;

            _ghost_point[num_points][0] =
                main_point;  // add the main point to the ghost point
            _ghost_point[main_point][0]++;  // add ghost point to the main point
            _ghost_point[main_point][_ghost_point[main_point][0]] = num_points;
            _disk_color[num_points] = disk_color;
            // num_points++;
            // cell_point[icell] = num_points;
            cell_point2[icell][0]++;
            cell_point2[icell][cell_point2[icell][0]] = num_points;
            _circumcenter[num_points] = false;
            num_points++;
            if (cell_point2[icell][0] > 2) {
                cout << "Error (8) at MirrorPoint().." << endl;
                system("pause");
            }
        }
    }
    _r /= 5.0;
    r_squared = _r * _r;
}

inline void DartThrowingWithColors(
    double _r,
    double r_in_sq,
    double r_out_sq) 
{

    size_t num_darts, idart, rand_index, icell, neighbor_cell, pointer, count;
    size_t c, r, i, ii, iii, j, jj, jjj, half_cell_x, half_cell_y, iactive;
    double rand_x, rand_y;
    double x_start, x_end, y_start, y_end, xx, yy;
    double xmin, ymin, xmax, ymax;
    double xp, yp, dist_x, dist_y, dist;
    double sssx, sssy;
    double xfar, yfar;
    double RF;

    size_t num_tmp_active;
    bool   check_point;
    bool   no_active_cells;
    bool   inside;
    double r_com_sq;

    double ssx = _sx;
    double ssy = _sy;

    size_t* r_b_covered =
        new size_t[2];  // array tells if the cell is covered by r_out of only
                        // red,only blue or both red and blue disks. if
                        // r_b_covered[0]=1 -->means covered by a blue disk
                        // other wise not covered by any blue if
                        // r_b_covered[1]=1 -->means covered by a red disk other
                        // wise not covered by any red

    size_t color;  // 0 blue 1 red.. used to decide the dart color randomly


    RF = 0.8;

    for (int iref = 0; iref < num_refinement_levels;
         iref++)  // loop for refining active cells
    {

        no_active_cells = false;
        num_darts =
            size_t(RF * num_active);  // number of the darts to be thrown


        // PLOT(_r_s,_r_b,iref+1,ssx,ssy);

        for (idart = 0; idart < num_darts; idart++) {


            rand_index = size_t((num_active - 1) * ((double)rand() / RAND_MAX));

            ii = active_i[rand_index];
            jj = active_j[rand_index];

            i = size_t(double(ii / pow(2, double(iref))));
            j = size_t(double(jj / pow(2, double(iref))));

            icell = i * _ny + j;

            rand_x = ((double)rand() / RAND_MAX);
            rand_y = ((double)rand() / RAND_MAX);

            if (rand_x == 1 || rand_y == 1) {
                continue;
            }

            xx = _xo + (ii + rand_x) * (ssx);
            yy = _yo + (jj + rand_y) * (ssy);

            if (xx > 1 || xx < 0 || yy > 1 || yy < 0)
                continue;

            if (rand_x > 0.5) {
                color = 0;
            } else {
                color = 1;
            }
            /* if(num_points==0){color=0;}
             else if(_disk_color[num_points-1]==0){color=1;}
             else if(_disk_color[num_points-1]==1){color=0;}*/

            // checking the conflict
            bool close = false;
            for (c = i - 4; c <= i + 4; c++) {
                for (r = j - 4; r <= j + 4; r++) {
                    if (c >= _nx || r >= _ny || c < 0 || r < 0) {
                        continue;
                    }

                    neighbor_cell = c * _ny + r;
                    // pointer = cell_point[neighbor_cell];
                    // if (pointer > 0){
                    if (cell_point2[neighbor_cell][0] > 1) {
                        cout << "Error(0) at DartThrowingWithColors(). "
                             << endl;
                        system("pause");
                    }
                    if (cell_point2[neighbor_cell][0] > 0) {
                        // pointer-=1;
                        pointer = cell_point2[neighbor_cell]
                                             [1];  // at this stage each cell
                                                   // contains one point only
                        xp = x[pointer];
                        yp = y[pointer];

                        // rand_x = ((double)rand() / RAND_MAX);

                        dist = (xp - xx) * (xp - xx) + (yp - yy) * (yp - yy);

                        if (color == _disk_color[pointer]) {
                            r_com_sq = r_out_sq;
                        }
                        if (color != _disk_color[pointer]) {
                            r_com_sq = r_in_sq;
                        }
                        if (dist < r_com_sq) {
                        
                            close = true;
                            break;
                        }
                    }
                }
                if (close) {
                    break;
                }
            }

            if (!close) {

                // cell_point[icell] = num_points + 1;  // adding the new dart
                // and updating x,y arraies
                cell_point2[icell][0]++;  // adding the new dart and updating
                                          // x,y arraies
                cell_point2[icell][cell_point2[icell][0]] = num_points;
                if (cell_point2[icell][0] > 1) {
                    cout << "Error(1) at DartThrowingWithColors(). " << endl;
                    system("pause");
                }


                x[num_points] = xx;
                y[num_points] = yy;

                _disk_color[num_points] = color;
                /*if(num_points==0){_disk_color[num_points]=0;}// first disk is
                blue else{

                    if(_disk_color[num_points-1]==0){_disk_color[num_points]=1;}//
                if previous disk is blue, then this one is red
                    else{_disk_color[num_points]=0;} // esle, make it blue
                }*/

                num_points++;


                MirrorPoint(num_points - 1, xx, yy, _r_input,
                            _disk_color[num_points - 1]);

                // PLOT(_r_s,_r_b,0,0,0);

                num_active--;

                if (num_active == 0) {
                    no_active_cells = true;
                    break;
                }

                active_i[rand_index] = active_i[num_active];
                active_j[rand_index] = active_j[num_active];
            }        
        }


        // PLOT(_r_s,_r_b,0,0,0);


        // PLOT(_r_s,_r_b,iref+1,ssx,ssy);

        if (no_active_cells == true) {
            // delete[] active_i;
            // delete[] active_j;
            break;
        }

        // tmp_active_i = new size_t[4*num_active];
        // tmp_active_j = new size_t[4*num_active];
        num_tmp_active = 0;
        no_active_cells = false;

        sssx = ssx / 2.0;
        sssy = ssy / 2.0;

        for (iactive = 0; iactive < num_active; iactive++) {

            ii = active_i[iactive];
            jj = active_j[iactive];


            i = size_t(double(ii / pow(2.0, double(iref))));
            j = size_t(double(jj / pow(2.0, double(iref))));

            xmin = _xo + ii * ssx;
            ymin = _yo + jj * ssy;
            xmax = xmin + ssx;
            ymax = ymin + ssy;

            check_point = true;

            r_b_covered[0] = 0;
            r_b_covered[1] = 0;

            for (c = i - 4; c <= i + 4; c++) {
                for (r = j - 4; r <= j + 4; r++) {
                    if (c >= _nx || r >= _ny || c < 0 || r < 0)
                        continue;
                    neighbor_cell = c * _ny + r;
                    // pointer = cell_point[neighbor_cell];

                    // if(pointer>0){
                    if (cell_point2[neighbor_cell][0] > 0) {
                        // pointer--;
                        pointer = cell_point2[neighbor_cell]
                                             [1];  // at this stage cell
                                                   // contains one point only
                        xp = x[pointer];
                        yp = y[pointer];
                        FarPoint(xp, yp, xmin, xmax, ymin, ymax, xfar, yfar);
                        dist_x = (xfar - xp);
                        dist_y = (yfar - yp);

                        dist = (dist_x * dist_x) + (dist_y * dist_y);
                        if (dist < r_in_sq) {
                            check_point = false;
                        }

                        if (dist < r_out_sq) {
                            if (_disk_color[pointer] == 0) {
                                r_b_covered[0] = 1;
                            } else {
                                r_b_covered[1] = 1;
                            }
                        }
                    }
                }
                if (check_point == false) {
                    break;
                }
            }
            if (check_point == false) {
                continue;
            }
            if (r_b_covered[0] == 1 && r_b_covered[1] == 1) {
                continue;
            }  // covered by a blue and red disk


            for (half_cell_x = ii * 2; half_cell_x <= (ii * 2 + 1);
                 half_cell_x++) {
                iii = half_cell_x;
                for (half_cell_y = jj * 2; half_cell_y <= (jj * 2 + 1);
                     half_cell_y++) {
                    jjj = half_cell_y;

                    xmin = _xo + iii * sssx;
                    ymin = _yo + jjj * sssy;
                    xmax = xmin + sssx;
                    ymax = ymin + sssy;
                    if (xmin >= 1 || ymin >= 1 || xmax <= 0 || ymax <= 0) {
                        continue;
                    }

                    check_point = true;

                    r_b_covered[0] = 0;
                    r_b_covered[1] = 0;

                    for (c = i - 4; c <= i + 4; c++) {
                        for (r = j - 4; r <= j + 4; r++) {
                            if (c >= _nx || r >= _ny || c < 0 || r < 0) {
                                continue;
                            }
                            neighbor_cell = c * _ny + r;
                            // pointer = cell_point[neighbor_cell];

                            // if(pointer != 0){
                            if (cell_point2[neighbor_cell][0] > 0) {
                                // pointer--;
                                pointer = cell_point2[neighbor_cell][1];
                                xp = x[pointer];
                                yp = y[pointer];
                                FarPoint(xp, yp, xmin, xmax, ymin, ymax, xfar,
                                         yfar);
                                dist_x = (xfar - xp);
                                dist_y = (yfar - yp);
                                dist = (dist_x * dist_x) + (dist_y * dist_y);
                                if (dist < r_in_sq) {
                                    check_point = false;
                                }

                                if (dist < r_out_sq) {
                                    if (_disk_color[pointer] == 0) {
                                        r_b_covered[0] = 1;
                                    } else {
                                        r_b_covered[1] = 1;
                                    }
                                }
                            }
                        }
                        if (check_point == false) {
                            break;
                        }
                    }
                    if (check_point == false) {
                        continue;
                    }
                    if (r_b_covered[0] == 1 && r_b_covered[1] == 1) {
                        continue;
                    }  // covered by a blue and red disk

                    tmp_active_i[num_tmp_active] = iii;
                    tmp_active_j[num_tmp_active] = jjj;
                    num_tmp_active++;
                    if (num_tmp_active > num_expected) {
                        cout << "Do something (0).. We are losing the baby "
                             << endl;
                        system("pause");
                    }
                }
            }
        }


        // delete [] active_i;
        // delete [] active_j;

        num_active = num_tmp_active;
        if (num_active != 0) {

            // active_i = new size_t [num_active];
            // active_j = new size_t [num_active];

            for (count = 0; count < num_active; count++) {
                active_i[count] = tmp_active_i[count];
                active_j[count] = tmp_active_j[count];
            }
            // delete[] tmp_active_i;
            // delete[] tmp_active_j;
        } else {
            break;
        }
        ssx /= 2.0;
        ssy /= 2.0;
    }
}