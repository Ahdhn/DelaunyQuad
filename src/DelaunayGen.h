#pragma once

inline size_t cell_index(size_t i, size_t j)
{
    size_t icell = i * _ny + j;
    return icell;
}

inline bool CheckConnectivity(size_t index_A, size_t index_B)
{
    bool check = false;
    for (size_t i = 1; i <= _node_nodes[index_A][0]; i++) {
        if (_node_nodes[index_A][i] == index_B) {
            check = true;
            break;
        }
    }
    return check;
}

inline double GetAngle360(double dy, double dx)
{
    double theta = atan2(dy, dx);
    if (theta < 0)
        theta += TwoPI;
    //	return theta*180/PI;
    return theta;
}
inline void Partition_2(int& i, int& j, size_t ip)
{
    int    pivot = (i + j) / 2;
    size_t tmp;
    double tmp1;
    double pivot_value = _angles[pivot];

    while (i <= j) {
        while (pivot_value < _angles[j]) {
            j--;
        }
        while (pivot_value > _angles[i]) {
            i++;
        }
        if (i <= j) {
            tmp = _node_nodes[ip][i + 1];
            _node_nodes[ip][i + 1] = _node_nodes[ip][j + 1];
            _node_nodes[ip][j + 1] = tmp;
            tmp1 = _angles[i];
            _angles[i] = _angles[j];
            _angles[j] = tmp1;
            j--;
            i++;
        }
    }
}
inline void Sort_2(int left, int right, size_t ip)
{
    int i, j;
    i = left;
    j = right;
    Partition_2(i, j, ip);
    if (left < j) {
        Sort_2(left, j, ip);
    }
    if (i < right) {
        Sort_2(i, right, ip);
    }
}
inline void SortDelaunay2(size_t ip)
{
    size_t num_nodes = _node_nodes[ip][0];
    double dx, dy;
    size_t ip1;

    double xp = x[ip];
    double yp = y[ip];

    for (size_t i = 1; i <= num_nodes; i++) {

        ip1 = _node_nodes[ip][i];
        dx = x[ip1] - xp;
        dy = y[ip1] - yp;

        _angles[i - 1] = GetAngle360(dy, dx);
    }
    Sort_2(0, num_nodes - 1, ip);
}

inline void Trimming(double  xi,
                     double  yi,
                     double& x1,
                     double& y1,
                     double& x2,
                     double& y2,
                     double  x3,
                     double  y3,
                     double  x4,
                     double  y4)
{
    double n_x = -(y4 - y3);
    double n_y = (x4 - x3);


    double dot_1 = n_x * (x3 - x1) + n_y * (y3 - y1);

    double dot_2 = n_x * (x2 - x1) + n_y * (y2 - y1);


    double u = dot_1 / dot_2;

    if (u > 0 && u < 1)

    {
        double x = x1 + u * (x2 - x1);  // intersection points
        double y = y1 + u * (y2 - y1);

        if (dot_1 < 0) {
            x1 = x;
            y1 = y;
        }

        else {
            x2 = x;
            y2 = y;
        }


    }

    else if (dot_1 < 0)  // the degenerated edge case
    {
        x1 = xi;
        x2 = xi;
        y1 = yi;
        y2 = yi;
    }
}

inline void VerticalLine(double  x1,
                         double  y1,
                         double  x2,
                         double  y2,
                         double& qo_x,
                         double& qo_y,
                         double& qn_x,
                         double& qn_y)
{
    double x_mid = (x1 + x2) / 2.0;
    double y_mid = (y1 + y2) / 2.0;

    double dy = -(x2 - x1);  // the slope of the horizontal line is equal ( -1 /
                             // slope of vertical line )
    double dx = (y2 - y1);

    double t = sqrt(pow(dy, 2) + pow(dx, 2));  // Normalize the vector
    dy = dy / t;
    dx = dx / t;
    /*
    if(dy < 10E-6 && dy > -10E-6)
    {
        dx = dx * -1;
    }
    if(dx < 10E-6 && dx > -10E-6)
    {
        dy = dy * -1;
    }
    */
    qo_x = x_mid - 6 * _r_input * dx;
    qo_y = y_mid - 6 * _r_input * dy;
    qn_x = x_mid + 6 * _r_input * dx;
    qn_y = y_mid + 6 * _r_input * dy;
}

inline bool DelaunayEdgeExist(size_t point1,
                              size_t point2,
                              size_t i_node_i,
                              size_t j_node_i,
                              size_t i_node_j,
                              size_t j_node_j,
                              double xi,
                              double yi,
                              double xj,
                              double yj)
{
    double eps = _r_input * 1E-11;

    double po_x, po_y, pn_x,
        pn_y;  // start and end points of the initial voronoi edge

    double qo_x, qo_y, qn_x,
        qn_y;  // end and start poitns of the trimming voronoi line

    double pcell_x, pcell_y;

    VerticalLine(xi, yi, xj, yj, po_x, po_y, pn_x,
                 pn_y);  // construct initial voronoi edge

    double m = (pn_y - po_y) / (pn_x - po_x);

    bool   outside_edge = false;
    size_t ip;

    for (int c = int(i_node_i - 15); c <= int(i_node_i + 15); c++) {
        for (int r = int(j_node_i - 15); r <= int(j_node_i + 15); r++) {
            outside_edge = false;

            if (c >= _nx || r >= _ny || c < 0 || r < 0) {
                continue;
            }

            size_t pcell = cell_index(c, r);  // getting the cell index

            // ip = cell_point[pcell];
            if (_cell_point2[pcell][0] == 0) {
                continue;
            }

            // if(ip!=0){
            for (size_t AM = 1; AM <= _cell_point2[pcell][0]; AM++) {
                // ip--;
                ip = _cell_point2[pcell][AM];
                if (ip == point1 || ip == point2) {
                    continue;
                }

                if (_active[ip] == 0) {
                    cout << "Error at ReconstructDelaunay().. inactive point "
                            "in a cell "
                         << endl;
                    system("pause");
                }
                pcell_x = x[ip];  // getting x and y poisition of the point
                                  // inside the cell (pcell)
                pcell_y = y[ip];

                VerticalLine(xi, yi, pcell_x, pcell_y, qo_x, qo_y, qn_x,
                             qn_y);  // constructing the trimmign voronoi line

                Trimming(xi, yi, po_x, po_y, pn_x, pn_y, qo_x, qo_y, qn_x,
                         qn_y);  // trimming of line po,pn usinf line qo,qn

                if (po_x == xi && po_y == yi && pn_x == xi &&
                    pn_y == yi)  // if degenerated voronoi edge
                {
                    outside_edge = true;
                    break;
                }
            }
        }

        if (outside_edge == true) {
            break;
        }
    }

    double l = sqrt((po_x - pn_x) * (po_x - pn_x) +
                    (po_y - pn_y) * (po_y - pn_y));  // length of voronoi edge

    double l2 = sqrt((xi - xj) * (xi - xj) +
                     (yi - yj) * (yi - yj));  // length of delaunay edge

    if (l > eps) {
        return true;
    } else {
        return false;
    }
}

inline void DelaunayMeshGenerator()
{
    size_t i, j;
    size_t icell;
    double node_x, node_y;
    size_t ip, ipoint;
    double xp, yp;
    int    c, r;

    for (i = 0; i < _s_d + _s_d; i++) {
        _node_nodes[i][0] = 0;
    }

    // loop over each point
    for (ipoint = 0; ipoint < _num_points; ipoint++) {

        node_x = x[ipoint];
        node_y = y[ipoint];

        i = size_t((node_x - _xo) / _sx);
        j = size_t((node_y - _yo) / _sy);

        bool check_point;

        for (c = int(i - 15); c <= int(i + 15); c++) {

            for (r = int(j - 15); r <= int(j + 15); r++) {

                if (c >= _nx || r >= _ny || c < 0 || r < 0) {
                    continue;
                }

                // if(c == i && r == j) { continue; }

                icell = cell_index(c, r);  // getting the cell index

                check_point = false;

                // ip = cell_point[icell];
                if (_cell_point2[icell][0] == 0) {
                    continue;
                }
                // if(ip != 0){
                for (size_t AM = 1; AM <= _cell_point2[icell][0]; AM++) {
                    // ip--;
                    ip = _cell_point2[icell][AM];
                    if (ip == ipoint) {
                        continue;
                    }
                    xp = x[ip];  // getting the point x and y poisition
                    yp = y[ip];


                    if (CheckConnectivity(ipoint, ip)) {
                        continue;
                    }

                    if (DelaunayEdgeExist(ipoint, ip, i, j, c, r, node_x,
                                          node_y, xp, yp)) {
                        _node_nodes[ipoint][0]++;
                        _node_nodes[ipoint][_node_nodes[ipoint][0]] = ip;
                        _node_nodes[ip][0]++;
                        _node_nodes[ip][_node_nodes[ip][0]] = ipoint;

                        if (_node_nodes[ipoint][0] > 15 ||
                            _node_nodes[ip][0] > 15) {
                            cout << "Error.. way to many points connected to "
                                    "point in delaunay "
                                 << endl;
                            system("pause");
                        }

                        // Nodes_Connector(ipoint,ip-1);
                        // Nodes_Connector(ip-1,ipoint);
                    }
                }
            }
        }
    }

     for (size_t i = 0; i < _num_points; i++) {
        SortDelaunay2(i);
    }

}