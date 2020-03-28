#pragma once

inline void Partition(int& i, int& j, size_t ip)
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
            tmp = _quad_nodes[ip][i + 1];
            _quad_nodes[ip][i + 1] = _quad_nodes[ip][j + 1];
            _quad_nodes[ip][j + 1] = tmp;
            tmp1 = _angles[i];
            _angles[i] = _angles[j];
            _angles[j] = tmp1;
            j--;
            i++;
        }
    }
}

inline void Sort(int left, int right, size_t ip)
{
    int i, j;
    i = left;
    j = right;
    Partition(i, j, ip);
    if (left < j) {
        Sort(left, j, ip);
    }
    if (i < right) {
        Sort(i, right, ip);
    }
}

inline void SortQuadMesh2(size_t ip)
{
    size_t num_nodes = _quad_nodes[ip][0];
    double dx, dy;
    size_t ip1;

    double xp = x[ip];
    double yp = y[ip];

    for (size_t i = 1; i <= num_nodes; i++) {

        ip1 = _quad_nodes[ip][i];
        dx = x[ip1] - xp;
        dy = y[ip1] - yp;

        _angles[i - 1] = GetAngle360(dy, dx);
    }
    Sort(0, num_nodes - 1, ip);
}

inline void QuadMesh_Constructor2()
{
    size_t V, dnum, ip, ip1;

    for (ip = 0; ip < num_points; ip++) {
        for (V = 1; V <= node_nodes[ip][0]; V++) {
            ip1 = node_nodes[ip][V];
            if (_disk_color[ip] == _disk_color[ip1]) {
                continue;
            }  // same color
            _quad_nodes[ip][0]++;
            _quad_nodes[ip][_quad_nodes[ip][0]] = ip1;
        }
    }
}