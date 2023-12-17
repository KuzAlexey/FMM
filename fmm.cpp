#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <thread>

using namespace std;

int P = 2;
const double EPS = 10;

struct point {
    double x, y, z;
};

struct field {
    double sum_x{0}, sum_y{0}, sum_z{0};
    double sum_xy{0}, sum_xz{0}, sum_yz{0};
    double sum_x2{0}, sum_y2{0}, sum_z2{0};
    double x_c{0}, y_c{0}, z_c{0};
    double x_d{0}, y_d{0}, z_d{0};
    vector<point > vec;
    double Fx, Fy;
    int Nc{0};
    field* nxt[8];

    double mm1() {
        return sum_x - Nc * x_c;
    }

    double mm2() {
        return sum_y - Nc * y_c;
    }

    double mm3() {
        return sum_z - Nc * z_c;
    }

    double mm4() {
        return Nc * x_c * x_c + sum_x2 - 2 * x_c * sum_x;
    }

    double mm5() {
        return Nc * y_c * y_c + sum_y2 - 2 * y_c * sum_y;
    }

    double mm6() {
        return Nc * z_c * z_c + sum_z2 - 2 * z_c * sum_z;
    }

    double mm7() {
        return Nc * (x_c * y_c) - x_c * sum_y - y_c * sum_x + sum_xy;
    }

    double mm8() {
        return Nc * (x_c * z_c) - x_c * sum_z - z_c * sum_x + sum_xz;
    }

    double mm9() {
        return Nc * (y_c * z_c) - y_c * sum_z - z_c * sum_y + sum_yz;
    }
};

vector<point> get(vector<point> vec, double x, double y, double z, bool upx, bool upy, bool upz) {
    vector<point> res;
    for (int i = 0; i < vec.size(); ++i) {
        if (((vec[i].x < x) ^ upx) && ((vec[i].y < y) ^ upy) && ((vec[i].z < z) ^ upz)) {
            res.push_back(vec[i]);
        }
    }
    return res;
}

void init_mulitpoles(field* v, vector<point> &points) {
    if (points.size() == 0) {
        return;
    }
    double min_x = 1e9, max_x = -1e9;
    double min_y = 1e9, max_y = -1e9;
    double min_z = 1e9, max_z = -1e9;
    for (int i = 0; i < points.size(); ++i) {
        min_x = min(min_x, points[i].x);
        max_x = max(max_x, points[i].x);
        min_y = min(min_y, points[i].y);
        max_y = max(max_y, points[i].y);
        min_z = min(min_y, points[i].z);
        max_z = max(max_y, points[i].z);
    }
    double x_d = (min_x + max_x) / 2.0, y_d = (min_y + max_y) / 2.0, z_d = (min_z + max_z) / 2.0;
    v->x_d = x_d;
    v->y_d = y_d;
    v->z_d = z_d;
    if (points.size() < P) {
        v->vec = points;
        for (auto xyz: points) {
            double x = xyz.x;
            double y = xyz.y;
            double z = xyz.z;
            v->sum_x += x;
            v->sum_y += y;
            v->sum_z += z;
            v->sum_x2 += x * x;
            v->sum_y2 += y * y;
            v->sum_z2 += z * z;
            v->sum_xy += x * y;
            v->sum_xz += x * z;
            v->sum_yz += y * z;
        }
        v->Nc = points.size();
        if (v->Nc != 0) {
            v->x_c = v->sum_x / v->Nc;
            v->y_c = v->sum_y / v->Nc;
            v->z_c = v->sum_z / v->Nc;
        }
    } else {
        vector<point> p1, p2, p3, p4, p5, p6, p7, p8;
        p1 = get(points, x_d, y_d, z_d, 0, 0, 0);
        p2 = get(points, x_d, y_d, z_d, 0, 0, 1);
        p3 = get(points, x_d, y_d, z_d, 0, 1, 0);
        p4 = get(points, x_d, y_d, z_d, 0, 1, 1);
        p5 = get(points, x_d, y_d, z_d, 1, 0, 0);
        p6 = get(points, x_d, y_d, z_d, 1, 0, 1);
        p7 = get(points, x_d, y_d, z_d, 1, 1, 0);
        p8 = get(points, x_d, y_d, z_d, 1, 1, 1);
        for (int i = 0; i < 8; ++i) {
            v->nxt[i] = new field;
        }
        init_mulitpoles(v->nxt[0], p1);
        init_mulitpoles(v->nxt[1], p2);
        init_mulitpoles(v->nxt[2], p3);
        init_mulitpoles(v->nxt[3], p4);
        init_mulitpoles(v->nxt[4], p5);
        init_mulitpoles(v->nxt[5], p6);
        init_mulitpoles(v->nxt[6], p7);
        init_mulitpoles(v->nxt[7], p8);
        for (int i = 0; i < 8; ++i) { // Multipole to Multipole
            v->sum_x += v->nxt[i]->sum_x;
            v->sum_y += v->nxt[i]->sum_y;
            v->sum_z += v->nxt[i]->sum_z;

            v->sum_x2 += v->nxt[i]->sum_x2;
            v->sum_y2 += v->nxt[i]->sum_y2;
            v->sum_z2 += v->nxt[i]->sum_z2;

            v->sum_xy += v->nxt[i]->sum_xy;
            v->sum_xz += v->nxt[i]->sum_xz;
            v->sum_yz += v->nxt[i]->sum_yz;
        }
        v->Nc = points.size();
        v->x_c = v->sum_x / v->Nc;
        v->y_c = v->sum_y / v->Nc;
        v->z_c = v->sum_z / v->Nc;
    }
}

int get_next(field* v, double x, double y, double z) {
    bool bx = x < v->x_d;
    bool by = y < v->y_d;
    bool bz = z < v->z_d;
    return ((!bx) << 2) + ((!by) << 1) + ((!bz) << 0);
}

double dist(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
}

double MultipoleToPoint(field *v, double x, double y, double z, double &FX, double &FY, double &FZ) { // Multipole to Point
    if (v == nullptr) {
        return 0;
    }
    if (v->Nc == 0) {
        return 0;
    }
    double x_c = v->x_c, y_c = v->y_c, z_c = v->z_c;
    double R = dist(v->x_c, v->y_c, v->z_c, x, y, z);
    double XR = x - x_c;
    double YR = y - y_c;
    double ZR = z - z_c;
    int Nc = v->Nc;

    if (R < EPS) {
        return 0;
    }

    //cerr << v->sum_x2 << " " << v->sum_y2 << '\n';

    double P = Nc / (R * R) + v->mm1() * (-(x_c - x) / (R * R * R * R) * 2) + v->mm2() * (-(y_c - y) / (R * R * R * R) * 2) + \
            v->mm3() * (-(z_c - z) / (R * R * R * R) * 2);
//    ((4 * (x_c - x)*(x_c - x) - R * R) * v->mm4() / (R * R * R * R * R * R)) + \
//    ((4 * (y_c - y)*(y_c - y) - R * R) * v->mm5() / (R * R * R * R * R * R)) + \
//    ((4 * (z_c - z)*(z_c - z) - R * R) * v->mm6() / (R * R * R * R * R * R)) + \
//    2 * (4 * (x_c - x)*(y_c - y) / (R * R * R * R * R * R) * v->mm7()) + \
//    2 * (4 * (x_c - x)*(z_c - z) / (R * R * R * R * R * R) * v->mm8()) + \
//    2 * (4 * (y_c - y)*(z_c - z) / (R * R * R * R * R * R) * v->mm9());
    // P = abs(P);
    FX += XR / R * P;
    FY += YR / R * P;
    FZ += ZR / R * P;

    return P;
}

void calculate_potencial(field *v, double x, double y, double z, double &res, double &FX, double &FY, double &FZ) {
    if (v == nullptr) {
        return;
    }
    int ind = get_next(v, x, y, z);
    for (int i = 0; i < 8; ++i) {
        if (i != ind) {
            res += MultipoleToPoint(v->nxt[i], x, y, z, FX, FY, FZ);
        }
    }
    if (v->nxt[ind] != nullptr && v->nxt[ind]->Nc > 0) {
        calculate_potencial(v->nxt[ind], x, y, z, res, FX, FY, FZ);
    } else {
        // should be direct summing
        for (auto point : v->vec) {
            if (abs(x - point.x) > EPS || abs(y - point.y) > EPS || abs(z - point.z) > EPS) {
                double R = dist(point.x, point.y, point.z, x, y, z);
                double XR = x - point.x;
                double YR = y - point.y;
                double ZR = y - point.z;
                double P = 1 / (R * R);
                res += 1 / (R * R);
                FX += XR / R * P;
                FY += YR / R * P;
                FZ += ZR / R * P;
            }
        }
    }
}

vector<double> direct_sum(vector<point> vec) {
    vector<double> vec_p;
    cerr << "SIZE: " << vec.size() << '\n';
    for (auto xyz: vec) {
        double x = xyz.x;
        double y = xyz.y;
        double z = xyz.z;
        double res = 0;
        for (auto point : vec) {
            double R = dist(point.x, point.y, point.z, x, y, z);
            if (R >= EPS) {
                res += 1 / (R * R);
            }
        }
        vec_p.push_back(res);
        // cout << "x = " << x << " y = " << y << " P = " << res << '\n';
    }
    return vec_p;
}

void calcMAE(int N, vector<point > vec, vector<double> p_mm) {
    cerr << "START EXECUTION CALC MAE\n";
    vector<double> p;
    p = direct_sum(vec);
    double MAE = 0;
    for (int i = 0; i < N; ++i) {
        if (i < 10) {
            cerr << p[i] << " " << p_mm[i] << '\n';
        }
        MAE += abs(p[i] - p_mm[i]) / p[i];
    }
    MAE /= N;
    cerr << "MAE = " << MAE << '\n';
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main() {
    freopen("input.txt", "r", stdin);
    int N;
    cin >> N;
    std::cerr << "CALCULATE FOR N = " << N << '\n';
    vector<point > vec;
    for (int i = 0; i < N; ++i) {
        double x, y, z;
        cin >> x >> y >> z;
        // cerr << "COORD = {" << x << " " << y << " " << z << "}\n";
//         x = fRand(0, 1000);
//         y = fRand(0, 1000);
//         z = fRand(0, 1000);
        vec.push_back({x, y, z});
    }

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    field* root = new field;
    init_mulitpoles(root, vec); // nlog(n)
    vector<double> p_mm;
    double DT = 3;
    freopen("input.txt", "w", stdout);
    cout << N << '\n';
    for (int i = 0; i < vec.size(); ++i) { // n
        double res = 0, FX = 0, FY = 0, FZ = 0;
        calculate_potencial(root, vec[i].x, vec[i].y, vec[i].z, res, FX, FY, FZ); // log(n)
        p_mm.push_back(res);

        cout << vec[i].x - DT * FX << ' ' << vec[i].y - DT * FY << ' ' << vec[i].z - DT * FZ << '\n';
        // cout << "x = " << vec[i].first << " y = " << vec[i].second << " P = " << res << "FX, FY = " << FX << ' ' << FY << '\n';
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // std::cerr << "Execution Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;
    // calcMAE(N, vec, p_mm);
}
