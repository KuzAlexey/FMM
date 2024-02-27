#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdio>
#include <thread>
#include <complex>
#include <map>


using namespace std;

namespace AlgorithmParams {
    int p = 30;
    int D = 1;
//    double MIN_X = -(1 << 20);
//    double MAX_X = (1 << 20);
//    double MIN_Y = -(1 << 20);
//    double MAX_Y = (1 << 20);
    double MIN_X = -(1 << 7);
    double MAX_X = (1 << 7);
    double MIN_Y = -(1 << 7);
    double MAX_Y = (1 << 7);
    double eps = 1e-10;
};

struct Point {
    complex<double> z, pot;
    double q;

    double x() {
        return z.real();
    }
    double y() {
        return z.imag();
    }

    Point(double x, double y, double q_p) {
        q = q_p;
        z = complex<double>(x, y);
    }
    Point() {

    }
};

struct Field {
    vector<Point> field_points;
    vector<complex<double> > A, B;
    vector<Field*> child_fields;
    vector<Field*> neighbours;
    vector<Field*> interaction_list;

    Point p_c;
    Point z_0;
    double Q = 0;
    double l_x, r_x, l_y, r_y;

    Field() {

    }

    void init_z_0_Q_A() {
        Q = 0;
        double sum_qx = 0, sum_qy = 0;
        for (auto p: field_points) {
            sum_qx += p.q * p.x();
            sum_qy += p.q * p.y();
            Q += p.q;
        }
        z_0 = Point(sum_qx / Q, sum_qy / Q, 0);
        A[0] = Q;
        for (int k = 1; k < AlgorithmParams::p; k++) {
            A[k] = 0;
            for (auto & field_point : field_points) {
                complex<double> zi(field_point.x(), field_point.y());
                complex<double> z0(z_0.x(), z_0.y());
                complex<double> zi_z0 = zi - z0;
                complex<double> zi_z0_k = pow(zi_z0, k);
                A[k] += -field_point.q / k * zi_z0_k;
            }
        }
    }

    bool isLast() {
        return int(r_x - l_x) == 1;
    }
};

map<pair<pair<int, int>, pair<int, int> >, Field*> mp;

int C(int n, int k) {
    if (k == 0 || k == n)
        return 1;
    return C(n - 1, k - 1) * n / k;
}

void init_interaction_list(Field *f, Field *parent) {
    if (f == nullptr) {
        return;
    }
    if (parent != nullptr) {
        int D = parent->r_x - parent->l_x;
        int Dp = f->r_x - f->l_x;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                double nl_x = parent->l_x + i * D;
                double nr_x = parent->r_x + i * D;
                double nl_y = parent->l_y + j * D;
                double nr_y = parent->r_y + j * D;
                auto par_neigh = mp[{{nl_x, nr_x},
                                     {nl_y, nr_y}}];
                if (par_neigh == nullptr) {
                    continue;
                }
                for (auto f_neigh: par_neigh->child_fields) {
                    if (f_neigh == nullptr) {
                        continue;
                    }
                    int dx = abs(f_neigh->l_x - f->l_x);
                    int dy = abs(f_neigh->l_y - f->l_y);
                    if (dx >= 2 * Dp || dy >= 2 * Dp) {
                        f->interaction_list.push_back(f_neigh);
                    } else if (dx != 0 || dy != 0) {
                        f->neighbours.push_back(f_neigh);
                    }
                }
            }
        }
    }
    if (f->isLast()) {
        return;
    }
    init_interaction_list(f->child_fields[0], f);
    init_interaction_list(f->child_fields[1], f);
    init_interaction_list(f->child_fields[2], f);
    init_interaction_list(f->child_fields[3], f);
}

void init_B(Field* f) {
    if (f == nullptr) {
        return;
    }
    for (auto neigh: f->interaction_list) {
        if (neigh == nullptr || neigh->field_points.empty()) {
            continue;
        }
        if (f->l_x == 0 && f->r_x == 2 && f->l_y == 0) {
            cerr << "";
        }
        auto z_00 = f->z_0.z - neigh->z_0.z;
        complex<double> delta = 0;
        for (int k = 1; k < AlgorithmParams::p; k++) {
            delta += neigh->A[k] / pow(z_00, k);
            // f->B[0] += neigh->A[k] / pow(z_00, k);
        }
        delta += neigh->A[0] * log(z_00);
        f->B[0] += delta;
        //f->B[0] += neigh->A[0] * log(z_00);
        for (int l = 1; l < AlgorithmParams::p; l++) {
            delta = 0;
            for (int k = 1; k < AlgorithmParams::p; k++) {
                delta += neigh->A[k] / pow(z_00, k) * double(C(l + k - 1, k - 1));
                // f->B[l] += neigh->A[k] / pow(z_00, k) * double(C(l + k - 1, k - 1));
            }
            delta /= pow(z_00, l) * pow(-1, l);
            delta += neigh->A[0] / (double(l) * pow(z_00, l)) * pow(-1, l + 1);
            f->B[l] += delta;
            // f->B[l] /= pow(z_00, l); /// ???
            // f->B[l] += neigh->A[0] / (double(l) * pow(z_00, l)) * pow(-1, l + 1);
        }
    }

    if (f->isLast() || f->field_points.empty()) {
        return;
    }

    for (int l = 0; l < AlgorithmParams::p; l++) {
        for (int k = l; k < AlgorithmParams::p; k++) {
            f->child_fields[0]->B[l] += double(C(k, l)) * f->B[k] * pow(-(f->z_0.z - f->child_fields[0]->z_0.z), k - l);
            f->child_fields[1]->B[l] += double(C(k, l)) * f->B[k] * pow(-(f->z_0.z - f->child_fields[1]->z_0.z), k - l);
            f->child_fields[2]->B[l] += double(C(k, l)) * f->B[k] * pow(-(f->z_0.z - f->child_fields[2]->z_0.z), k - l);
            f->child_fields[3]->B[l] += double(C(k, l)) * f->B[k] * pow(-(f->z_0.z - f->child_fields[3]->z_0.z), k - l);
        }
    }
    for (auto child: f->child_fields) {
        init_B(child);
    }
}

void init_field(Field* f, double l_x, double r_x, double l_y, double r_y) {
    mp[{{l_x, r_x}, {l_y, r_y}}] = f;
    f->l_x = l_x;
    f->l_y = l_y;
    f->r_x = r_x;
    f->r_y = r_y;
    f->child_fields = {new Field, new Field, new Field, new Field};

    f->p_c = Point((r_x + l_x) / 2, (r_y + l_y) / 2, 0);
    f->A.resize(AlgorithmParams::p);
    f->B.resize(AlgorithmParams::p);

    f->init_z_0_Q_A();

    if (f->isLast()) {
        return;
    }

    for (auto p: f->field_points) {
        if (p.x() < f->p_c.x() && p.y() < f->p_c.y()) {
            f->child_fields[0]->field_points.push_back(p);
        } else if (p.x() >= f->p_c.x() && p.y() < f->p_c.y()) {
            f->child_fields[1]->field_points.push_back(p);
        } else if (p.x() < f->p_c.x() && p.y() >= f->p_c.y()) {
            f->child_fields[2]->field_points.push_back(p);
        } else {
            f->child_fields[3]->field_points.push_back(p);
        }
    }
    init_field(f->child_fields[0], l_x, f->p_c.x(), l_y, f->p_c.y());
    init_field(f->child_fields[1], f->p_c.x(), r_x, l_y, f->p_c.y());
    init_field(f->child_fields[2], l_x, f->p_c.x(), f->p_c.y(), r_y);
    init_field(f->child_fields[3], f->p_c.x(), r_x, f->p_c.y(), r_y);
}

void calc_p(Field *f, Point *p) {
    if (f->field_points.size() == 1 && f->isLast()) {
        p->pot += p->q * f->B[0];
        for (auto neigh: f->neighbours) {
            if (neigh == nullptr) {
                continue;
            }
            for (auto pn: neigh->field_points) {
                p->pot += p->q * pn.q * log(p->z - pn.z);
            }
        }
        return;
    }
    if (p->x() < f->p_c.x() && p->y() < f->p_c.y()) {
        calc_p(f->child_fields[0], p);
    } else if (p->x() >= f->p_c.x() && p->y() < f->p_c.y()) {
        calc_p(f->child_fields[1], p);
    } else if (p->x() < f->p_c.x() && p->y() >= f->p_c.y()) {
        calc_p(f->child_fields[2], p);
    } else {
        calc_p(f->child_fields[3], p);
    }
}

vector<Point> direct_calc_power(vector<Point> points) {
    vector<Point> updated_points;
    for (auto cur_p: points) {
        for (auto p: points) {
            if (abs(cur_p.z - p.z) < AlgorithmParams::eps) {
                continue;
            }
            cur_p.pot += cur_p.q * p.q * log(cur_p.z - p.z);
        }
        updated_points.emplace_back(cur_p);
    }
    return updated_points;
}

void cerr_errors(vector<Point> correct, vector<Point> my_answer) {
    double F_ERR = 0;
    for (int i = 0; i < correct.size(); i++) {
        if (i < 50) {
            cerr << correct[i].pot << " " << my_answer[i].pot << "\n";
        }
        if (abs(correct[i].pot.real() - my_answer[i].pot.real()) > 0.05) {
            cerr << i << "\n";
            cerr << correct[i].pot.real() << " " << my_answer[i].pot.real() << " " << correct[i].pot.real() - my_answer[i].pot.real() << "\n";
        }
        F_ERR += abs(correct[i].pot.real() - my_answer[i].pot.real()) / abs(correct[i].pot.real());
    }
    F_ERR /= correct.size();

    cerr << "F error = " << F_ERR * 100 << "%\n";
}

int main() {
    freopen("C:\\FMM\\input.txt", "r", stdin);

    int N;
    cin >> N;

    Field *f = new Field;
    vector<Point> points;
    for (int i = 0; i < N; ++i) {
        double m, v_x, v_y, v_z, x, y, z;
        cin >> x >> y >> m;
        f->field_points.emplace_back(x, y, m);
    }
    points = f->field_points;


    freopen("C:\\FMM\\output.txt", "w", stdout);
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    init_field(f, AlgorithmParams::MIN_X, AlgorithmParams::MAX_X, AlgorithmParams::MIN_Y, AlgorithmParams::MAX_Y);
    init_interaction_list(f, nullptr);
    init_B(f);

    for (int i = 0; i < N; i++) {
        calc_p(f, &points[i]);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Execution Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    vector<Point> correct_points = direct_calc_power(f->field_points);

    cerr_errors(correct_points, points);
    return 0;
}

/*
6
0.5 0.5 1
1.5 1.5 1
-1.5 1.5 1
-3.5 -2.5 1
3.4 0.4 1
-3.5 1.5 1
 */