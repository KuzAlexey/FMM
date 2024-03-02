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

// ALGORITHM PARAMS
const int MIN_Y = -(1 << 29), MAX_Y = (1 << 29);
const int MIN_X = -(1 << 29), MAX_X = (1 << 29);
const long double EPS = 1e-15;
const int P = 15;
const int D = 3;
//

struct Point {
    complex<long double> z = 0, pot = 0, F = 0;
    long double q = 0;
    long double x() {
        return z.real();
    }
    long double y() {
        return z.imag();
    }
    Point(long double x, long double y, long double q_p) {
        q = q_p;
        z = complex<long double>(x, y);
    }
    int get_octet(int x, int y) {
        if (this->x() <= x && this->y() >= y) {
            return 0;// r_x = x_m, l_y = y_m
        } else if (this->x() > x && this->y() >= y) {
            return 1;// l_x = x_m, l_y = y_m
        } else if (this->x() <= x && this->y() < y) {
            return 2;// r_x = x_m, r_y = y_m
        } else if (this->x() > x && this->y() < y) {
            return 3;// l_x = x_m, r_y = y_m
        } else {
            exit(0);
        }
    }
    Point() {}
};

struct Field {
    vector<Point> points;
    vector<complex<long double> > A, B;
    Field* child_fields[4] = {nullptr, nullptr, nullptr, nullptr};
    Field* parent = nullptr;
    int x_m{}, y_m{};
    vector<Field*> neighbours, interaction_list;
    int l_y{}, r_y{}, l_x{}, r_x{};
    int level{};

    bool is_last() {
        return points.empty() || r_x - l_x == 2;
    }

    //DEBUG
    complex<long double> DEBUG_calc_p(long double x, long double y, long double q) {
        Point p(x, y, q);
        complex<long double> z_0(x_m, y_m), kek(1, 0);
        kek *= p.z - z_0;
        for (int i = 1; i < P; i++) {
            p.pot += A[i] / kek;
            kek *= p.z - z_0;
        }
        p.pot += A[0] * log(p.z - z_0);
        return p.q * p.pot;
    }

    complex<long double> DEBUG_calc_p_direct(long double x, long double y, long double q) {
        Point p(x, y, q);
        for (auto point: points) {
            p.pot += p.q * point.q * log(p.z - point.z);
        }
        return p.pot;
    }
};

vector<map<pair<int, int>, Field* > > MAP(40);
vector<vector<Field*>> fields(40);

vector<Point> input() {
    freopen("C:\\fmm_refactor\\input.txt", "r", stdin);
    int N;
    cin >> N;
    vector<Point> points;
    for (int i = 0; i < N; ++i) {
        long double x, y, m;
        cin >> x >> y >> m;
        points.push_back(Point(x, y, m));
    }
    return points;
}

void init_KDTREE(Field* f) {
    // init f
    if (f->parent == nullptr) {
        f->level = 0;
    } else {
        f->level = f->parent->level + 1;
    }
    f->x_m = (f->l_x + f->r_x) / 2;
    f->y_m = (f->l_y + f->r_y) / 2;
    MAP[f->level][{f->l_x, f->l_y}] = f;
    fields[f->level].push_back(f);
    f->A.resize(P);
    f->B.resize(P);

    if (f->is_last()) {
        return;
    }

    // init childs
    for (auto point: f->points) {
        int octet = point.get_octet(f->x_m, f->y_m);
        if (f->child_fields[octet] == nullptr) {
            f->child_fields[octet] = new Field;
            f->child_fields[octet]->parent = f;
            f->child_fields[octet]->l_x = f->l_x;
            f->child_fields[octet]->r_x = f->r_x;
            f->child_fields[octet]->l_y = f->l_y;
            f->child_fields[octet]->r_y = f->r_y;
            if (octet == 0) {
                f->child_fields[octet]->r_x = f->x_m; f->child_fields[octet]->l_y = f->y_m;
            } else if (octet == 1) {
                f->child_fields[octet]->l_x = f->x_m; f->child_fields[octet]->l_y = f->y_m;
            } else if (octet == 2) {
                f->child_fields[octet]->r_x = f->x_m; f->child_fields[octet]->r_y = f->y_m;
            } else {
                f->child_fields[octet]->l_x = f->x_m; f->child_fields[octet]->r_y = f->y_m;
            }
        }
        f->child_fields[octet]->points.push_back(point);
    }
    for (int i = 0; i < 4; i++) {
        if (f->child_fields[i] != nullptr) {
            init_KDTREE(f->child_fields[i]);
        }
    }
}

void init_near() {
    for (int level = 0; level < fields.size(); level++) {
        for (auto f: fields[level]) {
            Field* parent = f->parent;
            if (parent == nullptr) {
                continue;
            }
            int D_parent = parent->r_x - parent->l_x;
            int D_f = f->r_x - f->l_x;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    int nl_x = parent->l_x + i * D_parent;
                    int nl_y = parent->l_y + j * D_parent;
                    Field* neigh_parent = MAP[parent->level][{nl_x, nl_y}];
                    if (neigh_parent == nullptr) {
                        continue;
                    }
                    for (auto neigh: neigh_parent->child_fields) {
                        if (neigh == nullptr || neigh == f) {
                            continue;
                        }
                        int dx = abs(neigh->l_x - f->l_x);
                        int dy = abs(neigh->l_y - f->l_y);
                        if (dx >= 2 * D_f || dy >= 2 * D_f) {
                            f->interaction_list.push_back(neigh);
                        } else {
                            f->neighbours.push_back(neigh);
                        }
                    }
                }
            }
        }
    }
}

complex<long double> my_pow(complex<long double> z, int k) {
    complex<long double> zzz = z;
    if (k == 0) {
        return 1;
    }
    for (int i = 1; i < k; i++) {
        z *= zzz;
    }
    return z;
}

void init_a() {
    for (int level = 0; level < fields.size(); level++) {
        for (auto f: fields[level]) {
            f->A.resize(P);

            // init a
            f->A[0] = 0;
            for (auto point: f->points) {
                f->A[0] += point.q;
            }
            complex<long double> z_0(f->x_m, f->y_m);
            for (int k = 1; k < P; k++) {
                f->A[k] = 0;
                for (auto point: f->points) {
                    f->A[k] += -point.q * my_pow(point.z - z_0, k) / ((long double) k);
                }
            }
        }
    }
}

long double CNK[2 * P][2 * P];
bool USED[2 * P][2 * P];
long double C(int n, int k) {
    if (k == 0 || k == n)
        return 1;
    if (USED[n][k]) {
        return CNK[n][k];
    }
    USED[n][k] = true;
    CNK[n][k] = C(n - 1, k - 1) * n / k;
    return CNK[n][k];
}

void init_b() {
    for (int level = 0; level < fields.size(); level++) {
        for (auto f: fields[level]) {

            // interaction_list A -> B
            complex<long double> delta = 0, z_0(f->x_m, f->y_m);
            for (auto near: f->interaction_list) {
                complex<long double> z_n(near->x_m, near->y_m);
                complex<long double> z_0_z_n_l(1, 0);
                for (int l = 0; l < P; l++) {
                    delta = 0;
                    if (l == 0) {
                        complex<long double> z_0_z_n = z_0 - z_n;
                        for (int k = 1; k < P; k++) {
                            delta += near->A[k] / z_0_z_n;
                            z_0_z_n *= z_0 - z_n;
                        }
                        delta += near->A[0] * log(z_0 - z_n);
                    } else {
                        complex<long double> z_0_z_n = z_0 - z_n;
                        for (int k = 1; k < P; k++) {
                            delta += near->A[k] / z_0_z_n * (C(l + k - 1, k - 1));
                            z_0_z_n *= z_0 - z_n;
                        }
                        delta /= z_0_z_n_l * (long double) pow(-1, l % 2);
                        delta += near->A[0] / ((long double)(l) * z_0_z_n_l) * (long double) pow(-1, l % 2 + 1);
                    }
                    z_0_z_n_l *= z_0 - z_n;
                    f->B[l] += delta;
                }
            }

            // B -> B to childs
            for (auto child: f->child_fields) {
                if (child == nullptr) {
                    continue;
                }
                complex<long double> z_n(child->x_m, child->y_m), z_0_z_n_l = 1, z_0_z_n_k = 1;
                for (int l = 0; l < P; l++) {
                    z_0_z_n_k = z_0_z_n_l;
                    for (int k = l; k < P; k++) {
                        if (k == l) {
                            child->B[l] += C(k, l) * f->B[k];
                        } else {
                            child->B[l] += C(k, l) * f->B[k] * z_0_z_n_k / z_0_z_n_l;
                        }
                        z_0_z_n_k *= z_n - z_0;
                    }
                    z_0_z_n_l *= z_n - z_0;
                }
            }
        }
    }
}

void cacl_pf(Field *f, Point *p) {
    if (f->is_last()) {
        p->pot = 0;
        p->F = 0;
        complex<long double> z_0(f->x_m, f->y_m), fff(0, 0);
        for (int l = 0; l < P; l++) {
            p->pot += p->q * f->B[l] * pow(p->z - z_0, l);
            if (l != 0) {
                complex<long double> ff = l * p->q * f->B[l] * pow(p->z - z_0, l - 1);
                ff.real(-ff.real());
                p->F -= ff;
                fff += abs(ff);
            }
        }
        for (auto neigh: f->neighbours) {
            if (neigh == nullptr) {
                continue;
            }
            for (auto pn: neigh->points) {
                if (abs(p->z - pn.z) < EPS) {
                    continue;
                }
                p->pot += p->q * pn.q * log(p->z - pn.z);

                auto r2 = pow(abs(p->z - pn.z),2);
                complex<long double> v(p->z.real() - pn.z.real(), p->z.imag() - pn.z.imag());
                p->F += p->q * pn.q * v / r2;
            }
        }
        for (auto inp: f->points) {
            if (abs(p->z - inp.z) < EPS) {
                continue;
            }
            p->pot += p->q * inp.q * log(p->z - inp.z);

            auto r2 = pow(abs(p->z - inp.z),2);
            auto r3 = pow(abs(p->z - inp.z),3);
            complex<long double> v(p->z.real() - inp.z.real(), p->z.imag() - inp.z.imag());
            p->F += p->q * inp.q * v / r2;
        }
    } else {
        cacl_pf(f->child_fields[p->get_octet(f->x_m, f->y_m)], p);
    }
}

vector<Point> direct_calc_pf(vector<Point> points) {
    vector<Point> updated_points;
    for (auto cur_p: points) {
        for (auto p: points) {
            if (abs(cur_p.z - p.z) < EPS) {
                continue;
            }
            cur_p.pot += cur_p.q * p.q * log(cur_p.z - p.z);

            auto r2 = pow(abs(cur_p.z - p.z),2);
            complex<long double> v(cur_p.z.real() - p.z.real(), cur_p.z.imag() - p.z.imag());
            cur_p.F += cur_p.q * p.q * v / r2;
        }
        updated_points.emplace_back(cur_p);
    }
    return updated_points;
}


void cerr_errors_pf(vector<Point> correct, vector<Point> my_answer) {
    long double F_ERR = 0, POT_ERR = 0;
    for (int i = 0; i < correct.size(); i++) {
        F_ERR += abs(correct[i].F - my_answer[i].F);
        POT_ERR += abs(correct[i].pot.real() - my_answer[i].pot.real());
    }
    F_ERR /= correct.size();
    POT_ERR /= correct.size();

    cerr << "F absolute error = " << F_ERR << "\n";
    cerr << "POT_X absolute error = " << POT_ERR << "\n";
}

int main() {
    Field *f = new Field;
    vector<Point> points = input();
    int N = points.size();
    f->points = points;
    f->l_y = MIN_Y; f->r_y = MAX_Y; f->l_x = MIN_X; f->r_x = MAX_X;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(), beginALL =  std::chrono::steady_clock::now();
    init_KDTREE(f);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cerr << "Init KD Tree Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_near();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init Interaction List Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_a();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init A-multipoles Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    init_b();
    end = std::chrono::steady_clock::now();
    std::cerr << "Init B-multipoles Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;

    begin = std::chrono::steady_clock::now();
    for (int i = 0; i < N; i++) {
        cacl_pf(f, &points[i]);
    }
    end = std::chrono::steady_clock::now();
    std::cerr << "Calcualtions Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;
    std::cerr << "ALL FMM Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - beginALL).count()) / 1e6 << "[s]" << std::endl;



    begin = std::chrono::steady_clock::now();
    vector<Point> correct_points = direct_calc_pf(f->points);
    end = std::chrono::steady_clock::now();
    std::cerr << "Direct calcualtion Time = " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1e6 << "[s]" << std::endl;
    cerr_errors_pf(correct_points, points);

    return 0;
}
