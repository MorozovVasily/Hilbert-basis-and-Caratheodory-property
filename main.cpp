#include <iostream>
#include "libnormaliz/cone.h"
#include <chrono>

typedef libnormaliz::Cone<long long> Cone;
using libnormaliz::InputType;
using libnormaliz::operator<<;

using namespace std;

typedef vector<long long> vect;
typedef vector<vect> matrix;

constexpr int d = 6;

vector<vector<long long>> v5_1((1 << d) - 1, vector<long long>(d));

const vector<vector<long long>> contr_example7 = {
    {0, 1, 1, 0, 0, 0},
    {0, 1, 0, 1, 1, 0},
    {0, 1, 0, 0, 0, 1},
    {1, 0, 0, 0, 1, 0},
    {1, 0, 0, 1, 0, 0},
    {0, 1, 1, 1, 0, 0},
    {0, 1, 0, 0, 1, 1},
    {1, 0, 0, 1, 0, 1},
    {1, 0, 1, 0, 0, 1},
    {1, 0, 1, 0, 1, 0}
};

void create_5_01_vectors(vector<vector<long long>>& v5_1) {
    for (int i = 0; i < (1 << d) - 1; ++i) {
        auto i2 = i+1;
        for (int j = 0; i2 != 0; ++j) {
            v5_1[i][j] = i2 & 1;
            i2 >>= 1;
        }
    }
    sort(v5_1.begin(), v5_1.end());
}

matrix whereWeStop = {{0, 0, 0, 0, 0, 0}};

auto now() {
    return std::chrono::high_resolution_clock::now();
}

typedef std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1, 1000000000> > > tp;

std::chrono::duration<double> copy_vector, create_cone, hilbert_basis, optimize, combine, prepare, optimize_rank, optimize2;


map<string, pair<pair<std::chrono::duration<double>, tp>, int>> myTimers;
map<string, int> skipCone;

inline void startTimer(string s) {
    myTimers[s].first.second = now();
}

inline void endTimer(string s) {
    myTimers[s].first.first += now() - myTimers[s].first.second;
    myTimers[s].second += 1;
}

void printTime(pair<string, pair<pair<std::chrono::duration<double>, tp>, int>> t) {
    cout << t.first << "\t time: " << t.second.first.first.count() << "\t " << t.second.first.first.count() / t.second.second << endl;
}

int skip_optimize = 0, skip_optimize2 = 0, skip_rank = 0, skip_extrimerays = 0, skip_hb = 0;


void getUnimodulaSubcones_helper(vector<Cone> &result, const matrix& hb, int id, libnormaliz::Matrix<long long>& currentResult) {
    if (id == hb.size()) {
        return;
    }
    int n = currentResult.nr_of_rows();
    currentResult.append(hb[id]);
    if (currentResult.rank() == n+1) {
        if (currentResult.nr_of_rows() == currentResult.nr_of_columns()) {
            if (currentResult.vol() == 1) {
                result.push_back(Cone(InputType::cone, currentResult));
            }
        } else {
            getUnimodulaSubcones_helper(result, hb, id+1, currentResult);
        }
    }
    currentResult.remove_row(n);
    getUnimodulaSubcones_helper(result, hb, id+1, currentResult);
}

vector<Cone> getUnimodulaSubcones(const matrix& hb) {
    if (hb.empty()) {
        return {};
    }
    vector<Cone> result;
    libnormaliz::Matrix<long long> currentResult(0, hb[0].size());
    getUnimodulaSubcones_helper(result, hb, 0, currentResult);
    return result;
}

bool unicover(Cone& d, int n, vector<Cone>& U) {
    matrix erd = d.getExtremeRays();
    int dimension = erd[0].size();
    sort(erd.begin(), erd.end());
    for (int i = n; i < U.size(); ++i) {
        bool brk = false;
        for (const vect& hiperplane : U[i].getSupportHyperplanes()) {
            bool inside = false;
            bool outside = false;
            for (auto& er : erd) {
                int sum = 0;
                for (int x = 0; x < er.size(); ++x) {
                    sum += er[x] * hiperplane[x];
                }
                if (sum > 0) {
                    inside = true;
                }
                if (sum < 0) {
                    outside = true;
                }
            }
            if (outside && !inside) {
                brk = true;
                break;
            }
            if (outside && inside) {
                Cone c1 = Cone(InputType::cone, d.getExtremeRays(), InputType::inequalities, {hiperplane});
                assert (c1.getRank() == dimension);
                vect hiperplane_negative = hiperplane;
                for (auto& x : hiperplane_negative) {
                    x *= -1;
                }
                Cone c2 = Cone(InputType::cone, d.getExtremeRays(), InputType::inequalities, {hiperplane_negative});
                assert (c2.getRank() == dimension);
                bool bSucceed = true;
                #pragma omp parallel for reduction(&& : bSucceed)
                for (int z = 0; z <= 1; ++z) {
                    if (z == 0)
                        bSucceed &= unicover(c1, n, U);
                    else
                        bSucceed &= unicover(c2, n+1, U);
                }
                return bSucceed;
            }
        }
        if (!brk) {
            return true;
        }
    }
    cout << "finish\n";
    cout << d.getHilbertBasis();
    return false;
}

bool UHC(Cone& c) {
    startTimer("UHC1");
    auto [p1, p2] = c.getTriangulation();
    vector<Cone> U = getUnimodulaSubcones(c.getHilbertBasis());
    cout << "U";
    cout << U.size() << endl;

    // Посчитаем вспомогательные гиперплоскости зарание, 
    // чтобы когда потом они нам понадобятся данный метод вызывался константно
    // И у нас не было ошибок при распараллеливание
    // Структура не менялась из разных потоков
    for (auto& x : U) {
        x.getSupportHyperplanes();
    }
    cout << "p1.size() = " << p1.size() << '\n';

    result = true;
    #pragma omp parallel for
    for (auto& subcone : p1) {
        auto& keys = subcone.key;
        matrix coneBasis;
        for (auto x : keys) {
            coneBasis.push_back(p2[x]);
        }
        auto triangulationCone = Cone(InputType::cone, coneBasis);
        if (!unicover(triangulationCone, 0, U)) {
            cout << "subcone\n";
            cout << triangulationCone.getHilbertBasis();

            // Можно было бы написать
            // return false;
            // Но так как мы распараллелили через omp
            // Мы не можем выходить внутри цикла
            result = false;
        }
    }
    endTimer("UHC1");
    printTime(*myTimers.find("UHC1"));
    return result;
}


bool isTight(Cone& c) {
    auto hb = c.getHilbertBasis();
    sort(hb.begin(), hb.end());
    auto er = c.getExtremeRays();
    sort(er.begin(), er.end());
    vect tmp = vector<long long>(hb[0].size());
    for (int i = 0, j = 0; i < er.size(); ++i, ++j) {
        while (hb[j] != er[i]) {
            ++j;
            assert(j < hb.size());
        }
        swap(tmp, er[i]);
        auto hb2 = Cone(InputType::cone, er).getHilbertBasis();
        swap(tmp, er[i]);
        if (hb2.size() != hb.size() - 1) {
            continue;
        }
        sort(hb2.begin(), hb2.end());
        bool result = false;
        for (int k = 0; k < hb2.size(); ++k) {
            if (hb2[k] != hb[k + (k >= j)]) {
                result = true;
                break;
            }
        }
        if (result == true) {
            continue;
        }
        return false;
    }
    return true;
}

constexpr int factorial(int x) {
    if (x == 0) {
        return 1;
    } else {
        return x * factorial(x - 1);
    }
}

set<vect> smx[factorial(d)];

int qq = 0;

void combinate(matrix& m, int i, matrix& current) {
    if (i == m.size()) {
        return;
    }
    ++qq;
    if (qq % 100000 == 0) {
        cout << qq << '\n';
        cout << current << endl;
    }
    bool zz = true;
    current.push_back(m[i]);
    for (int i = 0; i < min(current.size(), whereWeStop.size()); ++i) {
        if (current[i] < whereWeStop[i]) {
            zz = false;
            break;
        } else if (current[i] > whereWeStop[i]) {
            break;
        }
    }
    auto c = Cone(InputType::cone, current);
    if (zz && c.getExtremeRays().size() == current.size()) {
        int perm[d];
        for (int i = 0; i < d; ++i) {
            perm[i] = i;
        }
        bool brk = false;
        int k = 0;
        for (; k < factorial(d); ++k) {
            next_permutation(perm, perm + d);
            vect v(d);
            for (int j = 0; j < d; ++j) {
                v[j] = m[i][perm[j]];
            }
            smx[k].insert(v);

            int j = 0;
            for (const vect& v2 : smx[k]) {
                if (v2 < current[j]) {
                    brk = true;
                    break;
                } else if (v2 > current[j]) {
                    break;
                }
                ++j;
            }
            if (brk) {
                break;
            }
        }
        if (!brk) {
            if (libnormaliz::Matrix(current).rank() == d) {
                if (isTight(c)) {
                    cout << "isTight\n";
                    cout << "hb";
                    cout << c.getHilbertBasis();
                    cout << "er";
                    cout << c.getExtremeRays();
                    bool cont = false;
                    for (auto& x : c.getHilbertBasis()) {
                        for (auto y : x) {
                            if (y > 1) {
                                cont = true;
                                break;
                            }
                        }
                        if (cont) {
                            break;
                        }
                    }
                    if (!cont) {
                        if (!UHC(c)) {
                            cout << "disprove UHC\n";
                            cout << c.getHilbertBasis();
                            exit(1);
                        }}
                }
            }
            combinate(m, i+1, current);
        }
        for (int i = 0; i < d; ++i) {
            perm[i] = i;
        }
        for (int z = 0; z <= k; ++z) {
            next_permutation(perm, perm + d);
            vect v(d);
            for (int j = 0; j < d; ++j) {
                v[j] = m[i][perm[j]];
            }
            smx[z].erase(v);
        }
    }
    current.pop_back();
    combinate(m, i+1, current);
}

#define print(elapsed) std::cout << #elapsed << "\t time: " << elapsed.count() << " s\n";
#define print2(a) cout << #a << "\t " << a << '\n';


int main() {
    matrix m;
    combinate(v5_1, 0, m);
    cout << "End\n" << endl;
}
