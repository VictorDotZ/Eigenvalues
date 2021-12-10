#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#define M(i, j) matrix[(i)*n + (j)]
#define EPSILON 1e-100

double f1(size_t n, size_t i, size_t j) { return n - std::max(i, j); }

double f2(size_t n, size_t i, size_t j)
{
    if (i == j)
        return 2;
    if (std::abs(0.0 + i - j) == 1)
        return -1;
    return 0;
}

double f3(size_t n, size_t i, size_t j)
{
    if (i == n - 1)
        return 1.0 + j;
    if (j == n - 1)
        return 1.0 + i;
    return static_cast<double>(i == j);
}

double f4(size_t n, size_t i, size_t j) { return 1.0 / (i + j + 1.0); }

std::vector<double (*)(size_t n, size_t i, size_t j)> a = { nullptr, &f1, &f2, &f3, &f4 };

void fillFromFunction(size_t n, double* matrix, double (*f)(size_t n, size_t i, size_t j))
{
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            M(i, j) = f(n, i, j);
}

bool fillFromFile(size_t n, double* matrix, char* path)
{
    std::ifstream file(path);
    if (!file.is_open())
        return false;

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            file >> M(i, j);

    if (!file.eof()) {
        file.close();
        return false;
    }

    file.close();
    return true;
}

void print(size_t n, double* vector)
{
    for (size_t i = 0; i < n; ++i)
        std::cout << std::setprecision(3) << std::scientific << vector[i] << ' ';
    std::cout << std::endl;
}

void print(size_t n, size_t m, double* matrix)
{
    for (size_t i = 0; i < m; ++i)
        print(m, matrix + (i * n));
    std::cout << std::endl;
}

void rotate(size_t n, double* matrix)
{
    double tmp1;
    double tmp2;
    double r;
    double cosPhi;
    double sinPhi;

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            tmp1 = M(i, i - 1);
            tmp2 = M(j, i - 1);

            r = std::sqrt(tmp1 * tmp1 + tmp2 * tmp2);

            if (r < EPSILON)
                continue;

            cosPhi = tmp1 / r;
            sinPhi = -tmp2 / r;

            M(i, i - 1) = r;
            M(j, i - 1) = 0.0;

            for (size_t k = i; k < n; ++k) {
                tmp1 = M(i, k);
                tmp2 = M(j, k);

                M(i, k) = tmp1 * cosPhi - tmp2 * sinPhi;
                M(j, k) = tmp1 * sinPhi + tmp2 * cosPhi;
            }

            for (size_t k = 0; k < n; ++k) {
                tmp1 = M(k, i);
                tmp2 = M(k, j);

                M(k, i) = tmp1 * cosPhi - tmp2 * sinPhi;
                M(k, j) = tmp1 * sinPhi + tmp2 * cosPhi;
            }
        }
    }
}

void QRDecompose(size_t n, double* matrix, size_t k, double* x1, double* x2)
{
    double tmp1;
    double tmp2;
    for (size_t i = 0; i < k - 1; ++i) {
        tmp1 = M(i + 1, i) * M(i + 1, i);

        if (tmp1 < EPSILON) {
            tmp2 = std::fabs(M(i, i));
            x1[i] = (M(i, i) > 0.0 ? 1.0 : -1.0);
            x2[i] = 0.0;
        } else {
            tmp2 = std::sqrt(M(i, i) * M(i, i) + tmp1);

            M(i, i) -= tmp2;

            tmp1 = std::sqrt(M(i, i) * M(i, i) + tmp1);
            x1[i] = M(i, i) / tmp1;
            x2[i] = M(i + 1, i) / tmp1;
        }

        for (size_t j = i + 1; j < n; ++j) {
            tmp1 = x1[i] * M(i, j);
            tmp1 += x2[i] * M(i + 1, j);

            tmp1 *= 2.0;

            M(i, j) -= tmp1 * x1[i];
            M(i + 1, j) -= tmp1 * x2[i];
        }

        M(i, i) = tmp2;
        M(i + 1, i) = 0.0;
    }
}

void shift(size_t n, double* matrix, size_t k, double shiftValue)
{
    for (size_t i = 0; i < k; ++i) {
        M(i, i) -= shiftValue;
    }
}

double norm(size_t n, double* matrix)
{
    double tmp;
    double res = 0.0;

    for (size_t i = 0; i < n; ++i) {
        tmp = 0.0;
        for (size_t j = 0; j < n; ++j)
            tmp += std::fabs(M(i, j));

        if (res < tmp)
            res = tmp;
    }

    return res;
}

double inv1(size_t n, double* matrix)
{
    double res = 0.0;
    for (size_t i = 0; i < n; ++i) {
        res -= M(i, i);
    }
    return res;
}

double inv2(size_t n, double* matrix)
{
    double res = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res -= M(i, j) * M(j, i);
        }
    }
    return res;
}

void mul(size_t n, double* matrix, size_t k, double* x1, double* x2)
{
    double tmp;

    for (size_t i = 0; i < k - 1; ++i) {
        for (size_t j = 0; j < i + 2; ++j) {
            tmp = M(j, i) * x1[i];
            tmp += M(j, i + 1) * x2[i];

            tmp *= 2.0;

            M(j, i) -= tmp * x1[i];
            M(j, i + 1) -= tmp * x2[i];
        }
    }
}

void findEigenvalues(size_t n, double* matrix, double* eigenvalues, double epsilon, double* x1, double* x2)
{
    double t = norm(n, matrix) * epsilon;
    double s = 0;
    rotate(n, matrix);

    for (size_t k = n; k > 2; --k) {
        size_t currentCycleStep = 0;
        while (std::fabs(M(k - 1, k - 2)) > t) {
            s = M(k - 1, k - 1);

            if (currentCycleStep++ > 10) {
                s = M(k - 1, k - 1) + 0.1;
                currentCycleStep = 0;
            }

            shift(n, matrix, k, s);

            QRDecompose(n, matrix, k, x1, x2);
            mul(n, matrix, k, x1, x2);

            shift(n, matrix, k, -s);
        }
    }

    if (n > 1) {
        t = M(0, 0) + M(1, 1);
        s = M(0, 0) * M(1, 1) - M(0, 1) * M(1, 0);
        s = std::sqrt(t * t - 4.0 * s);

        M(0, 0) = 0.5 * (t + s);
        M(1, 1) = 0.5 * (t - s);
    }

    for (size_t i = 0; i < n; ++i) {
        eigenvalues[i] = M(i, i);
    }
}

int main(int argc, char** argv)
{

    size_t n = std::atoi(argv[1]);
    size_t m = std::atoi(argv[2]);
    double epsilon = std::atof(argv[3]);
    size_t k = std::atoi(argv[4]);
    double* matrix = new double[n * n];

    if (k == 0) {
        if (!fillFromFile(n, matrix, argv[5])) {
            std::cout << "file error." << std::endl;
            return -1;
        }
    } else {
        if (a[k] != nullptr) {
            fillFromFunction(n, matrix, a[k]);
        } else {
            std::cout << "function does not exist." << std::endl;
            return -1;
        }
    }

    print(n, m, matrix);

    double* x1 = new double[n];
    double* x2 = new double[n];
    double* eigenvalues = new double[n];

    double _inv1 = inv1(n, matrix);
    double _inv2 = inv2(n, matrix);

    auto t1 = std::chrono::high_resolution_clock::now();

    findEigenvalues(n, matrix, eigenvalues, epsilon, x1, x2);

    auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> tms = t2 - t1;
    std::cout << "time: " << tms.count() / 1000.0 << " sec." << std::endl;

    print(m, eigenvalues);

    for (size_t i = 0; i < n; ++i) {
        _inv1 += eigenvalues[i];
        _inv2 += eigenvalues[i] * eigenvalues[i];
    }

    std::cout << "|Sum(x_i) - Sum(a_i)| = " << std::abs(_inv1) << std::endl
              << "|Sum(x_i ^ 2) - Sum(a_ij * a_ji)| = " << std::abs(_inv2) << std::endl;

    delete[] matrix;
    delete[] x1;
    delete[] x2;
    delete[] eigenvalues;
    return 0;
}