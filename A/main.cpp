#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <cstring>

using namespace std;

// Uogólnić?
class SquareMatrix
{
    int n;
    int *data;

public:
    SquareMatrix(int n);

    SquareMatrix(const SquareMatrix &other);

    ~SquareMatrix();

    int &operator()(int x, int y)
    {
        return data[x * n + y];
    }

    const int &operator()(int x, int y) const
    {
        return const_cast<SquareMatrix *>(this)->operator()(x, y);
    }

    int width() const
    {
        return n;
    }
};

SquareMatrix::SquareMatrix(int n) : n(n)
{
    data = new int[n * n];
    //memset(data, 0, n * n * sizeof(int));
    for (int i = 0; i < n * n; ++i)
        data[i] = 0;
}

SquareMatrix::SquareMatrix(const SquareMatrix &other) : n(other.n)
{
    data = new int[n * n];
    for (int i = 0; i < n * n; ++i)
        data[i] = other.data[i];
}

SquareMatrix::~SquareMatrix()
{
    delete[] data;
}

static SquareMatrix operator*(const SquareMatrix &lhs, const SquareMatrix &rhs)
{
    const int width = lhs.width();
    assert(width == rhs.width());
    SquareMatrix result(width);
    // przeplot pętli dla obliczenia cache-friendly
    for (int i = 0; i < width; ++i) {
        for (int k = 0; k < width; ++k) {
            for (int j = 0; j < width; ++j) {
                result(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return result;
}

static SquareMatrix operator*(int lhs, const SquareMatrix &rhs)
{
    const int width = rhs.width();
    SquareMatrix result(width);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            result(i, j) = lhs * rhs(i, j);
        }
    }
    return result;
}

static SquareMatrix operator-(const SquareMatrix &lhs, const SquareMatrix &rhs)
{
    const int width = lhs.width();
    assert(width == rhs.width());
    SquareMatrix result(width);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            result(i, j) = lhs(i, j) - rhs(i, j);
        }
    }
    return result;
}

/******************* ROZWIĄZANIE *****************/

// Kwadrat grafu (wykład2 6/26) (macierz sąsiedztwa -> macierz sąsiedztwa)
// wierzchołki są połączone, jeżeli odległość między nimi jest nie większa niż 2
static SquareMatrix squaredGraph(const SquareMatrix &graph, const SquareMatrix &squaredGraph)
{
    const int width = graph.width();
    SquareMatrix result(width);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            if (i == j) {
                continue;
            }
            result(i, j) = static_cast<int>(graph(i, j) > 0 || squaredGraph(i, j) > 0);
        }
    }
    return result;
}

static SquareMatrix squaredGraph(const SquareMatrix &graph)
{
    return squaredGraph(graph, graph * graph);
}

static SquareMatrix ADP(const SquareMatrix &A)
{
    const int n = A.width();
    const SquareMatrix Z = A * A;
    const SquareMatrix Ap = squaredGraph(A, Z);
    bool finished = true;
    for (int i = 0; finished && (i < n); ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            if (Ap(i, j) != 1) {
                finished = false;
                break;
            }
        }
    }
    if (finished) {
        return 2 * Ap - A;
    } else {
        const SquareMatrix Dp = ADP(Ap);
        const SquareMatrix S = A * Dp;
        SquareMatrix D(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                D(i, j) = 2 * Dp(i, j) - static_cast<int>(S(i, j) < Dp(i, j) * Z(i, i));
            }
        }
        return D;
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    int n;
    cin >> n;
    SquareMatrix input(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> input(i, j);
        }
    }
    const SquareMatrix output = ADP(input);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << output(i, j) << ' ';
        }
        cout << '\n';
    }
    return 0;
}