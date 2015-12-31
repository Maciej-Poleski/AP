#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <cstring>
#include <cmath>

using namespace std;

// Uogólnić?
class SquareMatrix
{
    int n;
    int *data;

public:
    explicit SquareMatrix(int n);

    SquareMatrix(const SquareMatrix &other);

    ~SquareMatrix();

    SquareMatrix &operator=(const SquareMatrix &o);

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

SquareMatrix &SquareMatrix::operator=(const SquareMatrix &o)
{
    SquareMatrix A(o);
    using std::swap;
    swap(data, A.data);
    swap(n, A.n);
    return *this;
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

// operator+= ?
static SquareMatrix &operator+=(SquareMatrix &lhs, const SquareMatrix &rhs)
{
    const int width = rhs.width();
    assert(lhs.width() == width);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            lhs(i, j) += rhs(i, j);
        }
    }
    return lhs;
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

static SquareMatrix operator-(const SquareMatrix &o)
{
    const int width = o.width();
    SquareMatrix result(width);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < width; ++j) {
            result(i, j) = -o(i, j);
        }
    }
    return result;
}

class IntMod
{
    int value;
    int mod;

public:
    IntMod(int value, int base);
};

IntMod::IntMod(int value, int base) : value(value), mod(base)
{ }

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

static int log2(int n) __attribute__((const));

static int log2(int n)
{
    int result = 0;
    for (; n > 1; n /= 2, ++result);
    return result;
}

static bool randRN(int r, int n)
{
    return (rand() % n) < r;
}

static SquareMatrix BPWM(const SquareMatrix &A, const SquareMatrix &B)
{
    SquareMatrix W = -(A * B);
    const int n = A.width();
    assert(n == B.width());
    for (int t = 0; t <= log2(n); ++t) {
        int r = 1 << t;
        int x = static_cast<int>(3.77 * log(n) / log(2) + 0.5);
        while (x--) {
            vector<int> R(n);
            for (int i = 0; i < n; ++i) {
                R[i] = randRN(r, n);
            }
            SquareMatrix Ar(n);
            for (int i = 0; i < n; ++i) {
                for (int k = 0; k < n; ++k) {
                    Ar(i, k) = R[k] * (k + 1) * A(i, k);
                }
            }
            SquareMatrix Z = Ar * B;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    int k = Z(i, j);
                    if ((k >= 1) && (k <= n) && (A(i, k - 1) != 0) && (B(k - 1, j) != 0)) {
                        W(i, j) = k;
                    }
                }
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (W(i, j) < 0) {
                for (int k = 0; k < n; ++k) {
                    if ((A(i, k - 1) != 0) && (B(k - 1, j) != 0)) {
                        W(i, j) = k;
                        break;
                    }
                }
            }
        }
    }
#ifdef DEBUG
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (W(i, j) != 0) {
                assert((A(i, W(i, j) - 1) != 0) && (B(W(i, j) - 1, j) != 0));
            } else {
                for (int k = 0; k < n; ++k) {
                    assert(A(i, k) * B(k, j) == 0);
                }
            }
#endif
    return W;
}

static SquareMatrix APSP(const SquareMatrix &A)
{
    const SquareMatrix D = ADP(A);
    int n = D.width();
    SquareMatrix S(n);
    for (int s = 0; s < 3; ++s) {
        SquareMatrix Ds(n);
        SquareMatrix Bs(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Ds(i, j) = ((D(i, j) % 3) == s);
                Bs(i, j) = (((D(i, j) + 1) % 3) == s);
            }
        }
        //Ws = świadkowie dla Bs i A
        SquareMatrix Ws = BPWM(Bs, A);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (Ds(i, j) != 1) {
                    Ws(i, j) = 0;
                }
            }
        }
        S += Ws;
    }
#ifdef DEBUG
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            if (D(i, j) != D(i, S(i, j) - 1) + 1) {
                cerr << i << " " << j << " S: " << S(i, j) << "\n";
            }
        }
#endif
    return S;
}

// find num^(-1) in Z_{base} (base is prime)
static int inverted(int num, int base)
{
    const long oldNum = num;
    const int oldBase = base;
    assert(num > 0);
    assert(base > 0);
    int p = 1, q = 0;
    int r = 0, s = 1;
    while (base != 0) {
        int c = num % base;
        int quot = num / base;
        num = base;
        base = c;
        int rTmp = r, sTmp = s;
        r = p - quot * r;
        s = q - quot * s;
        p = rTmp;
        q = sTmp;
    }
    while (p < 0) {
        p += oldBase;
    }
    assert(p * oldNum % oldBase == 1);
    return p;
}

static bool det(SquareMatrix matrix, int base)
{
    const int n = matrix.width();
    for (int i = 0; i < n - 1; ++i) {
        // i - wiersz pozostawiony
        // 1: ustlić czy dobry wybór
        if (matrix(i, i) == 0) {
            bool ok = false;
            for (int j = i + 1; j < n; ++j) {
                if (matrix(j, i) != 0) {
                    ok = true;
#ifdef DEBUG
                    for (int k = 0; k < i; ++k) {
                        assert(matrix(i, k) == 0);
                        assert(matrix(j, k) == 0);
                    }
#endif
                    for (int k = i; k < n; ++k) {
                        using std::swap;
                        swap(matrix(i, k), matrix(j, k));
                    }
                    break;
                }
            }
            if (!ok) {
                return false;
                continue; // wszystkie kolejne wiersze mają 0 w tej kolumnie - nic do robienia
            }
        }
        assert(matrix(i, i) != 0);
        for (int j = i + 1; j < n; ++j) {
            //j - wiersz aktywny
            //1: ustalić mnożnik
            if (matrix(j, i) == 0) {
                continue;
            }
            assert(matrix(i, i) > 0);
            assert(matrix(i, i) < base);
            assert(matrix(j, i) > 0);
            assert(matrix(j, i) < base);
            const long coeff = (base - (static_cast<long>(inverted(matrix(i, i), base)) * matrix(j, i)) % base) % base;
            assert(coeff > 0);
            assert(coeff < base);
            for (int k = i; k < n; ++k) {
                assert(matrix(i, k) >= 0);
                matrix(j, k) = (matrix(j, k) + coeff * matrix(i, k)) % base;
            }
            assert(matrix(j, i) == 0);
        }
    }
    return true;
}

// Znajdź pierwszy wiersz który ma nie 0 w danej kolumnie (przy założeniu że taki istnieje)
static int findRow(const SquareMatrix &matrix, int column, int startRow)
{
    int i;
    for (i = startRow; matrix(i, column) == 0; ++i);
    return i;
}

static void swapRows(SquareMatrix &leftPart, SquareMatrix &rightPart, const int startColumnt, const int i, const int j)
{
    if (i == j) {
        return;
    }
    const int n = leftPart.width();
    assert(n == rightPart.width());
    using std::swap;
#ifdef DEBUG
    for (int k = 0; k < startColumnt; ++k) {
        assert(leftPart(i, k) == 0);
        assert(leftPart(j, k) == 0);
    }
#endif
    for (int k = startColumnt; k < n; ++k)
        swap(leftPart(i, k), leftPart(j, k));
    for (int k = 0; k < n; ++k) {
        swap(rightPart(i, k), rightPart(j, k));
    }
}

static void multiplyRow(SquareMatrix &leftPart, SquareMatrix &rightPart, const int startColumn, int row, int multiply_,
                        int base)
{
    const long multiply = multiply_;
    const int n = leftPart.width();
    assert(n == rightPart.width());
    assert(multiply > 0);
    assert(multiply < base);
#ifdef DEBUG
    for (int i = 0; i < startColumn; ++i) {
        assert(leftPart(row, i) == 0);
    }
#endif
    for (int i = startColumn; i < n; ++i) {
        leftPart(row, i) = (leftPart(row, i) * multiply) % base;
    }
    for (int i = 0; i < n; ++i) {
        rightPart(row, i) = (rightPart(row, i) * multiply) % base;
    }
}

static void addRowMultiplied(SquareMatrix &leftPart, SquareMatrix &rightPart, int fromRow, int toRow, int startColumn,
                             int multiply_, int base)
{
    const long multiply = multiply_;
    const int n = leftPart.width();
    assert(rightPart.width() == n);
    assert(multiply > 0);
    assert(multiply < base);
#ifdef DEBUG
    // tylko wiersz źródłowy ma oczyszczoną lewą strone
    for (int i = 0; i < startColumn; ++i) {
        assert(leftPart(fromRow, i) == 0);
    }
#endif
    for (int i = startColumn; i < n; ++i) {
        leftPart(toRow, i) = (leftPart(toRow, i) + leftPart(fromRow, i) * multiply) % base;
    }
    for (int i = 0; i < n; ++i) {
        rightPart(toRow, i) = (rightPart(toRow, i) + rightPart(fromRow, i) * multiply) % base;
    }
}

static SquareMatrix inverted(SquareMatrix matrix, const int base)
{
    const int n = matrix.width();
    SquareMatrix result(n);
    for (int i = 0; i < n; ++i) {
        result(i, i) = 1;
    }
    for (int column = 0; column < n; ++column) {
        int selectedRow = findRow(matrix, column, column);
        swapRows(matrix, result, column, column, selectedRow);
#ifdef DEBUG
        for (int i = 0; i < column; ++i) {
            assert(matrix(column, i) == 0); // początek wiersza column = 0
        }
#endif
        assert(matrix(column, column) != 0);
        const int invCoeff = inverted(matrix(column, column), base);
        assert(invCoeff > 0);
        multiplyRow(matrix, result, column, column, invCoeff, base);
        for (int row = 0; row < n; ++row) {
            if (row == column) {
                continue;
            }
            if (matrix(row, column) != 0) {
                addRowMultiplied(matrix, result, column, row, column, base - matrix(row, column), base);
            }
        }
    }
#ifdef DEBUG
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            assert(result(i, j) >= 0);
            assert(result(i, j) < base);
        }
    }
#endif
    return result;
}

static const int bigPrime = 512000009;

static void removeCross(const SquareMatrix &A, const SquareMatrix &Ainv, int i, int j, int base, SquareMatrix &NA,
                        SquareMatrix &NAinv)
{
    const int n = A.width();
    assert(n == Ainv.width());
    assert(i >= 0);
    assert(j >= 0);
    assert(i < n);
    assert(j < n);
    SquareMatrix B(n - 1);
    SquareMatrix Ap(n - 1);
    int xDisp = 0;
    for (int x = 0; x < n; ++x) {
        if (x == i) {
            xDisp = 1;
            continue;
        }
        int yDisp = 0;
        for (int y = 0; y < n; ++y) {
            if (y == j) {
                yDisp = 1;
                continue;
            }
            B(x - xDisp, y - yDisp) = Ainv(x, y);
            Ap(x - xDisp, y - yDisp) = A(x, y);
        }
    }
    long factor = 0;
    xDisp = 0;
    int yDisp = 0;
    for (int k = 0; k < n - 1; ++k) {
        if (k == i) {
            xDisp = 1;
        }
        if (k == j) {
            yDisp = 1;
        }
        assert(Ainv(k + xDisp, j) >= 0);
        assert(Ainv(i, k + yDisp) >= 0);
        factor = (factor + Ainv(k + xDisp, j) * Ainv(i, k + yDisp)) % base;
        assert(factor >= 0);
    }
    factor = (factor * inverted(Ainv(i, j), base)) % base;
    assert(factor >= 0);
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            B(i, j) = (B(i, j) + base - factor) % base;
        }
    }
    NAinv = B;
    NA = Ap;
}

static vector<bool> removed;

static void printOutput(const int i, const int j)
{
    int iOff = 0;
    for (int x = i; x > 0 || removed[iOff]; ++iOff) {
        if (!removed[iOff]) {
            x -= 1;
        }
    }
    int jOff = 0;
    for (int x = j; x > 0 || removed[jOff]; ++jOff) {
        if (!removed[jOff]) {
            x -= 1;
        }
    }
    removed[iOff] = true;
    removed[jOff] = true;
    cout << iOff << ' ' << jOff << '\n';
}

static void SimplePerfectMatching(SquareMatrix A, int base)
{
    SquareMatrix Ainv = inverted(A, base);
    for (; ;) {
        const int n = A.width();
        assert(n == Ainv.width());
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                if (A(i, j) != 0 && Ainv(i, j) != 0) {
                    printOutput(i, j);
                    goto step;
                }
            }
        }
        step:
        if (n == 2) {
            break;
        }
        assert(n > 2);
        SquareMatrix NA(0);
        SquareMatrix NAInv(0);
        if (i < j) {
            removeCross(A, Ainv, i, j, base, NA, NAInv);
            removeCross(NA, NAInv, j - 1, i, base, A, Ainv);
        } else {
            assert(i > j);
            removeCross(A, Ainv, i, j, base, NA, NAInv);
            removeCross(NA, NAInv, j, i - 1, base, A, Ainv);
        }
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    int n, m;
    cin >> n >> m;
    SquareMatrix input(n);
    removed.resize(n, false);
    srand(time(0));
    for (int i = 0; i < m; ++i) {
        int a, b;
        cin >> a >> b;
        int x = (rand() % (bigPrime - 1)) + 1;
        if (a < b) {
            input(a, b) = x;
            input(b, a) = bigPrime - x;
        } else if (a > b) {
            input(a, b) = bigPrime - x;
            input(b, a) = x;
        }
    }
    SimplePerfectMatching(input, bigPrime);
    return 0;
}