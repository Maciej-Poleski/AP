//#define DEBUG
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
#include <set>

using namespace std;

/******************* ROZWIĄZANIE *****************/

static int randN(int n)
{
    return rand() % n;
}

int main()
{
    ios_base::sync_with_stdio(false);
    srand(time(0));
    int z;
    cin >> z;
    while (z--) {
        int n, k;
        cin >> n >> k;
        vector<int> input(n);
        for (int i = 0; i < n; ++i) {
            cin >> input[i];
        }
        if (n < 100) {
            // przypadek brzegowy - nie chce mi się nad tym zastanawiać
            sort(input.begin(), input.end());
            cout << input[k] << '\n';
            continue;
        }
        for (bool failed = true; failed;) {
            int reducedSize = static_cast<int>(pow(n, 0.75) + .5);
            vector<int> reduced;
            set<int> selected;
            for (int i = 0; i < reducedSize; ++i) {
                int selection;
                do {
                    selection = randN(n);
                } while (selected.find(selection) != selected.end());
                selected.insert(selection);
                reduced.push_back(input[selection]);
            }
            assert(reducedSize == reduced.size());
            sort(reduced.begin(), reduced.end());
            double sqrtN = sqrt(n);
            double position = reducedSize * (static_cast<double>(k) / n);
            double left = position - sqrtN;
            int begin = left > 0 ? static_cast<int>(left) : 0;
            double right = position + sqrtN;
            int end = (right + .5) < reducedSize ? static_cast<int>(right + .5) : reducedSize - 1;
            assert(0 < begin);
            assert(begin < end);
            assert(end < n);
            vector<int> extracted;
            int less = 0;
            for (int i = 0; i < n; ++i) {
                if (input[i] < reduced[begin]) {
                    less += 1;
                } else if (input[i] <= reduced[end]) {
                    extracted.push_back(input[i]);
                }
            }
            if ((k >= less) && (k - less < extracted.size())) {
                sort(extracted.begin(), extracted.end());
                cout << extracted[k - less] << '\n';
                failed = false;
            }
        }
    }
    return 0;
}