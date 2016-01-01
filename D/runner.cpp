#ifndef DEBUG
#define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <map>
#include <set>

using namespace std;

int main()
{
    for (; ;) {
        if (system("./gen > killer.in")) {
            cout << "Generator zdechł\n";
            break;
        }
        int result = system("./solution < killer.in > killer.tout");
        if (result) {
            cout << "Wzorcówka zdechła\n";
            break;
        } else {
            cout << "OK?\n";
        }
        if (system("./verify killer.in killer.tout")) {
            cout << "Błędna odpowiedź\n";
            break;
        }
    }
    return 0;
}