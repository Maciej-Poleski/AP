#include <iostream>
#include <fstream>
#include <set>
#include <vector>

using namespace std;


int main(int argc, char **argv)
{
    if (argc != 3) {
        cerr << argv[0] << " [input] [output]\n";
        return 1;
    }
    ifstream input(argv[1]);
    ifstream output(argv[2]);
    int n, m;
    input >> n >> m;
    set<pair<int, int> > edges;
    for (int i = 0; i < m; ++i) {
        int a, b;
        input >> a >> b;
        edges.insert(make_pair(a, b));
        edges.insert(make_pair(b, a));
    }
    vector<bool> seen(n);
    for (int i = 0; i < n; i += 2) {
        int a, b;
        output >> a >> b;
        if (edges.find(make_pair(a, b)) == edges.end()) {
            cerr << "Krawędź {" << a << ", " << b << "} nie istnieje\n";
            return 2;
        }
        if (seen[a]) {
            cerr << "Wierzchołek " << a << " widziany\n";
            return 3;
        }
        if (seen[b]) {
            cerr << "Wierzchołek " << b << " widziany\n";
            return 3;
        }
    }
    return 0;
}