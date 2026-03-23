#include "Graph.h"

void processedDynamic(string filename) {

    ifstream in(filename);
    int n, q;
    in >> n >> q;
    cerr << n << ' ' << q << '\n';
    Graph g;
    g.assign(n);
    auto startTime = std::chrono::system_clock::now();
    for (int i = 0; i < q; i++) {
        int countOperation;
        in >> countOperation;
        for (int j = 0; j < countOperation; j++) {
            char operation;
            int x, y;
            in >> operation >> x >> y;
            x--; y--;
            if (operation == '+') {
                g.add_edge(x, y);
            } else {
                // operation == '-'
                g.remove_edge(x, y);
            }
        }
        cerr << g << '\n';
    }
    auto endTime = std::chrono::system_clock::now();
    long long int time = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();
    cerr << "filename : " << filename << endl;
    cerr << "time compute : " << time << " ms " << endl;
}

int main() {

    vector<string> filenameData = {
        "../Graphs/ParsedAs-733/as1.txt",
        "../Graphs/ParsedAs-733/as2.txt",
        "../Graphs/ParsedAs-733/as3.txt",
    };
    processedDynamic(filenameData[0]);

    return 0;
}