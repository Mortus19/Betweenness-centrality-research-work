#include <bits/stdc++.h>

#include "../include/matrix.h"

using namespace std;
using ll = long long;
#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

struct Graph {
    int n, m;
    vector<vector<int>> g, parents;
    vector<double> betweenness, sigma, delta;
    vector<int> dist, order;
    const string graphPath;
    string graphName, format;

    Graph(const string &graphPath) : graphPath(graphPath) {
        int slashPos = sz(graphPath) - 1;
        int dotPos   = sz(graphPath) - 1;
        while (graphPath[slashPos] != '/')
            --slashPos;
        while (graphPath[dotPos] != '.')
            --dotPos;
        graphName = graphPath.substr(slashPos + 1, dotPos - slashPos - 1);
        format = graphPath.substr(dotPos + 1, graphPath.size() - dotPos);

        readGraph<double>();
    }
    void assign(int vertexCount) {
        n = vertexCount;
        g.clear();
        g.resize(n);
        m = 0;
    }

    void add_edge(int u, int v) {
        g[u].push_back(v);
        m++;
    }

    void bfs(int s) {
        dist.assign(n, -1);
        order.clear();
        order.reserve(n);
        parents.clear(); parents.resize(n);
        sigma.assign(n, 0.0); delta.assign(n, 0.0);
        sigma[s] = 1.0;
        dist[s] = 0;
        queue<int> q;
        q.push(s);
        while (!q.empty()) {
            auto x = q.front();
            q.pop();
            order.push_back(x);
            for (auto &y : g[x]) {
                if (dist[y] == -1) {
                    dist[y] = dist[x] + 1;
                    q.push(y);
                }
                if (dist[y] == dist[x] + 1) {
                    sigma[y] += sigma[x];
                    parents[y].push_back(x);
                }
            }
        }
        for (int i = sz(order) - 1; i >= 0; i--) {
            int x = order[i];
            for (auto &y : parents[x]) {
                delta[y] += (sigma[y] / sigma[x]) * (1 + delta[x]);
            }
            if (x != s) {
                betweenness[x] += delta[x];
            }
        }
    }

    auto computeBetweennessMs() {
        auto startTime = std::chrono::system_clock::now();
        betweenness.assign(n, 0);
        for (int i = 0; i < n; i++) {
            bfs(i);
        }
        auto endTime = std::chrono::system_clock::now();
        return chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count();;
    }

    template <typename ValType>
    void readGraph() {
        spMtx<ValType> gr(graphPath.c_str(), format.c_str());
        assign(gr.n);
        for (int u = 0; u < gr.n; ++u) {
            for (int id = gr.Rst[u]; id < gr.Rst[u + 1]; ++id) {
                int v = gr.Col[id];
                add_edge(u, v);
            }
        }
    }

    friend ostream& operator<<(ostream& os, Graph& g) {
        os << "Graph name: " << g.graphName << endl;
        os << "Vertices: " << g.n << endl;
        os << "Edges: " << g.m << endl;
        os << "Time Betweenness: " << g.computeBetweennessMs() << " ms" << endl;
        return os;
    }
};


int main() {
    vector<string> graphPaths = {
        "../Graphs/cage3.mtx",
        "../Graphs/G43.mtx",
        "../Graphs/iris_dataset_30NN.mtx",
        "../Graphs/yeast_30NN.mtx",
    };
    for (auto &graphPath : graphPaths) {
        Graph g(graphPath);
        cerr << g << '\n';
    }
    return 0;
}