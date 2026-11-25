#include <omp.h>

#include "StatisticTools.h"

vector<vector<int>> temporalGraph;
vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}
constexpr int MIN_TIME = 0;
constexpr int INF = 1e9;


struct BfsData {
    vector<int> dist, order, marked, blocked;
    vector<double> delta, sigma, weights;
    vector<vector<int>> parents;

    BfsData(int n)
        : dist(n, INF),
          marked(n, 0),
          blocked(n, 0),
          delta(n, 0),
          sigma(n, 0),
          weights(n, 0),
          parents(n)
    {
        order.reserve(n);
    }

    void clear() {
        fill(dist.begin(), dist.end(), INF);
        fill(marked.begin(), marked.end(), 0);
        fill(blocked.begin(), blocked.end(), 0);
        fill(delta.begin(), delta.end(), 0.0);
        fill(sigma.begin(), sigma.end(), 0.0);
        fill(weights.begin(), weights.end(), 0.0);

        order.clear();

        for (auto& p : parents) {
            p.clear();
        }
    }
};


void bfs(int v, vector<double> &BC, BfsData &data) {
    data.clear();

    auto& dist = data.dist;
    auto& order = data.order;
    auto& marked = data.marked;
    auto& blocked = data.blocked;
    auto& delta = data.delta;
    auto& sigma = data.sigma;
    auto& weights = data.weights;
    auto& parents = data.parents;


    int n = sz(temporalGraph);

    int startv = -1;
    for (auto [t1, idv] : vertexAllTimes[v]) {
        if (t1 == MIN_TIME) {
            startv = idv;
        } else {
            blocked[idv] = 1;
        }
    }
    dist[startv] = 0;
    delta[startv] = 0;
    sigma[startv] = 1;

    queue<int> q;
    q.push(startv);
    while (!q.empty()) {
        int x = q.front();
        order.push_back(x);
        q.pop();
        for (auto y : temporalGraph[x]) {
            if (blocked[y]) continue;
            if (dist[y] > dist[x] + 1) {
                dist[y] = dist[x] + 1;
                parents[y].clear();
                sigma[y] = sigma[x];
                q.push(y);
            }
            else if (dist[y] == dist[x] + 1) {
                sigma[y] += sigma[x];
            }
            if (dist[y] == dist[x] + 1) {
                parents[y].push_back(x);
            }
        }
    }
    for (int i = 0; i < sz(vertexAllTimes); i++) {
        if (i == v) continue;
        int mindist = INF;
        for (auto [t, num] : vertexAllTimes[i]) {
            mindist = min(mindist, dist[num]);
        }
        if (mindist == INF) continue;
        double w = 0;
        for (auto [t, num] : vertexAllTimes[i]) {
            if (dist[num] == mindist) {
                marked[num] = 1;
                w += sigma[num];
            }
        }
        for (auto [t, num] : vertexAllTimes[i]) {
            if (marked[num]) {
                weights[num] = w;
            }
        }
    }

    for (int i = sz(order) - 1; i >= 0; i--) {
        int x = order[i];
        for (auto y : parents[x]) {
            delta[y] += (sigma[y] / sigma[x]) * delta[x];
            if (marked[x]) {
                delta[y] += sigma[y] / weights[x];
            }
        }
    }

    for (int i = 0; i < sz(vertexAllTimes); i++) {
        if (i == v) continue;
        for (auto [t, num] : vertexAllTimes[i]) {
            if (dist[num] != INF)
                BC[i] += delta[num];
        }
    }
}

vector<double> calcBC(bool normalize = false, int threads = 1) {
    int n = sz(vertexAllTimes);
    vector<vector<double>> partialBC(
        threads, vector<double>(n, 0.0)
    );

#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();
        BfsData data(sz(temporalGraph));
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; i++) {
            bfs(i, partialBC[tid], data);
        }
    }

    vector<double> BC(n, 0.0);

    for (int tid = 0; tid < threads; tid++) {
        for (int i = 0; i < n; i++) {
            BC[i] += partialBC[tid][i];
        }
    }
    if (normalize) {
        for (int i = 0; i < n; i++) {
            double coef = 1.0 / ((double)(n) * (n - 1));
            BC[i] *= coef;
        }
    }
    return BC;
}

void processedDynamic(string filename, int threads = 1) {
    ifstream in(filename);
    string line;
    getline(in, line);
    vector<array<int, 3>> edges;
    vector<int> id_vertex;
    vector<int> id_time;
    while (getline(in,line)) {
        stringstream ss(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        array<int, 3> edge;
        ss >> edge[0] >> edge[1] >> edge[2];
        if (edge[0] == edge[1]) continue;
        edge[2]++;
        id_vertex.push_back(edge[0]);
        id_vertex.push_back(edge[1]);
        id_time.push_back(edge[2]);
        edges.push_back(edge);
    }
    id_time.push_back(MIN_TIME);

    sort(all(id_vertex));
    sort(all(id_time));
    sort(all(edges));

    edges.resize(unique(all(edges)) -  edges.begin());
    id_vertex.resize(unique(all(id_vertex)) - id_vertex.begin());
    id_time.resize(unique(all(id_time)) - id_time.begin());

    vector<int> min_time_in(sz(id_vertex), INF);
    for (auto& [x, y, t] : edges) {
        x = lower_bound(all(id_vertex), x) - id_vertex.begin();
        y = lower_bound(all(id_vertex), y) - id_vertex.begin();
        t = lower_bound(all(id_time), t) - id_time.begin();
        min_time_in[y] = min(min_time_in[y], t);
    }

    int N = 0; // number vertex in Temporal Graph
    int n = sz(id_vertex);
    vertexAllTimes.resize(n);
    vector<set<int>> distinctTimes(n);
    for (int i = 0; i < n; i++) {
        vertexAllTimes[i].push_back({MIN_TIME, N});
        distinctTimes[i].insert(MIN_TIME);
        N++;
    }
    for (auto [x, y, t] : edges) {
        if (distinctTimes[y].count(t) > 0) continue;
        distinctTimes[y].insert(t);
        vertexAllTimes[y].push_back({t, N});
        N++;
    }
    temporalGraph.resize(N);
    int numberEdges = 0;
    for (auto [x, y, t] : edges) {
        // проводим для всех вершин ребра
        int idy = -1;
        for (auto [t1, id_y] : vertexAllTimes[y]) {
            if (t == t1) {
                idy = id_y;
                break;
            }
        }
        for (auto [t1, id_x] : vertexAllTimes[x]) {
            if (t1 >= t) continue;
            temporalGraph[id_x].push_back(idy);
            numberEdges++;
        }
    }

    auto start = chrono::steady_clock::now();
    vector<double> BC = calcBC(true, threads);
    auto finish = chrono::steady_clock::now();
    double peakRAM = Statistic::getPeakRAM();
    // cerr << "Total visited : " << totalVisited << endl;
    // cerr << "Total scanned Edge: " << totalScannedEdges << endl;
    // cerr << "Total parent edges : " << totalParentEdges << endl;
    string nameAlgo = "TBC-base";
    Statistic::getAllStatistic(filename, nameAlgo, threads, peakRAM, std::chrono::duration<double>(finish - start).count());
}

int main(int argc, char* argv[]) {
    string filename;
    int threads = 1;
    if (argc < 3) {
        threads = 4;
        filename = "../Graphs/Temporary/1_06_college_msg.patg";
    } else {
        threads = atoi(argv[1]);
        filename = argv[2];
    }
    processedDynamic(filename, threads);
    return 0;
}