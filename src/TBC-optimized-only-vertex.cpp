#include <omp.h>

#include "StatisticTools.h"

vector<int> duplicate;
vector<vector<int>> reducedTimes;

vector<vector<int>> temporalGraph;
vector<vector<int>> reversedTemporalGraph;

vector<int> vToCompressed, compressedTov;

vector<vector<int>> g;
vector<vector<int>> rg;
vector<int> parentOffset;

vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}

constexpr int MIN_TIME = 0;
constexpr int MAX_TIME = 1e9;

constexpr int INF = 1e9;

struct BfsData {
    vector<int> dist, order, marked, layer, nextLayer;
    vector<double> delta, sigma, weights;
    vector<int> parentCount;
    vector<int> parent;
    BfsData(int n, int64_t numEdges)
        : dist(n, INF),
          marked(n, 0),
          delta(n, 0),
          sigma(n, 0),
          weights(n, 0),
          parentCount(n, 0)
    {
        layer.reserve(n);
        nextLayer.reserve(n);
        parent.resize(numEdges);
        order.reserve(n);
    }

    void clear() {
        for (int x : order) {
            dist[x] = INF;
            marked[x] = 0;
            delta[x] = 0;
            sigma[x] = 0;
            weights[x] = 0;
            parentCount[x] = 0;
        }
        order.clear();
    }
    void addParent(int x, int y) {
        parent[parentOffset[x] + parentCount[x]] = y;
        parentCount[x]++;
    }
};


void bfs(int v, vector<double> &BC, BfsData& data) {
    data.clear();

    auto& dist = data.dist;
    auto& order = data.order;
    auto& marked = data.marked;
    auto& delta = data.delta;
    auto& sigma = data.sigma;
    auto& weights = data.weights;

    auto& layer = data.layer;
    auto& nextLayer = data.nextLayer;

    int startv = reducedTimes[v][0];
    dist[startv] = 0;
    sigma[startv] = 1;

    layer.push_back(startv);

    int nd = 1;

    while (!layer.empty()) {
        for (int x : layer) {
            order.push_back(x);

            double waysX = sigma[x] * duplicate[x];

            for (int y : g[x]) {
                if (dist[y] == INF) {
                    dist[y] = nd;
                    sigma[y] = waysX;

                    nextLayer.push_back(y);
                    data.addParent(y, x);
                }
                else if (dist[y] == nd) {
                    sigma[y] += waysX;
                    data.addParent(y, x);
                }
            }
        }

        layer.swap(nextLayer);
        nextLayer.clear();
        nd++;
    }
    for (int i = 0; i < sz(reducedTimes); i++) {
        if (i == v) continue;
        int mindist = INF;
        for (auto id : reducedTimes[i]) {
            mindist = min(mindist, dist[id]);
        }
        if (mindist == INF) continue;
        double w = 0;
        for (auto num : reducedTimes[i]) {
            if (dist[num] == mindist) {
                marked[num] = 1;
                w += sigma[num] * duplicate[num];
            }
        }
        for (auto num : reducedTimes[i]) {
            if (marked[num]) {
                weights[num] = w;
            }
        }
    }

    for (int i = sz(order) - 1; i >= 0; i--) {
        int x = order[i];

        double factor = delta[x] / sigma[x];

        if (marked[x]) {
            factor += double(duplicate[x]) / weights[x];
        }

        for (int j = 0; j < data.parentCount[x]; j++) {
            int y = data.parent[parentOffset[x] + j];
            delta[y] += double(duplicate[y]) * sigma[y] * factor;
        }
    }

    for (int i = 0; i < sz(reducedTimes); i++) {
        if (i == v) continue;
        for (auto num : reducedTimes[i]) {
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
    int64_t numEdges = 0;
    parentOffset.resize(sz(rg));
    parentOffset[0] = 0;
    for (int i = 1; i < sz(rg); i++) {
        parentOffset[i] = parentOffset[i - 1] + sz(rg[i - 1]);
    }
    for (int i = 0; i < sz(g); i++) {
        numEdges += sz(g[i]);
    }
#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();

        BfsData data(sz(g), numEdges);
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
        // u v t
        ss >> edge[0] >> edge[1] >> edge[2];
        if (edge[0] == edge[1]) continue;
        edge[2]++;
        id_vertex.push_back(edge[0]);
        id_vertex.push_back(edge[1]);
        id_time.push_back(edge[2]);
        edges.push_back(edge);
    }

    id_time.push_back(MIN_TIME);
    id_time.push_back(MAX_TIME);
    sort(all(id_vertex));
    sort(all(id_time));
    sort(all(edges));

    edges.resize(unique(all(edges)) -  edges.begin());
    id_vertex.resize(unique(all(id_vertex)) - id_vertex.begin());
    id_time.resize(unique(all(id_time)) - id_time.begin());

    for (auto& [x, y, t] : edges) {
        x = lower_bound(all(id_vertex), x) - id_vertex.begin();
        y = lower_bound(all(id_vertex), y) - id_vertex.begin();
        t = lower_bound(all(id_time), t) - id_time.begin();
    }

    int N = 0; // number vertex in Temporal Graph
    int n = sz(id_vertex);
    vertexAllTimes.resize(n);
    vector<vector<int>> times(n);
    for (int i = 0; i < n; i++) {
        times[i].push_back(MIN_TIME);
    }
    for (auto [x, y, t] : edges) {
        times[y].push_back(t);
    }
    for (int i = 0; i < n; i++) {
        sort(all(times[i]));
        times[i].resize(unique(all(times[i])) - times[i].begin());
        for (auto t : times[i]) {
            vertexAllTimes[i].push_back({t, N});
            N++;
        }
    }
    temporalGraph.resize(N);
    reversedTemporalGraph.resize(N);
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
            // cerr << x << ' ' << y << ' ' << t1 << '\n';
            temporalGraph[id_x].push_back(idy);
            numberEdges++;
        }
    }
    for (auto &vec : temporalGraph) {
        sort(all(vec));
        vec.resize(unique(all(vec)) - vec.begin());
    }
    for (int i = 0; i < sz(temporalGraph); i++) {
        for (auto y : temporalGraph[i]) {
            reversedTemporalGraph[y].push_back(i);
        }
    }
    int reduceN = 0;
    int reduceE = 0;
    reducedTimes.resize(n);
    vToCompressed.resize(N);

    for (int i = 0; i < n; i++) {
        // TODO add Hash-s
        map<pair<vector<int>, vector<int>>, vector<int>> eqls;
        for (auto [t, idx] : vertexAllTimes[i]) {
            if (t == MIN_TIME) {
                continue;
            }
            eqls[{temporalGraph[idx], reversedTemporalGraph[idx]}].push_back(idx);
        }
        int startId = vertexAllTimes[i][0].second;

        vToCompressed[startId] = reduceN;
        reducedTimes[i].push_back(reduceN);
        compressedTov.push_back(i);
        duplicate.push_back(1);
        reduceN++;

        for (auto& [gr, ids] : eqls) {
            for (int idx : ids) {
                vToCompressed[idx] = reduceN;
            }

            reducedTimes[i].push_back(reduceN);
            compressedTov.push_back(i);
            duplicate.push_back(sz(ids));

            reduceN++;
        }
    }
    g.resize(reduceN);
    rg.resize(reduceN);
    for (int i = 0; i < sz(temporalGraph); i++) {
        for (auto y : temporalGraph[i]) {
            if (vToCompressed[i] == vToCompressed[y]) continue;
            g[vToCompressed[i]].push_back(vToCompressed[y]);
        }
    }
    for (int i = 0; i < sz(g); i++) {
        sort(all(g[i]));
        g[i].resize(unique(all(g[i])) - g[i].begin());
        reduceE += sz(g[i]);
    }
    for (int i = 0; i < sz(g); i++) {
        for (auto y : g[i]) {
            rg[y].push_back(i);
        }
    }
    auto start = chrono::steady_clock::now();
    vector<double> BC = calcBC(true, threads);
    auto finish = chrono::steady_clock::now();
    double peakRAM = Statistic::getPeakRAM();
    string nameAlgo = "TBC-optimized-only-vertex";

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