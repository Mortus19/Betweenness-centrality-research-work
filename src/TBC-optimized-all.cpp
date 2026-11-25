#include <omp.h>

#include  "StatisticTools.h"

// Reduce simple segment time [l, r]
struct vertexInfo {
    int instance, timeL;

    vertexInfo(int _intstance = 1, int _l = 0) {
        instance = _intstance;
        timeL = _l;
    }
};

struct Edge {
    int to, w;
    Edge(int _to = 0, int _w = 1) : to(_to), w(_w) {}
};

vector<vertexInfo> vertexData;
vector<vector<int>> reducedTimes;

vector<vector<int>> temporalGraph;
vector<vector<int>> reversedTemporalGraph;

vector<vector<Edge>> g;

vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}

constexpr int MIN_TIME = 0;
constexpr int MAX_TIME = 1e9;
constexpr int INF = 1e9;
vector<int> zeroTopoOrder;
vector<int> csrOffset, parentOffset;
vector<Edge> csrEdges;

void buildZeroTopoOrder() {
    int n = sz(g);

    vector<int> indeg(n, 0);
    vector<vector<int>> zeroGraph(n);

    for (int u = 0; u < n; u++) {
        for (auto e : g[u]) {
            if (e.w == 0) {
                zeroGraph[u].push_back(e.to);
                indeg[e.to]++;
            }
        }
    }

    queue<int> q;

    for (int i = 0; i < n; i++) {
        if (indeg[i] == 0) {
            q.push(i);
        }
    }

    zeroTopoOrder.clear();

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        zeroTopoOrder.push_back(v);

        for (int to : zeroGraph[v]) {
            indeg[to]--;

            if (indeg[to] == 0) {
                q.push(to);
            }
        }
    }
    if (sz(zeroTopoOrder) != n) {
        throw runtime_error("The zero-weight subgraph is not acyclic");
    }
}


// Prepare contiguous adjacency and a fixed slot range for each vertex's predecessors.
void buildCSR() {
    const int n = sz(g);
    csrOffset.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        csrOffset[i + 1] = csrOffset[i] + sz(g[i]);
    }
    csrEdges.resize(csrOffset[n]);
    for (int i = 0; i < n; ++i) {
        copy(g[i].begin(), g[i].end(), csrEdges.begin() + csrOffset[i]);
    }

    parentOffset.assign(n + 1, 0);
    for (const Edge& e : csrEdges) {
        ++parentOffset[e.to + 1];
    }
    for (int i = 0; i < n; ++i) {
        parentOffset[i + 1] += parentOffset[i];
    }

    vector<vector<Edge>>().swap(g);
}

struct BfsData {
    vector<int> dist, marked, blocked, order, visited, blockedNodes;
    vector<int> layer, nextLayer, parentCount, parent;
    vector<int> positions;
    vector<double> sigma, endSigma, weights, propagate, bcLocal;

    BfsData(int n, size_t edges)
        : dist(n, INF), marked(n, 0), blocked(n, 0),
          parentCount(n, 0), parent(edges), positions(n + 1),
          sigma(n, 0), endSigma(n, 0), weights(n, 0),
          propagate(n, 0), bcLocal(n, 0) {
        order.reserve(n);
        visited.reserve(n);
        layer.reserve(n);
        nextLayer.reserve(n);
        blockedNodes.reserve(n);
    }

    void clear() {
        for (int x : visited) {
            dist[x] = INF;
            marked[x] = 0;
            sigma[x] = 0;
            endSigma[x] = 0;
            weights[x] = 0;
            propagate[x] = 0;
            bcLocal[x] = 0;
            parentCount[x] = 0;
        }
        for (int x : blockedNodes) blocked[x] = 0;
        visited.clear();
        blockedNodes.clear();
        order.clear();
        layer.clear();
        nextLayer.clear();
    }
};

void bfs(int v, vector<double>& BC, BfsData& data) {
    data.clear();

    auto& dist = data.dist;
    auto& marked = data.marked;
    auto& blocked = data.blocked;
    auto& sigma = data.sigma;
    auto& endSigma = data.endSigma;
    auto& weights = data.weights;
    auto& propagate = data.propagate;
    auto& bcLocal = data.bcLocal;
    auto& order = data.order;
    auto& layer = data.layer;
    auto& nextLayer = data.nextLayer;

    int startv = -1;
    for (int id : reducedTimes[v]) {
        if (vertexData[id].timeL == MIN_TIME) {
            startv = id;
        } else {
            blocked[id] = 1;
            data.blockedNodes.push_back(id);
        }
    }

    dist[startv] = 0;
    sigma[startv] = 1;
    endSigma[startv] = 1;
    data.visited.push_back(startv);
    layer.push_back(startv);

    // 0/1 BFS by layers: zero-weight edges remain in this layer;
    // unit-weight edges go to the next one. Stale entries are skipped.
    int depth = 0;
    while (!layer.empty()) {
        for (size_t pos = 0; pos < layer.size(); ++pos) {
            const int x = layer[pos];
            if (dist[x] != depth) continue;

            for (int j = csrOffset[x]; j < csrOffset[x + 1]; ++j) {
                const Edge e = csrEdges[j];
                const int y = e.to;
                if (blocked[y]) continue;
                const int nd = depth + e.w;
                if (nd < dist[y]) {
                    if (dist[y] == INF) data.visited.push_back(y);
                    dist[y] = nd;
                    if (e.w == 0) layer.push_back(y);
                    else nextLayer.push_back(y);
                }
            }
        }
        layer.swap(nextLayer);
        nextLayer.clear();
        ++depth;
    }

    // Counting sort by distance, retaining zero-edge topological order
    // within each distance layer.
    int maxDist = 0;
    for (int x : data.visited) maxDist = max(maxDist, dist[x]);
    fill(data.positions.begin(), data.positions.begin() + maxDist + 2, 0);
    for (int x : data.visited) ++data.positions[dist[x] + 1];
    for (int d = 1; d <= maxDist + 1; ++d) {
        data.positions[d] += data.positions[d - 1];
    }
    order.resize(data.visited.size());
    for (int x : zeroTopoOrder) {
        if (dist[x] != INF) {
            order[data.positions[dist[x]]++] = x;
        }
    }

    // Dynamic programming for shortest-path multiplicities and endpoints;
    // parents use preallocated contiguous ranges, one range per destination.
    for (int x : order) {
        for (int j = csrOffset[x]; j < csrOffset[x + 1]; ++j) {
            const Edge e = csrEdges[j];
            const int y = e.to;
            if (blocked[y] || dist[y] != dist[x] + e.w) continue;

            data.parent[parentOffset[y] + data.parentCount[y]++] = x;
            const double ways = sigma[x] * (e.w ? vertexData[y].instance : 1);
            sigma[y] += ways;
            if (e.w) endSigma[y] += ways;
        }
    }

    for (int i = 0; i < sz(reducedTimes); ++i) {
        if (i == v) continue;
        int mindist = INF;
        for (int id : reducedTimes[i]) mindist = min(mindist, dist[id]);
        if (mindist == INF) continue;

        double totalWays = 0;
        for (int id : reducedTimes[i]) {
            if (dist[id] == mindist) {
                marked[id] = 1;
                totalWays += endSigma[id];
            }
        }
        for (int id : reducedTimes[i]) {
            if (marked[id]) weights[id] = totalWays;
        }
    }

    for (int i = sz(order) - 1; i >= 0; --i) {
        const int x = order[i];
        const double baseFactor = sigma[x] > 0 ? propagate[x] / sigma[x] : 0.0;
        const double endFactor = baseFactor + (marked[x] ? 1.0 / weights[x] : 0.0);
        const int begin = parentOffset[x];
        const int end = begin + data.parentCount[x];

        for (int j = begin; j < end; ++j) {
            const int y = data.parent[j];
            const int w = dist[x] - dist[y];
            const double contrib = sigma[y] * (w ? vertexData[x].instance : 1)
                                 * (w ? endFactor : baseFactor);
            propagate[y] += contrib;
            if (w) bcLocal[y] += contrib;
        }
    }

    for (int i = 0; i < sz(reducedTimes); ++i) {
        if (i == v) continue;
        for (int id : reducedTimes[i]) {
            if (dist[id] != INF) BC[i] += bcLocal[id];
        }
    }
}

vector<double> calcBC(bool normalize = false, int threads = 1) {
    const int n = sz(vertexAllTimes);
    vector<vector<double>> partialBC(threads, vector<double>(n, 0.0));

    #pragma omp parallel num_threads(threads)
    {
        const int tid = omp_get_thread_num();
        BfsData data(sz(csrOffset) - 1, csrEdges.size());
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; ++i) {
            bfs(i, partialBC[tid], data);
        }
    }

    vector<double> BC(n, 0.0);
    for (int t = 0; t < threads; ++t) {
        for (int i = 0; i < n; ++i) BC[i] += partialBC[t][i];
    }
    if (normalize) {
        const double coef = 1.0 / (double(n) * (n - 1));
        for (double& x : BC) x *= coef;
    }
    return BC;
}

void processedDynamic(string filename, int threads = 1, string output = "") {
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
            temporalGraph[id_x].push_back(idy);
            reversedTemporalGraph[idy].push_back(id_x);
            numberEdges++;
        }
    }
    for (auto &vec : temporalGraph) {
        sort(all(vec));
        vec.resize(unique(all(vec)) - vec.begin());
    }
    for (auto &vec : reversedTemporalGraph) {
        sort(all(vec));
        vec.resize(unique(all(vec)) - vec.begin());
    }

    int reduceN = 0;
    int reduceE = 0;
    reducedTimes.resize(n);
    map<pair<int, int>, int> index; // <number, time> -> vertex id
    map<int, vector<int>> N_IN; // vertex id -> IN
    map<pair<int, int>, int> root_virtual_edge_vertex; // <number, time> -> vertex root
    for (int i = 0; i < n; i++) {
        map<pair<vector<int>, vector<int>>, vector<int>> eqls;
        for (auto [t, idx] : vertexAllTimes[i]) {
            if (t == MIN_TIME) continue;
            eqls[{temporalGraph[idx], reversedTemporalGraph[idx]}].push_back(t);
        }
        index[{i, MIN_TIME}] = reduceN;
        vertexData.push_back({1, MIN_TIME});
        reducedTimes[i].push_back(reduceN);
        g.push_back({});
        reduceN++;
        map<vector<int>, int> lst;
        for (auto [gr, ts] : eqls) {
            for (auto t : ts) {
                index[{i, t}] = reduceN;
            }
            N_IN[reduceN] = gr.second;
            vertexData.push_back(vertexInfo(sz(ts), ts[0]));
            reducedTimes[i].push_back(reduceN);
            g.push_back({});
            reduceN++;
        }
        sort(all(reducedTimes[i]), [&](int a, int b) {
            if (vertexData[a].timeL == vertexData[b].timeL) {
                return a < b;
            }
            return vertexData[a].timeL < vertexData[b].timeL;
        });
        // virtual edges
        for (auto id : reducedTimes[i]) {
            auto &in = N_IN[id];
            if (lst.count(in) > 0) {
                reduceE++;
                g[lst[in]].push_back(Edge(id, 0));
                root_virtual_edge_vertex[{i, id}] = root_virtual_edge_vertex[{i, lst[in]}];
            } else {
                root_virtual_edge_vertex[{i, id}] = id;
            }
            lst[in] = id;
        }
    }

    map<pair<int, int>, int> bestFrom;

    for (auto [x, y, t] : edges) {
        int target = index[{y, t}];

        for (int i = sz(reducedTimes[x]) - 1; i >= 0; i--) {
            int from = reducedTimes[x][i];

            if (vertexData[from].timeL >= t) continue;

            int sourceRoot = root_virtual_edge_vertex[{x, from}];

            pair<int, int> key = {sourceRoot, target};

            if (!bestFrom.count(key) || vertexData[from].timeL > vertexData[bestFrom[key]].timeL) {
                bestFrom[key] = from;
            }
        }
    }

    for (auto [key, from] : bestFrom) {
        int target = key.second;
        g[from].push_back(Edge(target, 1));
        reduceE++;
    }
    buildZeroTopoOrder();
    buildCSR();

    auto start = chrono::steady_clock::now();
    vector<double> BC = calcBC(true, threads);
    auto finish = chrono::steady_clock::now();
    if (!output.empty()) {
        ofstream out(output);
        out << setprecision(17);
        for (int i = 0; i < sz(BC); ++i) out << i << "," << BC[i] << '\n';
    }
    double peakRAM = Statistic::getPeakRAM();
    string nameAlgo = "TBC-optimized-vertex-and-edges";
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
    processedDynamic(filename, threads, argc > 3 ? argv[3] : "");
    return 0;
}
