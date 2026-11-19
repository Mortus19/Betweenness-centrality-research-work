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
}

vector<int> buildProcessingOrder(const vector<int> &dist, const vector<int> &blocked) {
    int maxDist = 0;

    for (int v = 0; v < sz(dist); v++) {
        if (!blocked[v] && dist[v] != INF) {
            maxDist = max(maxDist, dist[v]);
        }
    }

    vector<vector<int>> buckets(maxDist + 1);

    for (int v : zeroTopoOrder) {
        if (!blocked[v] && dist[v] != INF) {
            buckets[dist[v]].push_back(v);
        }
    }

    vector<int> order;
    order.reserve(sz(dist));

    for (int d = 0; d <= maxDist; d++) {
        for (int v : buckets[d]) {
            order.push_back(v);
        }
    }

    return order;
}

void bfs(int v, vector<double> &BC) {

    int n = sz(g);

    vector<int> order, marked, blocked, dist;
    vector<double> bcLocal, propagate, delta, sigma, endSigma, weights;
    vector<vector<pair<int, int>>> parents;
    deque<int> q;

    dist.assign(n, INF);
    order.clear();
    marked.assign(n, 0); blocked.assign(n, 0);
    bcLocal.assign(n, 0); propagate.assign(n, 0);
    sigma.assign(n, 0); endSigma.assign(n, 0); weights.assign(n, 0);
    parents.assign(n, {});
    /*
    sigma[x]    — общее число способов дойти до сжатого экземпляра x для дальнейшего продолжения пути;

    endSigma[x] — число способов закончить путь именно в экземпляре x через настоящее temporal-ребро;

    weights[x]  — число оптимальных temporal-путей до исходной вершины, считается через endSigma;

    virtual edge — только передаёт возможность продолжения пути, но не создаёт новый endpoint-вклад.
     */
    int startv = -1;

    for (auto id : reducedTimes[v]) {
        if (vertexData[id].timeL == MIN_TIME) {
            startv = id;
        } else {
            blocked[id] = 1;
        }
    }

    dist[startv] = 0;
    sigma[startv] = 1;
    endSigma[startv] = 1;

    q.push_front(startv);

    while (!q.empty()) {
        int x = q.front();
        q.pop_front();

        for (auto e : g[x]) {
            int y = e.to;
            int w = e.w;

            if (blocked[y]) continue;

            if (dist[y] > dist[x] + w) {
                dist[y] = dist[x] + w;

                if (w == 0) {
                    q.push_front(y);
                } else {
                    q.push_back(y);
                }
            }
        }
    }

    order = buildProcessingOrder(dist, blocked);
    for (int i = 0; i < sz(order); i++) {
        int x = order[i];

        for (auto e : g[x]) {
            int y = e.to;
            int w = e.w;

            if (blocked[y]) continue;

            if (dist[y] == dist[x] + w) {
                parents[y].push_back({x, w});

                int mult = (w > 0 ? vertexData[y].instance : 1);

                double waysX = sigma[x] * mult;

                sigma[y] += waysX;

                if (w > 0) {
                    endSigma[y] += waysX;
                }
            }
        }
    }

    for (int i = 0; i < sz(reducedTimes); i++) {
        if (i == v) continue;

        int mindist = INF;

        for (auto id : reducedTimes[i]) {
            mindist = min(mindist, dist[id]);
        }

        if (mindist == INF) continue;

        double totalWays = 0;

        for (auto num : reducedTimes[i]) {
            if (dist[num] == mindist) {
                marked[num] = 1;
                totalWays += endSigma[num];
            }
        }

        for (auto num : reducedTimes[i]) {
            if (marked[num]) {
                weights[num] = totalWays;
            }
        }
    }

    for (int i = sz(order) - 1; i >= 0; i--) {
        int x = order[i];

        for (auto [y, w] : parents[x]) {
            int mult = (w > 0 ? vertexData[x].instance : 1);

            double contrib = 0.0;

            if (sigma[x] > 0) {
                contrib +=
                    (sigma[y] * mult / sigma[x]) * propagate[x];
            }

            if (w > 0 && marked[x]) {
                contrib +=
                    (sigma[y] * mult / weights[x]);
            }

            propagate[y] += contrib;

            if (w > 0) {
                bcLocal[y] += contrib;
            }
        }
    }

    for (int i = 0; i < sz(reducedTimes); i++) {
        if (i == v) continue;

        for (auto num : reducedTimes[i]) {
            if (dist[num] != INF) {
                BC[i] += bcLocal[num];
            }
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

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n; i++) {
            bfs(i, partialBC[tid]);
        }
    }

    vector<double> BC(n, 0.0);

    for (int tid = 0; tid < threads; tid++) {
        for (int i = 0; i < n; i++) {
            BC[i] += partialBC[tid][i];
        }
    }

    if (normalize) {
        double coef = 1.0 / (double(n) * (n - 1));

        for (int i = 0; i < n; i++) {
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

    auto start = chrono::steady_clock::now();
    vector<double> BC = calcBC(true, threads);
    auto finish = chrono::steady_clock::now();
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
    processedDynamic(filename, threads);
    return 0;
}