#include <omp.h>

#include "StatisticTools.h"

vector<int> duplicate;
vector<vector<int>> reducedTimes;

vector<vector<int>> temporalGraph;
vector<vector<int>> reversedTemporalGraph;

vector<int> vToCompressed, compressedTov;

vector<vector<int>> g;
vector<vector<int>> rg;
struct GraphCSR {
    vector<int> offset, to;
    void fromAdj(const vector<vector<int>>& adj) {
        offset.resize(adj.size() + 1);
        offset[0] = 0;
        for (size_t v=0;v<adj.size();++v) offset[v+1] = offset[v] + (int)adj[v].size();
        to.resize(offset.back());
        for (size_t v=0;v<adj.size();++v) {
            std::copy(adj[v].begin(),adj[v].end(),to.begin()+offset[v]);
        }
    }
    int size() const {return (int)offset.size()-1;}
};
GraphCSR csrG, csrRG;


vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}

constexpr int MIN_TIME = 0;
constexpr int MAX_TIME = 1e9;

constexpr int INF = 1e9;

struct BfsData {
    int depth, start;
    uint64_t frontCost = 0;
    vector<int> dist, changedVertex, layerVertex, nextLayer;
    vector<double> sigma;
    BfsData(int size = 0) {
        start = -1;
        dist.assign(size, INF);
        sigma.assign(size, 0);
        depth = 0;
        changedVertex.reserve(size);
        layerVertex.reserve(size);
        nextLayer.reserve(size);
    }
    void setTarget(int t) {
        depth = 0;
        frontCost = 0;
        for (int i = 1; i < sz(reducedTimes[t]); i++) {
            int id = reducedTimes[t][i];
            dist[id] = 0;
            sigma[id] = 1;
            layerVertex.push_back(id);
            frontCost += csrRG.offset[id+1]-csrRG.offset[id];
            changedVertex.push_back(id);
        }
    }
    void setStart(int x) {
        dist[x] = 0;
        sigma[x] = 1;
        start = x;
        depth = 0;
        layerVertex.clear(); layerVertex.push_back(x);
        frontCost = csrG.offset[x+1]-csrG.offset[x];
        changedVertex.push_back(x);
    }
    void clear() {
        depth = 0;
        start = -1;
        frontCost = 0;
        for (auto id : changedVertex) {
            dist[id] = INF;
            sigma[id] = 0;
        }
        layerVertex.clear();
        nextLayer.clear();
        changedVertex.clear();
    }
};

void expandLayer(BfsData &data, const GraphCSR &graph, const vector<int> &otherDist, int& bestDist) {
    vector<int>& dist = data.dist;
    vector<double>& sigma = data.sigma;
    data.depth++;
    auto& nlayer = data.nextLayer;
    nlayer.clear();
    uint64_t newCost = 0;
    for (auto x : data.layerVertex) {
        double waysX = sigma[x] * duplicate[x];
        for(int p=graph.offset[x];p<graph.offset[x+1];++p) {
            int y=graph.to[p];
            if (dist[y] > dist[x] + 1) {
                dist[y] = dist[x] + 1;
                sigma[y] = waysX;
                data.changedVertex.push_back(y);
                nlayer.push_back(y);
                newCost += graph.offset[y+1]-graph.offset[y];
                if (otherDist[y] != INF) {
                    bestDist = min(bestDist, otherDist[y] + dist[y]);
                }
            }
            else if (dist[y] == dist[x] + 1) {
                sigma[y] += waysX;
            }
        }
    }
    swap(nlayer, data.layerVertex);
    data.frontCost = newCost;
}


void uniformPathSampling(BfsData &source, BfsData &finish, int s, int t, int c, vector<double>& BC, mt19937_64& rng) {
    int sourceState = reducedTimes[s][0];

    if (source.start != sourceState) {
        source.clear();
        source.setStart(sourceState);
    }
    finish.clear();
    finish.setTarget(t);
    int bestDist = INF;
    for (auto x : finish.layerVertex) {
        bestDist = min(bestDist, source.dist[x]);
    }
    while (true) {
        if (bestDist <= source.depth + finish.depth) {
            break;
        }
        if (source.layerVertex.empty() || finish.layerVertex.empty()) {
            // no path
            return;
        }

        if (source.frontCost <= finish.frontCost) {
            expandLayer(source, csrG, finish.dist, bestDist);
        } else {
            expandLayer(finish, csrRG, source.dist, bestDist);
        }
    }
    if (bestDist == INF) return;

    // пересекаются нужно найти равномерно рандомный путь
    vector<int> bestV;
    vector<double> weight;
    int leftDepth = max(0, bestDist - finish.depth);
    int rightDepth = bestDist - leftDepth;
    for (auto x : finish.layerVertex) {
        if (source.dist[x] == leftDepth && finish.dist[x] == rightDepth) {
            bestV.push_back(x);
            weight.push_back(source.sigma[x] * finish.sigma[x] * duplicate[x]);
        }
    }
    discrete_distribution<int> dis(weight.begin(), weight.end());
    for (int iter = 0; iter < c; iter++) {
        int v = bestV[dis(rng)];
        {
            // path s -> ... -> v
            vector<int> cand;
            vector<double> wcand;
            int tv = v;
            while (tv != sourceState) {
                BC[compressedTov[tv]] += 1;
                cand.clear();
                wcand.clear();
                for(int p=csrRG.offset[tv]; p<csrRG.offset[tv+1]; ++p) {
                    int y = csrRG.to[p];
                    if (source.dist[y] + 1 == source.dist[tv]) {
                        cand.push_back(y);
                        wcand.push_back(source.sigma[y] * duplicate[y]);
                    }
                }
                discrete_distribution<int> tdis(wcand.begin(), wcand.end());
                tv = cand[tdis(rng)];
            }
        }

        {
            // path v -> ... -> t
            vector<int> cand;
            vector<double> wcand;
            int tv = v;
            while (finish.dist[tv] > 0) {
                cand.clear();
                wcand.clear();
                for(int p=csrG.offset[tv]; p<csrG.offset[tv+1]; ++p) {
                    int y = csrG.to[p];
                    if (finish.dist[y] + 1 == finish.dist[tv]) {
                        cand.push_back(y);
                        wcand.push_back(finish.sigma[y] * duplicate[y]);
                    }
                }
                discrete_distribution<int> tdis(wcand.begin(), wcand.end());
                tv = cand[tdis(rng)];
                BC[compressedTov[tv]] += 1;
            }
            BC[compressedTov[tv]] -= 1;
        }
    }
}

vector<double> calcBC(int k, int threads, uint64_t seed) {
    int n = sz(vertexAllTimes);
    if (n < 2) throw invalid_argument("The graph needs at least two vertices");
    vector<double> BC(n, 0);
    mt19937_64 rng(seed);
    vector<array<int, 3>> quer;
    uniform_int_distribution<int> dis(0, n - 1);

    // Large K on moderate-size graphs: count pairs without materializing K pairs.
    // Pairs are emitted in the same lexicographic order as the sorting approach.
    int64_t cells = int64_t(n) * n;
    if (k >= 250'000 && cells <= 10'000'000LL) {
        vector<int> counts(size_t(cells), 0);
        for (int i = 0; i < k; i++) {
            int s = dis(rng), t = dis(rng);
            while (s == t) t = dis(rng);
            counts[size_t(s) * n + t]++;
        }
        for (int s = 0; s < n; s++) {
            for (int t = 0; t < n; t++) {
                int c = counts[size_t(s) * n + t];
                if (c) quer.push_back({s, t, c});
            }
        }
    } else {
        vector<pair<int, int>> q;
        q.reserve(k);
        for (int i = 0; i < k; i++) {
            int s = dis(rng), t = dis(rng);
            while (s == t) t = dis(rng);
            q.emplace_back(s, t);
        }
        sort(all(q));
        for (int i = 0; i < sz(q); ) {
            int j = i + 1;
            while (j < sz(q) && q[j] == q[i]) j++;
            quer.push_back({q[i].first, q[i].second, j - i});
            i = j;
        }
    }

    vector<pair<int, int>> batches;

    for (int i = 0; i < sz(quer); ) {
        int j = i + 1;

        while (j < sz(quer) && quer[j][0] == quer[i][0])
            j++;

        batches.push_back({i, j});
        i = j;
    }

    vector<vector<double>> partialBC(threads, vector<double>(n, 0.0));

#pragma omp parallel num_threads(threads)
    {
        int tid = omp_get_thread_num();

        BfsData source(csrG.size());
        BfsData finish(csrRG.size());

#pragma omp for schedule(dynamic)
        for (int b = 0; b < sz(batches); b++) {
            auto [l, r] = batches[b];

            int s = quer[l][0];

            mt19937_64 rng(seed + 0x9e3779b97f4a7c15ULL * uint64_t(s + 1));

            for (int i = l; i < r; i++) {
                auto [s, t, c] = quer[i];
                uniformPathSampling(source, finish, s, t, c,partialBC[tid], rng);
            }
        }
    }

    for (int tid = 0; tid < threads; tid++) {
        for (int i = 0; i < n; i++) {
            BC[i] += partialBC[tid][i];
        }
    }

    for (int i = 0; i < n; i++) {
        BC[i] /= double(k);
    }
    return BC;
}

void bfs(int v, vector<double> &BC) {

    int n = csrG.size();

    vector<int> dist(n, INF);
    vector<int> order;
    vector<int> marked(n, 0), blocked(n, 0);
    vector<double> delta(n, 0), sigma(n, 0), weights(n, 0);
    vector<vector<int>> parents(n);

    int startv = reducedTimes[v][0];
    dist[startv] = 0;
    delta[startv] = 0;
    sigma[startv] = 1;

    queue<int> q;
    q.push(startv);
    while (!q.empty()) {
        int x = q.front();
        order.push_back(x);
        q.pop();
        double waysX = sigma[x] * duplicate[x];
        for(int p=csrG.offset[x];p<csrG.offset[x+1];++p) {
            int y=csrG.to[p];
            if (blocked[y]) continue;
            if (dist[y] > dist[x] + 1) {
                dist[y] = dist[x] + 1;
                parents[y].clear();
                sigma[y] = waysX;
                q.push(y);
            }
            else if (dist[y] == dist[x] + 1) {
                sigma[y] += waysX;
            }
            if (dist[y] == dist[x] + 1) {
                parents[y].push_back(x);
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
        for (auto y : parents[x]) {
            delta[y] += duplicate[y] * (sigma[y] / sigma[x]) * delta[x];
            if (marked[x]) {
                delta[y] += duplicate[x] * duplicate[y] * (sigma[y] / weights[x]);
            }
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

vector<double> calcExactBC(bool normalize = false, int threads = 1) {
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
        for (int i = 0; i < n; i++) {
            double coef = 1.0 / ((double)(n) * (n - 1));
            BC[i] *= coef;
        }
    }
    return BC;
}

void processedDynamic(const string& filename, int k, int threads, uint64_t seed, bool noExact, const string& outputCSV) {
    ifstream in(filename);
    if (!in) throw runtime_error("Cannot open graph: " + filename);
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
    csrG.fromAdj(g);
    csrRG.fromAdj(rg);
    vector<vector<int>>().swap(g);
    vector<vector<int>>().swap(rg);
    vector<vector<int>>().swap(temporalGraph);
    vector<vector<int>>().swap(reversedTemporalGraph);
    vector<int>().swap(vToCompressed);
    auto start = chrono::steady_clock::now();
    vector<double> approximateBC = calcBC(k, threads, seed);
    auto finish = chrono::steady_clock::now();
    double peakRAM = Statistic::getPeakRAM();
        vector<double> exactBC;
    if (!noExact) exactBC = calcExactBC(true, threads);
    if (!outputCSV.empty()) {
        ofstream fout(outputCSV);
        if (!fout) throw runtime_error("Cannot write result file: " + outputCSV);
        fout << setprecision(17);
        for (double val : approximateBC) fout << val << '\n';
    }

    const double elapsed = chrono::duration<double>(finish - start).count();
    const string nameAlgo = "TBC-uniform-path-sampling-optimized-vertex-CSR";
    if (noExact) {
        cout << setprecision(10)
             << "Algo                : " << nameAlgo << '\n'
             << "Dataset             : " << filename << '\n'
             << "Threads             : " << threads << '\n'
             << "Samples (K)         : " << k << '\n'
             << "Execution time (s)  : " << elapsed << '\n'
             << "Peak RAM (MiB)      : " << peakRAM << '\n'
             << "BC sum normalized   : " << accumulate(approximateBC.begin(), approximateBC.end(), 0.0) << '\n';
    } else {
        Statistic::getAllStatistic(filename, nameAlgo, threads, exactBC, approximateBC, peakRAM, elapsed, k);
    }
}

// Backward-compatible wrapper: keep the original three-argument interface.
void processedDynamic(const string& filename, int k, int threads) {
    random_device rd;
    uint64_t seed = (uint64_t(rd()) << 32) ^ uint64_t(rd());
    processedDynamic(filename, k, threads, seed, false, "");
}

int main(int argc, char* argv[]) {
    string filename = "";
    int k, threads;
    if (argc < 4) {
        threads = 4;
        k = 1000000;
        filename = "../Graphs/Temporary/1_06_college_msg.patg";
    } else {
        filename = argv[1];
        k = stoi(argv[2]);
        threads = stoi(argv[3]);
    }
    processedDynamic(filename, k, threads);
    return 0;
}
