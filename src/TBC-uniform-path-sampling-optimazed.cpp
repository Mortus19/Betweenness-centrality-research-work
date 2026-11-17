#include <bits/stdc++.h>
using namespace std;

#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

vector<vector<int>> temporalGraph, reverseTemporalGraph;

vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}
vector<int> originalVertex;
constexpr int MIN_TIME = 0;
constexpr int INF = 1e9;

struct BfsData {
    int depth, start;
    vector<int> dist, changedVertex, layerVertex;
    vector<double> sigma;
    BfsData(int size = 0) {
        start = -1;
        dist.assign(size, INF);
        sigma.assign(size, 0);
        depth = 0;
        changedVertex.reserve(size);
    }
    void setTarget(int t) {
        depth = 0;
        for (auto [time, id] : vertexAllTimes[t]) {
            if (time == MIN_TIME) continue;
            if (dist[id] == INF) {
                dist[id] = 0;
                sigma[id] = 1;
                layerVertex.push_back(id);
                changedVertex.push_back(id);
            }
        }
    }
    void setStart(int x) {
        dist[x] = 0;
        sigma[x] = 1;
        start = x;
        depth = 0;
        layerVertex = {x};
        changedVertex.push_back(x);
    }
    void clear() {
        depth = 0;
        start = -1;
        for (auto id : changedVertex) {
            dist[id] = INF;
            sigma[id] = 0;
        }
        layerVertex.clear();
        changedVertex.clear();
    }
};

void expandLayer(BfsData &data, vector<vector<int>> &graph, vector<int> &otherDist, int& bestDist) {
    vector<int>& dist = data.dist;
    vector<double>& sigma = data.sigma;
    data.depth++;
    vector<int> nlayer;
    for (auto x : data.layerVertex) {
        for (auto y : graph[x]) {
            if (dist[y] > dist[x] + 1) {
                dist[y] = dist[x] + 1;
                sigma[y] = sigma[x];
                data.changedVertex.push_back(y);
                nlayer.push_back(y);
                bestDist = min(bestDist, otherDist[y] + dist[y]);
            }
            else if (dist[y] == dist[x] + 1) {
                sigma[y] += sigma[x];
            }
        }
    }
    swap(nlayer, data.layerVertex);
}

BfsData source, finish;

void uniformPathSampling(int s, int t, int c, vector<double>& BC, mt19937_64& rng) {
    int sourceState = s;
    for (auto [time, id] : vertexAllTimes[s]) {
        if (time == MIN_TIME) {
            sourceState = id;
            break;
        }
    }

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

        if (source.depth <= finish.depth) {
            expandLayer(source, temporalGraph, finish.dist, bestDist);
        } else {
            expandLayer(finish, reverseTemporalGraph, source.dist, bestDist);
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
            weight.push_back(source.sigma[x] * finish.sigma[x]);
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
                BC[originalVertex[tv]] += 1;
                cand.clear();
                wcand.clear();
                for (auto y : reverseTemporalGraph[tv]) {
                    if (source.dist[y] + 1 == source.dist[tv]) {
                        cand.push_back(y);
                        wcand.push_back(source.sigma[y]);
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
                for (auto y : temporalGraph[tv]) {
                    if (finish.dist[y] + 1 == finish.dist[tv]) {
                        cand.push_back(y);
                        wcand.push_back(finish.sigma[y]);
                    }
                }
                discrete_distribution<int> tdis(wcand.begin(), wcand.end());
                tv = cand[tdis(rng)];
                BC[originalVertex[tv]] += 1;
            }
            BC[originalVertex[tv]] -= 1;
        }
    }


}

vector<double> calcBC(int k = 1) {
    int n = sz(vertexAllTimes);
    vector<double> BC(n, 0);
    mt19937_64 rng(42);
    vector<array<int, 3>> quer;
    {
        uniform_int_distribution<int> dis(0, n - 1);
        vector<pair<int, int>> q;
        for (int i = 0; i < k; i++) {
            int s = dis(rng);
            int t = dis(rng);
            while (s == t) {
                t = dis(rng);
            }
            q.push_back({s, t});
        }
        sort(all(q));
        int cnt = 1;
        for (int i = 1; i < sz(q); i++) {
            if (q[i] == q[i - 1]) {
                cnt++;
            } else {
                quer.push_back({q[i - 1].first, q[i - 1].second, cnt});
                cnt = 1;
            }
        }
        quer.push_back({q.back().first, q.back().second, cnt});
    }
    source = BfsData(sz(temporalGraph));
    finish = BfsData(sz(reverseTemporalGraph));
    for (auto [s, t, c] : quer) {
        uniformPathSampling(s, t, c, BC, rng);
    }

    for (int i = 0; i < n; i++) {
        BC[i] /= (double)k;
    }

    return BC;
}


void bfs(int v, vector<double> &BC, double coef = 1.0) {

    int n = sz(temporalGraph);
    vector<int> dist(sz(temporalGraph), INF);

    vector<int> order, marked(n, 0), blocked(n, 0);
    vector<double> delta(n, 0);
    vector<double> sigma(n, 0);
    vector<double> weights(n, 0);
    vector<vector<int>> parents(n);

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
                BC[i] += delta[num] * coef;
        }
    }
}

vector<double> calcExactBC(bool normalize = false) {
    int n = sz(vertexAllTimes);
    vector<double> BC(n, 0);
    for (int i = 0; i < n; i++) {
        bfs(i, BC);
    }
    if (normalize) {
        for (int i = 0; i < n; i++) {
            double coef = 1.0 / ((double)(n) * (n - 1));
            BC[i] *= coef;
        }
    }
    return BC;
}

double getMaxError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        res = max(res, abs(exactBC[i] - approximateBC[i]));
    }
    return res;
}

double getMaeError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        res += abs(exactBC[i] - approximateBC[i]);
    }
    res /= (double)(exactBC.size());
    return res;
}

double getRmseError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        double delta = abs(exactBC[i] - approximateBC[i]);
        res += delta * delta;
    }
    res /= (double)(exactBC.size());
    res = sqrt(res);
    return res;
}


void processedDynamic(string filename, int k) {
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
        originalVertex.push_back(i);
        N++;
    }
    for (auto [x, y, t] : edges) {
        if (distinctTimes[y].count(t) > 0) continue;
        distinctTimes[y].insert(t);
        vertexAllTimes[y].push_back({t, N});
        originalVertex.push_back(y);
        N++;
    }
    temporalGraph.resize(N);
    reverseTemporalGraph.resize(N);
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
            reverseTemporalGraph[y].push_back(i);
        }
    }

    cout << "base version" << endl;
    cout << "vertex graph :" << N << endl;
    cout << "edge graph :" << numberEdges << endl;
    auto start = chrono::steady_clock::now();
    vector<double> approximateBC = calcBC(k);
    auto finish = chrono::steady_clock::now();
    vector<double> exactBC = calcExactBC(true);
    setprecision(6);
    cout << fixed << "Time propability approximate BC in second : " << std::chrono::duration<double>(finish - start).count() << '\n';
    cout << fixed << "K : " << k << '\n';
    cout << fixed << "Max absolute error: " << getMaxError(exactBC, approximateBC) << '\n';
    cout << fixed << "MAE error: " << getMaeError(exactBC, approximateBC) << '\n';
    cout << fixed << "RMSE error: " << getRmseError(exactBC, approximateBC) << '\n';
}

int main(int argc, char* argv[]) {
    string filename = "";
    int k;
    if (argc < 3) {
        filename = "../Graphs/Temporary/1_06_college_msg.patg";
        k = 10000;
    } else {
        filename = argv[1];
        k = stoi(argv[2]);
    }
    processedDynamic(filename, k);

    return 0;
}

/*
Max absolute error: 205412.121997
MAE error: 1160.501111
RMSE error: 6662.075202

==================================
Optimized
Time propability approximate BC in second : 4.972880
K : 10000
Max absolute error: 0.003211
MAE error: 0.000112
RMSE error: 0.000266
===================================
Base
Time propability approximate BC in second : 21.149011
K : 10000
Max absolute error: 0.003749
MAE error: 0.000109
RMSE error: 0.000262

 */