#include <bits/stdc++.h>
using namespace std;

#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

vector<vector<int>> temporalGraph;
vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}
vector<int> originalVertex;
constexpr int MIN_TIME = 0;
constexpr int INF = 1e9;

struct BfsData {
    vector<int> dist;
    vector<double> sigma;
    vector<vector<int>> parents;
    BfsData(int size = 0) {
        dist.assign(size, INF);
        sigma.assign(size, 0);
        parents.resize(size);
    }
};


BfsData shortesPath;

void calcShortesPath(int s) {
    int n = sz(temporalGraph);
    auto& obj = shortesPath;
    obj = BfsData(n);
    auto& dist = obj.dist;
    auto& sigma = obj.sigma;
    auto& parents = obj.parents;
    vector<int> blocked(n, 0);
    int startv = -1;
    for (auto [t1, idv] : vertexAllTimes[s]) {
        if (t1 == MIN_TIME) {
            startv = idv;
        } else {
            blocked[idv] = 1;
        }
    }
    dist[startv] = 0;
    sigma[startv] = 1;
    queue<int> q;
    q.push(startv);
    while (!q.empty()) {
        int x = q.front();
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
}


void uniformPathSampling(int s, int t, vector<double> &BC, mt19937_64& rng) {
    auto& parents = shortesPath.parents;
    auto& dist = shortesPath.dist;
    auto& sigma = shortesPath.sigma;
    int start = -1;
    for (auto [t, ids] : vertexAllTimes[s]) {
        if (t == MIN_TIME) {
            start = ids;
        }
    }
    vector<double> p;
    auto uniformParentChoose = [&](int x) {
        p.clear();
        p.reserve(sz(parents[x]));
        for (auto y : parents[x]) {
            p.push_back(sigma[y]);
        }
        discrete_distribution<int> dis(p.begin(), p.end());
        return parents[x][dis(rng)];
    };

    vector<int> starts;
    int minDist = INF;
    for (auto [t, id] : vertexAllTimes[t]) {
        minDist = min(minDist, dist[id]);
    }
    if (minDist == INF)
        return;
    for (auto [t, id] : vertexAllTimes[t]) {
        if (dist[id] == minDist) {
            starts.push_back(id);
            p.push_back(sigma[id]);
        }
    }
    discrete_distribution<int> dis(p.begin(), p.end());
    int cur = starts[dis(rng)];

    while (cur != start) {
        cur = uniformParentChoose(cur);
        if (cur != start)
            BC[originalVertex[cur]] += 1;
    }

}

vector<double> calcBC(int k = 1) {
    int countBFS = 0;
    int n = sz(vertexAllTimes);

    vector<double> BC(n, 0);

    mt19937_64 rng(42);
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
    int pred = -1;
    for (auto [s, t] : q) {
        if (s != pred) {
            calcShortesPath(s);
        }
        pred = s;
        uniformPathSampling(s, t, BC, rng);
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
 */