#include <bits/stdc++.h>
using namespace std;

#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

// Reduce simple segment time [l, r]
struct vertexInfo {
    int instance, timeL;

    vertexInfo(int _intstance = 1, int _l = 0) {
        instance = _intstance;
        timeL = _l;
    }
};

vector<vertexInfo> vertexData;
vector<vector<int>> reducedTimes;

vector<vector<int>> temporalGraph;
vector<vector<int>> reversedTemporalGraph;

vector<vector<int>> g;

vector<vector<pair<int, int>>> vertexAllTimes; //для вершины храним {время входа, индекс в графе}

constexpr int MIN_TIME = 0;
constexpr int MAX_TIME = 1e9;

constexpr int INF = 1e9;

vector<int> dist;

vector<int> order, marked, blocked;
vector<double> delta, sigma, weights;
vector<vector<int>> parents;

void bfs(int v, vector<double> &BC) {
    int n = sz(g);
    dist.assign(n, INF);
    order.clear();
    marked.assign(n, 0); blocked.assign(n, 0);
    delta.assign(n, 0); sigma.assign(n, 0); weights.assign(n, 0);
    parents.assign(n, {});

    int startv = -1;
    for (auto id : reducedTimes[v]) {
        if (vertexData[id].timeL == MIN_TIME) {
            startv = id;
        } else {
            blocked[id] = 1;
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
        double waysX = sigma[x] * vertexData[x].instance;
        for (auto y : g[x]) {
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
                w += sigma[num] * vertexData[num].instance;
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
            delta[y] += vertexData[y].instance * (sigma[y] / sigma[x]) * delta[x];
            if (marked[x]) {
                delta[y] += vertexData[x].instance * vertexData[y].instance * (sigma[y] / weights[x]);
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

vector<double> calcBC(bool normalize = false) {
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

void processedDynamic(string filename) {
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
    map<pair<int, int>, int> index;
    for (int i = 0; i < n; i++) {
        map<pair<vector<int>, vector<int>>, vector<int>> eqls;
        for (auto [t, idx] : vertexAllTimes[i]) {
            if (t == MIN_TIME) continue;
            eqls[{temporalGraph[idx], reversedTemporalGraph[idx]}].push_back(t);
        }
        index[{i, MIN_TIME}] = reduceN;
        vertexData.push_back({1, MIN_TIME});
        reducedTimes[i].push_back(reduceN);
        reduceN++;

        for (auto [gr, ts] : eqls) {
            for (auto t : ts) {
                index[{i, t}] = reduceN;
            }
            vertexData.push_back(vertexInfo(sz(ts), ts[0]));
            reducedTimes[i].push_back(reduceN);
            reduceN++;
        }
    }

    map<pair<int, int>, int> wasEdge;
    g.resize(reduceN);
    for (auto [x, y, t] : edges) {
        // проводим для всех вершин ребра
        int idy = index[{y, t}];
        for (auto idx : reducedTimes[x]) {
            if (vertexData[idx].timeL < t && wasEdge.count({idy, idx}) == 0) {
                g[idx].push_back(idy);
                wasEdge[{idy, idx}] = 1;
                reduceE++;
            }
        }
    }
    cout << "version vertex only optimize" << endl;
    cout << "reduced vertex graph :" << reduceN << endl;
    cout << "reduced edge graph :" << reduceE << endl;
    auto start = chrono::steady_clock::now();
    vector<double> BC = calcBC();
    auto finish = chrono::steady_clock::now();
    setprecision(6);
    cout << fixed << "Time in second : " << std::chrono::duration<double>(finish - start).count() << '\n';
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "No filename for launch" << endl;
        return 0;
    }
    string filename = argv[1];

    // vector<string> filenameData = {
    //     "../Graphs/Temporary/1_06_college_msg.patg",
    //     "../Graphs/Temporary/1_02_hospital_ward.patg",
    // };

    processedDynamic(filename);

    return 0;
}