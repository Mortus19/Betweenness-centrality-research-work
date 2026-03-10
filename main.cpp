#include <bits/stdc++.h>
using namespace std;
using ll = long long;
#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

vector<vector<int>> g;
vector<double> betweeness;

void bfs(int s) {
    vector<int> dist(sz(g), -1);
    vector<vector<int>> parents(sz(g));
    vector<double> num_paths(sz(g), 0);
    vector<double> delta(sz(g), 0);
    num_paths[s] = 1;
    dist[s] = 0;
    vector<int> stack;
    queue<int> q;
    q.push(s);
    cerr << s << '\n';
    while (!q.empty()) {
        auto x = q.front();
        q.pop();
        stack.push_back(x);
        for (auto &y : g[x]) {
            if (dist[y] == -1) {
                dist[y] = dist[x] + 1;
                q.push(y);
            }
            if (dist[y] == dist[x] + 1) {
                num_paths[y] += num_paths[x];
                parents[y].push_back(x);
            }
        }
    }
    for (int i = sz(stack) - 1; i >= 0; i--) {
        int x = stack[i];
        for (auto &y : parents[x]) {
            delta[y] += (num_paths[y] / num_paths[x]) * (1 + delta[x]);
        }
        if (x != s) {
            betweeness[x] += delta[x];
        }
    }
}

void solve() {
    int n, m;
    cin >> n >> m;
    g.resize(n);
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        // u--; v--;
        g[u].push_back(v);
        g[v].push_back(u);
    }
    betweeness.assign(n, 0);
    for (int i = 0; i < n; i++) {
        bfs(i);
    }
    for (int i = 0; i < n; i++) {
        betweeness[i] /= 2.0l;
    }
    for (int i = 0; i < n; i++) {
        cout << betweeness[i] << '\n';
    }
    cout << '\n';
}

signed main() {
    cerr << "HERE!\n";
    // auto can_in = freopen("input.in", "r", stdin);
    // auto can_out = freopen("output.out", "w", stdout);
    // if (can_in == NULL || can_out == NULL) {
    //     cerr << "Error with files";
    //     return 0;
    // }
    // solve();
    return 0;
}