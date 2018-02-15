#include <bits/stdc++.h>
using namespace std;

int n,m;
vector<vector<int> > adj, adjInv;
vector<int> scc;
int sccID;
vector<bool> visited;

void dfs1(int curr, stack<int> &seen) {
    visited[curr] = true;
    for(int x : adj[curr]) {
        if(!visited[x]) {
            dfs1(x, seen);
        }
    }
    seen.push(curr);
}

void dfs2(int curr) {
    visited[curr] = true;
    scc[curr] = sccID;
    for(int x : adjInv[curr]) {
        if(!visited[x]) {
            dfs2(x);
        }
    }
}

void calcSCC() {
    visited.resize(n+1,false);
    stack<int> seen;
    for(int i = 1; i <= n; ++i) {
        if(!visited[i]) {
            dfs1(i, seen);
        }
    }
    visited.clear();
    visited.resize(n+1,false);
    sccID = 0;
    while(!seen.empty()) {
        while(!seen.empty() && visited[seen.top()]) seen.pop();
        if(!seen.empty()) {
            dfs2(seen.top());
            sccID++;
        }
    }
}

int main() {ios::sync_with_stdio(false);cin.tie(0);cout.tie(0);
    cin >> n >> m;
    scc.resize(n+1);
    adj.resize(n+1);
    adjInv.resize(n+1);
    int a,b;
    vector<pair<int, int> > edges;
    for(int i = 0; i < m; ++i) {
        cin >> a >> b;
        adj[a].push_back(b);
        adjInv[b].push_back(a);
        edges.push_back(make_pair(a,b));
    }
    calcSCC();
    //now sccID = number of SCC's
    for(int i = 1; i <= n; ++i) {
        cout << i << ' ' << scc[i] << '\n';
    }
    return 0;
}














