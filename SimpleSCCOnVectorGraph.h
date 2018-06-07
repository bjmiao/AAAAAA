//
// Created by liuzhijun on 2018/6/7.
//

#ifndef GRAPHANDSET_SCCONVECTORGRAPH_H
#define GRAPHANDSET_SCCONVECTORGRAPH_H

#include <iostream>
#include <vector>
#include <stack>
using namespace std;

vector<int> graph[100], inv[100];
int n, m;
stack<int> s;
bool been[100] {0};

void dfs(int pos, vector<int> graph [], bool been [], stack<int> & s, bool print = false) {
    int len = graph[pos].size();
    for (int i = 0; i < len; ++i) {
        if (! been[graph[pos][i]]) {
            been[graph[pos][i]] = true;
            dfs(graph[pos][i], graph, been, s, print);
        }
    }
    if (print) {
        cout << pos << " ";
    } else s.push(pos);
}

int main() {
    cin >> n >> m;
    for (int i = 0; i < m; ++i) {
        int a, b;
        cin >> a >> b;
        graph[a].push_back(b);
        inv[b].push_back(a);
    }

    // first round dfs
    for (int j = 1; j <= n; ++j)
        if (!been[j]) dfs(j, graph, been, s);

    // clear memory
    memset(been, 0, sizeof(been));
    // second round dfs
    while (!s.empty()) {
        int top = s.top();
        s.pop();
        if (!been[top]) {
            been[top] = true;
            dfs(top, inv, been, s, true);
            cout << endl;
        }
    }
}
/*
5 7
1 2
2 1
2 3
3 4
4 3
5 3
4 5
*/

/*
4 4
1 2
2 1
3 4
4 3

*/


#endif //GRAPHANDSET_SCCONVECTORGRAPH_H
