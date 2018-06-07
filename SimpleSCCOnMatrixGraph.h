//
// Created by liuzhijun on 2018/6/7.
//

#ifndef GRAPHANDSET_SIMPLESCCONMATRIXGRAPH_H
#define GRAPHANDSET_SIMPLESCCONMATRIXGRAPH_H
#include <iostream>
#include <stack>
using namespace std;
bool map[100][100];

stack<int> s;
int n, m;
bool been[100];

void dfs(int pos, bool reverse=false, bool print=false) {
    for (int i = 0; i < n; ++i) {
        if (!been[i] && i != pos && (reverse ? map[i][pos] : map[pos][i])) {
            been[i] = 1;
            dfs(i, reverse, print);
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
        map[a][b] = true;
    }
    // first round
    for (int j = 0; j < n; ++j) {
        if (!been[j]) {
            been[j] = 1;
            dfs(j);
        }
    }

    // clear memory
    memset(been, 0, sizeof(been));

    // second round
    while(!s.empty()) {
        int top = s.top(); s.pop();
        if (!been[top]) {
            been[top] = true;
            dfs(top, true, true);
        }
        cout << endl;
    }
}

#endif //GRAPHANDSET_SIMPLESCCONMATRIXGRAPH_H
