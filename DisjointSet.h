//
// Created by liuzhijun on 2018/5/26.
//
#ifndef SETS_DISJOINTSET_H
#define SETS_DISJOINTSET_H

#include <stack>

struct StackedDisjointSet {
    int size;
    int * parent;
    StackedDisjointSet(int size): size(size), parent(new int [size]) {
        for (int i = 0; i < size; ++i) parent[i] = -1;
    }

    int find(int x) {
        static std::stack <int> s;
        int pt = x;
        while (parent[pt] > 0){
            s.push(pt);
            pt = parent[pt];
        }
        while (!s.empty()) {
            parent[s.top()] = pt;
            s.pop();
        }
        return pt;
    }

    void merge_root(int a, int b) {
        if (a == b) return;
        if (parent[a] > parent[b]) std::swap(a, b);
        parent[a] += parent[b];
        parent[b] = a;
    }

    void merge(int a, int b) {
        merge_root(find(a), find(b));
    }

    ~StackedDisjointSet() {
        delete [] parent;
    }
};


struct DisjointSet {
    int size;
    int * parent;
    DisjointSet(int size): size(size), parent(new int [size]) {
        for (int i = 0; i < size; ++i) parent[i] = -1;
    }

    int find(int x) {
        if (parent[x] < 0) return x;
        else return parent[x] = find(parent[x]);
    }

    void merge_root(int a, int b) {
        if (a == b) return;
        if (parent[a] > parent[b]) std::swap(a, b);
        parent[a] += parent[b];
        parent[b] = a;
    }

    void merge(int a, int b) {
        merge_root(find(a), find(b));
    }

    ~DisjointSet() {
        delete [] parent;
    }
};

#endif //SETS_DISJOINTSET_H
