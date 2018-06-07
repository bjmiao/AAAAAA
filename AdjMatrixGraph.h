//
// Created by liuzhijun on 2018/6/7.
//

#ifndef GRAPHANDSET_GRAPH_H
#define GRAPHANDSET_GRAPH_H

#include <stdexcept>
#include <iostream>
#include <queue>
#include "DisjointSet.h"
using namespace std;

template <class vertex, class edge>
class AdjMatrixGraph {
private:
    edge ** edges;
    edge ** distance;
    int ** prev;
    vertex * vertices;
    edge noedge;
    int edgecnt, vercnt;
    int find(vertex v) const {
        for (int i = 0; i < vercnt; ++i) {
            if (vertices[i] == v) return i;
        }
        throw std::runtime_error("Such vertex does not exists");
    }
    struct eulerNode {
        int id;
        eulerNode * next;
        eulerNode(int n, eulerNode * next = nullptr): id(n), next(next) {}
    };
    struct kruskalEdge {
        int beg, end;
        edge weight;
        bool operator < (const kruskalEdge & rp) const {
            return weight > rp.weight;
        }
    };
public:
    AdjMatrixGraph(int size, const vertex d[], edge noedge);
    void insert(vertex x, vertex y, edge w);
    void insert_undirected(vertex x, vertex y, edge w) {
        insert(x, y, w);
        insert(y, x, w);
    }
    void remove(vertex x, vertex y);
    bool exist(vertex x, vertex y) const;
    ~AdjMatrixGraph();
    void printGraph() const;
    int numOfVer() const { return vercnt; }
    int numOfEdge() const { return edgecnt; }
    void dfs() const;
    void dfs(int start, bool visited []) const;
    void bfs() const;
    int getOutDegree(int v) const;
    /**
     * @return This function returns whether the graph has a ring.
     */
    bool topoSort() const {
        int cnt = 0;
        queue<int> q;
        int current;
        int * degree = new int [vercnt] {0};
        for (int i = 0; i < vercnt; ++i) {
            for (int j = 0; j < vercnt; ++j) {
                if (edges[j][i] != noedge && j != i) ++degree[i];
            }
        }
        for (int j = 0; j < vercnt; ++j)
            if (degree[j] == 0) q.push(j);

        cout << "TopoSort : " << endl;
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            ++cnt;
            cout << vertices[cur] << " ";
            for (int i = 0; i < vercnt; ++i) {
                if (edges[cur][i] != noedge && i != cur) {
                    --degree[i];
                    if (degree[i] == 0) {
                        q.push(i);
                    }
                }
            }
        }
        cout << endl;
        return cnt == vercnt;
    }
    bool hasNext(int v, bool ** been) {
        for (int i = 0; i < vercnt; ++i)
            if(edges[v][i] && !been[v][i] && i != v) return true;
        return false;
    }
    void eulerCircuit(vertex start) {

        for (int k = 0; k < vercnt; ++k) {
            int deg = getOutDegree(k);
            if (deg < 2 || deg % 2) {
                cout << "No Euler Circuit" << vertices [k]<< endl;
                return;
            }
        }

        bool ** been = new bool * [vercnt];
        for (int i = 0; i < vercnt; ++i)
            been[i] = new bool [vercnt] {false};

        int i = find(start);
        eulerNode * beg, * end;
        eulerWalk(i, beg, end, been);
        while (true) {
            eulerNode * pt = beg;
            while (pt -> next) {
                int ptn = pt->next->id;
                if (hasNext(ptn, been)) break;
                pt = pt->next;
            }
            if (pt -> next == nullptr) break;
            eulerNode * q = pt->next;
            eulerNode * a, * b;
            eulerWalk(q->id, a, b, been);
            pt->next = a;
            b->next = q->next;
            delete q;
        }


        // remove been array
        for (int j = 0; j < vercnt; ++j)
            delete [] been[j];
        delete [] been;

        // print the circuit and release space:
        cout << "Euler Circuit : " << endl;
        while(beg != nullptr) {
            cout << vertices[beg->id] << " ";
            eulerNode * ptmp = beg;
            beg = beg->next;
            delete ptmp;
        }
        cout << endl;
    }
    void eulerWalk(int start, eulerNode * & b, eulerNode * & e, bool ** been) {
        int cur = start;
        b = e = new eulerNode(cur);
        while (true) {
            int next = -1;
            for (int i = 0; i < vercnt; ++i)
                if (edges[cur][i] != noedge && !been[cur][i]) {
                    next = i;
                    been[cur][i] = been[i][cur] = true;
                    break;
                }
            if (next == -1) break;
            e -> next = new eulerNode(next);
            e = e->next;
            cur = next;
        }

    }
    void kruskal() const {
        int edgesAccepted = 0;
        DisjointSet disjointSet(vercnt);
        priority_queue<kruskalEdge> q;
        for (int i = 0; i < vercnt; ++i)
            for (int j = 0; j < vercnt; ++j) {
                if (j > i && edges[i][j] != noedge)
                    q.push(kruskalEdge{i, j, edges[i][j]});
            }

        while(edgesAccepted < vercnt - 1) {
            const kruskalEdge & e = q.top();
            int u = disjointSet.find(e.beg);
            int v = disjointSet.find(e.end);
            if (u != v) {
                edgesAccepted++;
                disjointSet.merge_root(u, v);
                cout << vertices[e.beg] << " - " << vertices[e.end]<< endl;
            }
            q.pop();
        }
    }
    void floyd() {
        // this algorithm demands that noedge is a very large number
        // the book tries to be generic but failed to do so.
        // we assume that distance and prev matrices are nonexistent before calling the function;
        distance = new edge * [vercnt];
        prev = new int * [vercnt];
        for (int i = 0; i < vercnt; ++i) {
            distance[i] = new edge [vercnt];
            prev[i] = new int [vercnt];
            for (int j = 0; j < vercnt; ++j) {
                distance[i][j] = edges[i][j];
                prev[i][j] = (edges[i][j] == noedge) ? i : -1;
            }
        }
        // prove this !
        for (int k = 0; k < vercnt; ++k) {
            for (int i = 0; i < vercnt; ++i) {
                for (int j = 0; j < vercnt; ++j) {
                    if (distance[i][k] + distance[k][j] < distance[i][j]) {
                        distance[i][j] = distance[i][k] + distance[k][j];
                        prev[i][j] = prev[k][j];
                    }
                }
            }
        }
        for (int l = 0; l < vercnt; ++l) {
            for (int i = 0; i < vercnt; ++i) {
                cout << distance[l][i] << "\t";
            }
            cout << endl;
        }
        // the matrix is not destroied so that it can be used later.
    }

    void floydPathPrint(int beg, int end) {
        if ()
    }
};

template<class vertex, class edge>
AdjMatrixGraph<vertex, edge>::~AdjMatrixGraph() {
    delete [] vertices;
    for (int i = 0; i < vercnt; ++i) delete [] edges[i];
    delete [] edges;
}

template<class vertex, class edge>
/**
 * insert assumes that none exists yet;
 * You can add a exist check before but it should be called by the user.
 */
void AdjMatrixGraph<vertex, edge>::insert(vertex x, vertex y, edge w) {
    edges[find(x)][find(y)] = w;
    ++edgecnt;
}

template<class vertex, class edge>
/**
 * This function does not check for nonexistent edges
 */
void AdjMatrixGraph<vertex, edge>::remove(vertex x, vertex y) {
    edges[find(x)][find(y)] = noedge;
    --edgecnt;
}

template<class vertex, class edge>
bool AdjMatrixGraph<vertex, edge>::exist(vertex x, vertex y) const {
    return edges[find(x)][find(y)] != noedge;
}

template<class vertex, class edge>
AdjMatrixGraph<vertex, edge>::AdjMatrixGraph(int size, const vertex d[], const edge noedge) {
    vercnt = size;
    edgecnt = 0;
    this -> noedge = noedge;
    vertices = new vertex[size];
    for (int i = 0; i < vercnt; ++i) {
        vertices[i] = d[i];
    }
    edges = new edge * [size];
    for (int j = 0; j < size; ++j) {
        edges[j] = new edge[size];
        for (int i = 0; i < size; ++i) {
            edges[j][i] = noedge;
        }
        edges[j][j] = 0; // do not forget to make edges[i][i] 0, it must be true;
    }
}

template<class vertex, class edge>
void AdjMatrixGraph<vertex, edge>::printGraph() const {
    cout << "=========print=graph========" << endl;
    for (int i = 0; i < vercnt; ++i) {
        cout << vertices[i] << "\t -> ";
        for (int j = 0; j < vercnt; ++j) {
            if (edges[i][j] != noedge && i != j)
                cout << vertices[j] << "(" << edges[i][j] << ") ";
        }
        cout << endl;
    }
    cout << endl;
}

template<class vertex, class edge>
void AdjMatrixGraph<vertex, edge>::dfs() const {
    bool * visited = new bool [vercnt];
    for (int i = 0; i < vercnt; ++i)
        visited[i] = false;
    cout << "DFS Sequences : " << endl;
    for (int j = 0; j < vercnt; ++j)
        if (!visited[j]) {
            dfs(j, visited);
            cout << endl;
        }
}

template<class vertex, class edge>
void AdjMatrixGraph<vertex, edge>::dfs(int start, bool *visited) const {
    visited[start] = true;
    cout << vertices[start] << " ";
    for (int i = 0; i < vercnt; ++i) {
        if (edges[start][i] != noedge && !visited[i])
            dfs(i, visited);
    }
}

template<class vertex, class edge>
void AdjMatrixGraph<vertex, edge>::bfs() const {
    bool * visited = new bool [vercnt];
    for (int i = 0; i < vercnt; ++i)
        visited[i] = false;
    queue<int> q;
    cout << "BFS Sequences : " << endl;
    for (int j = 0; j < vercnt; ++j) {
        if (!visited[j]) {
            visited[j] = true;
            q.push(j);
            while ( !q.empty() ) {
                int current = q.front();
                q.pop();
                cout << vertices[current] << " ";
                for (int i = 0; i < vercnt; ++i)
                    if (edges[current][i] != noedge && !visited[i]) {
                        q.push(i);
                        visited[i] = true;
                    }
            }
            cout << endl;
        }
    }
}

template<class vertex, class edge>
int AdjMatrixGraph<vertex, edge>::getOutDegree(int v) const {
    int cnt = 0;
    for (int i = 0; i < vercnt; ++i) {
        if (i != v && edges[v][i] != noedge) ++cnt;
    }
    return cnt;
}


#endif //GRAPHANDSET_GRAPH_H
