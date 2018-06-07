//
// Created by liuzhijun on 2018/6/7.
//

#ifndef GRAPHANDSET_ADJLISTGRAPH_H
#define GRAPHANDSET_ADJLISTGRAPH_H

#include <iostream>
#include <stdexcept>
#include <queue>
#include <algorithm>
#include "DisjointSet.h"
using namespace std;

template <class vertex, class edge>
class AdjListGraph {
private:
    int vercnt, edgecnt{};
    struct edgeNode {
        int end;
        edge weight;
        edgeNode * next;
        edgeNode(int end, edge w, edgeNode * next = nullptr): end(end), weight(w), next(next) {}
    };
    struct verNode {
        vertex ver;
        edgeNode *head;

        verNode(edgeNode *head = nullptr) : head(head) {}
    } *verList;
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
    int find(vertex v) const {
        for (int i = 0; i < vercnt; ++i)
            if (verList[i].ver == v) return i;
        throw runtime_error("No such vertex");
    }
public:
    AdjListGraph(int size, const vertex d [], edge noedge = edge());
    void insert(vertex x, vertex y, edge w);
    void insert_undirected(vertex x, vertex y, edge w);
    int numOfVer() const { return vercnt; }
    int numOfEdge() const { return edgecnt; }
    void remove(vertex x, vertex y);
    void remove_id(int u, int v);
    bool exist(vertex x, vertex y) const;
    void printGraph() const;
    void dfs() const;
    void dfs(int start, bool visited []) const;
    void bfs() const;
    int getOutDegree(int v) const;
    verNode * clone() const;
    void eulerCircuit(vertex start);
    int * getInDegree() const;
    bool topoSort() const;
    /**
     * This function only accepts a DAG as input,
     * and the edge type must be a number(int or double);
     */
    void criticalPath() const {
        edge * early = new edge [vercnt]{0};
        edge * late = new edge [vercnt];
        int * top = new int [vercnt]{0};
        int * indegree = getInDegree();
        queue<int> q;
        for (int i = 0; i < vercnt; ++i) {
            if (indegree[i] == 0) q.push(i);
        }
        int cnt = 0;
        while (!q.empty()) {
            top[cnt] = q.front();
            q.pop();
            for(edgeNode * p = verList[top[cnt]].head; p != nullptr; p=p->next) {
                if (--indegree[p->end] == 0) q.push(p->end);
            }
            ++cnt;
        }
        if (cnt != vercnt) throw runtime_error("Graph not DAG");

        // Get the early time
        for (int j = 0; j < vercnt; ++j) {
            for(edgeNode * p = verList[top[j]].head; p != nullptr; p=p->next) {
                early[p->end] = max<int>(early[p->end], early[top[j]] + p->weight);
            }
        }

        for (int l = 0; l < vercnt; ++l) {
            late[l] = early[vercnt - 1];
        }

        // Get the late time
        for (int k = vercnt - 1; k >= 0; --k) {
            for(edgeNode * p = verList[top[k]].head; p != nullptr; p=p->next) {
                if (late[p->end] - p->weight < late[top[k]]) late[top[k]] = late[p->end] - p->weight;
            }
        }
        for (int n = 0; n < vercnt; ++n) {
            cout << early[n] << " ";
        }
        cout << endl;
        for (int n = 0; n < vercnt; ++n) {
            cout << late[n] << " ";
        }
        cout << endl;
        // Print the critical path
        for (int m = 0; m < vercnt; ++m) {
            if (early[top[m]] == late[top[m]])
                cout << "(" << verList[top[m]].ver << ", " << early[top[m]] <<")";
        }
    }
    void eulerWalk(int start, eulerNode * & b, eulerNode * & e);
    void kruskal() const;
    void prim(edge noedge) const;
    void printPath(int beg, int end, int * prev) const;
    void unweightedShortDistance(vertex start, int inf = INT32_MAX);
    void dijkstra(vertex start, edge noedge) const;
    bool bellmanFord(vertex start, edge noedge) const {
        edge * distance = new edge [vercnt];
        int * prev = new int [vercnt];
        int begin = find(start);
        for (int i = 0; i < vercnt; ++i) {
            distance[i] = noedge;
        }
        distance[begin ] = 0;
        prev[begin] = begin;

        for (int j = 1; j < vercnt; ++j) {
            for (int i = 0; i < vercnt; ++i) {
                for (edgeNode * p = verList[i].head; p!= nullptr; p=p->next) {
                    int u = i;
                    int v = p->end;
                    int w = p->weight;
                    if (distance[v] > distance[u] + w)
                    {
                        distance[v] = distance[u] + w;
                        prev[v] = u;
                    }
                }
            }
        }
        // another round to check for negative loops;
        for (int i = 0; i < vercnt; ++i) {
            for (edgeNode * p = verList[i].head; p!= nullptr; p=p->next) {
                int u = i;
                int v = p->end;
                int w = p->weight;
                if (distance[v] > distance[u] + w)
                {
                    cout << "Graph has negative loop" << endl;
                    delete [] distance;
                    delete [] prev;
                    return false;
                }
            }
        }

        // print path
        for (int k = 0; k < vercnt; ++k) {
            printPath(begin, k, prev);
            cout << " : " << distance[k] << endl;
        }
        delete [] distance;
        delete [] prev;
        return true;
    }
    bool spfa(vertex start, edge noedge) const {
        edge * distance = new edge [vercnt];
        int * prev = new int [vercnt];
        int * cnt = new int [vercnt]{0};
        bool * inqueue = new bool [vercnt]{false};
        queue<int> q;
        int begin = find(start);
        q.push(begin);
        for (int i = 0; i < vercnt; ++i) {
            distance[i] = noedge;
        }
        distance[begin] = 0;
        inqueue[begin] = true;
        prev[begin] = begin;
        q.push(begin);
        while(!q.empty()) {
            int now = q.front(); q.pop();
            inqueue[now] = false;
            if(++cnt[now] >= vercnt) {
                cout << "Negative loop found" << endl;
                delete [] distance;
                delete [] prev;
                delete [] cnt;
                delete [] inqueue;
                return false;
            }
            for (edgeNode * p = verList[now].head; p!=nullptr ; p=p->next) {
                int u = now;
                int v = p->end;
                int w = p->weight;
                if (distance[u] + w < distance[v]) {
                    distance[v] = distance[u] + w;
                    prev[v] = u;
                    if (!inqueue[v]) {
                        q.push(v);
                        inqueue[v] = true;
                    }
                }
            }

        }
        for (int k = 0; k < vercnt; ++k) {
            printPath(begin, k, prev);
            cout << " : " << distance[k] << endl;
        }
        delete [] distance;
        delete [] prev;
        delete [] cnt;
        delete [] inqueue;
        return true;
    }
    ~AdjListGraph();
};

template <class vertex, class edge>
AdjListGraph<vertex, edge>::AdjListGraph(int size, const vertex *d, edge noedge) {
    vercnt = size;
    edgecnt = 0;
    verList = new verNode [size];
    for (int i = 0; i < size; ++i) verList[i].ver = d[i];
}

template<class vertex, class edge>
AdjListGraph<vertex, edge>::~AdjListGraph() {
    for (int i = 0; i < vercnt; ++i) {
        edgeNode * tmp = verList[i].head;
        while (tmp) {
            verList[i].head = tmp -> next;
            delete tmp;
            tmp = verList[i].head;
        }
    }
    delete [] verList;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::insert(vertex x, vertex y, edge w) {
    int p = find(x);
    int q = find(y);
    verList[p].head = new edgeNode(q, w, verList[p].head);
    ++edgecnt;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::remove(vertex x, vertex y) {
    int u = find(x);
    int v = find(y);
    remove_id(u, v);
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::remove_id(int u, int v) {
    edgeNode * pt = verList[u].head;
    if (pt->end == v) {
        verList[u].head = pt -> next;
        delete pt; --edgecnt;
        return;
    }
    while (pt -> next != nullptr && pt->next->end != v) pt = pt->next;
    if (pt->next) {
        edgeNode * tmp = pt->next;
        pt->next = pt->next->next;
        delete tmp;
        --edgecnt;
    }
};

template<class vertex, class edge>
bool AdjListGraph<vertex, edge>::exist(vertex x, vertex y) const {
    int u = find(x);
    int v = find(y);
    edgeNode * pt = verList[u].head;
    while (pt != nullptr && pt->end != v) pt = pt->next;
    return pt != nullptr;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::printGraph() const {
    for (int i = 0; i < vercnt; ++i) {
        cout << verList[i].ver << "\t-> ";
        edgeNode * pt = verList[i].head;
        while (pt) {
            cout << verList[pt->end].ver << "(" << pt->weight << ") ";
            pt = pt->next;
        }
        cout << endl;
    }
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::dfs(int start, bool *visited) const {
    visited[start] = true;
    cout << verList[start].ver << " ";
    edgeNode * pt = verList[start].head;
    while (pt) {
        if (!visited[pt->end]) dfs(pt->end, visited);
        pt = pt->next;
    }
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::dfs() const {
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
void AdjListGraph<vertex, edge>::insert_undirected(vertex x, vertex y, edge w) {
    insert(x, y, w);
    insert(y, x, w);
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::eulerWalk(int start, AdjListGraph::eulerNode *&b, AdjListGraph::eulerNode *&e) {
    int cur = start;
    b = e = new eulerNode(cur);
    while(verList[cur].head) {
        int next = verList[cur].head->end;
        remove_id(cur, next);
        remove_id(next, cur);
        cur = next;
        e -> next = new eulerNode(cur);
        e = e -> next;
    }
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::bfs() const {
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
                cout << verList[current].ver << " ";
                edgeNode * pt = verList[current].head;
                while (pt) {
                    if (!visited[pt->end]) {
                        q.push(pt->end);
                        visited[pt->end] = true;
                    }
                    pt = pt-> next;
                }
            }
            cout << endl;
        }
    }
}

template<class vertex, class edge>
int AdjListGraph<vertex, edge>::getOutDegree(int v) const {
    edgeNode * pt = verList[v].head;
    int cnt = 0;
    while (pt) {
        ++cnt;
        pt = pt->next;
    }
    return cnt;
}

template<class vertex, class edge>
typename AdjListGraph<vertex, edge>::verNode *AdjListGraph<vertex, edge>::clone() const {
    verNode * ret = new verNode [vercnt];
    for (int i = 0; i < vercnt; ++i) {
        ret[i].ver = verList[i].ver;
        edgeNode * pt = verList[i].head;
        while (pt) {
            ret[i].head = new edgeNode(pt->end, pt->weight, ret[i].head);
            pt = pt->next;
        }
    }
    return ret;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::eulerCircuit(vertex start) {
    for (int i = 0; i < vercnt; ++i) {
        int deg = getOutDegree(i);
        if (deg % 2 || (deg < 2)) {
            cout << "No Euler Circuit";
            return;
        }
    }
    int i = find(start);
    verNode * tmp = clone();

    // find the circuit:
    eulerNode * beg, *end;
    eulerWalk(i, beg, end);
    while (true) {
        eulerNode * p = beg;
        while (p->next) {
            if (verList[p->next->id].head) break;
            else p=p->next;
        }
        if (p->next == nullptr) break;
        eulerNode * q = p->next;
        eulerNode * a, * b;
        eulerWalk(q->id, a, b);
        p->next = a;
        b->next = q->next;
        delete q;
    }
    // restore the graph
    delete [] verList;
    verList = tmp;

}

template<class vertex, class edge>
bool AdjListGraph<vertex, edge>::topoSort() const {
    int cnt = 0;
    queue<int> q;
    int current;
    int *degree = getInDegree();
    for (int j = 0; j < vercnt; ++j)
        if (degree[j] == 0) q.push(j);

    cout << "TopoSort : " << endl;
    while (!q.empty()) {
        int cur = q.front(); q.pop();
        ++cnt;
        cout << verList[cur].ver << " ";
        edgeNode * p = verList[cur].head;
        while (p) {
            if (--degree[p->end] == 0)
                q.push(p -> end);
            p = p->next;
        }
    }
    cout << endl;
    delete [] degree;
    return cnt == vercnt;
}

template<class vertex, class edge>
int *AdjListGraph<vertex, edge>::getInDegree() const {
    int * degree = new int [vercnt] {0};
    for (int i = 0; i < vercnt; ++i) {
        edgeNode * p = verList[i].head;
        while (p) {
            ++degree[p->end];
            p = p->next;
        }
    }
    return degree;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::kruskal() const {
    int edgesAccepted = 0;
    DisjointSet disjointSet(vercnt);
    priority_queue<kruskalEdge> q;
    for (int i = 0; i < vercnt; ++i)
        for (edgeNode * p = verList[i].head; p != nullptr; p = p->next)
            if (p->end > i) // so we won't add the same edge twice in an undirected graph.
                q.push(kruskalEdge{i, p->end, p->weight});

    while(edgesAccepted < vercnt - 1) {
        const kruskalEdge & e = q.top();
        int u = disjointSet.find(e.beg);
        int v = disjointSet.find(e.end);
        if (u != v) {
            edgesAccepted++;
            disjointSet.merge_root(u, v);
            cout << verList[e.beg].ver << " - " << verList[e.end].ver << endl;
        }
        q.pop();
    }
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::prim(edge noedge) const {
    bool * flag = new bool [vercnt];
    edge * lowCost = new edge [vercnt];
    int * startNode = new int [vercnt];
    for (int i = 0; i < vercnt; ++i) {
        flag[i] = false;
        lowCost[i] = noedge;
    }
    int start = 0;
    for (int j = 1; j < vercnt; ++j) {
        for (edgeNode * p = verList[start].head; p != nullptr ; p=p->next) {
            if (!flag[p->end] && lowCost[p->end] > p->weight) {
                lowCost[p->end] = p->weight;
                startNode[p->end] = start;
            }
        }

//            for (int k = 0; k < vercnt; ++k) {
//                cout << lowCost[k] << " ";
//            }
//            cout << endl;


        flag[start] = true;
        edge min = noedge;
        for (int i = 0; i < vercnt; ++i) {
            if (lowCost[i] < min && !flag[i]) {
                min = lowCost[i];
                start = i;
            }
        }
        cout << "(" << verList[startNode[start]].ver << ", " << verList[start].ver << " = " << lowCost[start] << ") ";
        lowCost[start] = noedge;
    }
    delete [] flag;
    delete [] startNode;
    delete [] lowCost;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::printPath(int beg, int end, int *prev) const {
    if (beg == end) {
        cout <<verList[end].ver;
        return;
    }
    printPath(beg, prev[end], prev);
    cout << " - " << verList[end].ver;
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::unweightedShortDistance(vertex start, int inf) {
    queue<int> q;
    int * distance = new int [vercnt];
    int * prev = new int [vercnt];
    for (int i = 0; i < vercnt; ++i) {
        distance[i] = inf;

    }
    int pt = find(start);
    distance[pt] = 0;
    prev[pt] = pt;
    q.push(pt);

    while (!q.empty()) {
        int now = q.front(); q.pop();
        for (edgeNode * p = verList[now].head; p != nullptr; p=p->next) {
            if (distance[p->end] > distance[now] + 1) {
                distance[p->end] = distance[now] + 1;
                prev[p->end] = now;
                q.push(p->end);
            }
        }
    }
    for (int j = 1; j < vercnt; ++j) {
        printPath(find(start), j, prev);
        cout << endl;
    }
}

template<class vertex, class edge>
void AdjListGraph<vertex, edge>::dijkstra(vertex start, edge noedge) const {
    edge * distance = new edge[vercnt];
    int * prev = new int [vercnt];
    bool * known = new bool [vercnt]{0};
    for (int i = 0; i < vercnt; ++i) distance[i] = noedge;
    int begin = find(start);
    distance[begin] = 0;
    prev[begin] = begin;
    for (int j = 1; j < vercnt; ++j) {
        edge min = noedge;
        int n = 0;
        for (int i = 0; i < vercnt; ++i) {
            if (!known[i] && distance[i] <= min) {
                n = i;
                min = distance[i];
            }
        }
        known[n] = true;
        for (edgeNode * p = verList[n].head; p!= nullptr ; p=p->next) {
            if (!known[p->end] && distance[p->end] > min + p->weight) {
                distance[p->end] = min + p->weight;
                prev[p->end] = n;
            }
        }
    }
    for (int k = 0; k < vercnt; ++k) {
        printPath(begin, k, prev);
        cout << " : " << distance[k] << endl;
    }
}


//struct edge {
//    int a, b, c;
//};
//edge p [] = {
//        {1, 2, 6},
//        {1, 3, 1},
//        {1, 4, 5},
//        {2, 5, 3},
//        {2, 3, 5},
//        {3, 4, 5},
//        {3, 5, 6},
//        {3, 6, 4},
//        {4, 6, 2},
//        {5, 6, 6}
//};
//int main() {
//    int name[] = {1,2,3,4,5,6};
//    AdjListGraph<int, int> testGraph(6, name, -1);
//    for (auto &&o : p)
//        testGraph.insert_undirected(o.a, o.b, o.c);
//    testGraph.dijkstra(1, 100000);
//}
#endif //GRAPHANDSET_ADJLISTGRAPH_H
