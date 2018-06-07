#include "AdjMatrixGraph.h"

struct edge {
    int a, b, c;
};
edge p [] = {
        {1, 2, 6},
        {1, 3, 1},
        {1, 4, 5},
        {2, 5, 3},
        {2, 3, 5},
        {3, 5, 6},
        {3, 6, 4},
        {4, 6, 2},
        {5, 6, 6}
};
int main() {
    int name[] = {1,2,3,4,5,6};
    AdjMatrixGraph<int, int> testGraph(6, name, 10000000);
    for (auto &&o : p)
        testGraph.insert_undirected(o.a, o.b, o.c);
    testGraph.insert(3, 4, 5);
    testGraph.insert(4, 3, -5);
    testGraph.floyd();
}