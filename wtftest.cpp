#include <iostream>
#include "AdjListGraph.h"
using namespace std;

int main()
{
  AdjListGraph<int,int> G(10,"0123456789");
  
  for(int i=1;i<10;i++)
    for(int j=1;j<10;j++)
     G.insert(i,j,(i*j)%11);
   G.kruskal();

  return 0;

}
