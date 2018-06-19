#include "Edge.h"

Edge::Edge()
{
    //ctor
}

Edge::Edge(int a, int b, int n, int i)
{
    originNode = a;
    targetNode = b;
    edgeNumber = n;
    indiceInTargetNode = i;
}

Edge::Edge(int a, int b, int i)
{
    originNode = a;
    targetNode = b;
    indiceInTargetNode = i;
}

Edge::~Edge()
{
    //dtor
}

