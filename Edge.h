#ifndef EDGE_H
#define EDGE_H
#include <cmath>
#define PI 3.14159
#include "Position.h"

class Node;

class Edge
{
    friend class Graf;
    public:


        Edge();

        Edge(int a, int b, int n, int i);

        Edge(int a, int b, int i);

        virtual ~Edge();
    protected:
    private:
        int color;
        int edgeNumber;
        int originNode;
        int targetNode;
        int indiceInTargetNode;
};
#include "Node.h"
#endif // EDGE_H
