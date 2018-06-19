#ifndef GRAF_H
#define GRAF_H
#define NMax 50000
#define MMax 300000
#include "Edge.h"
#include "Node.h"
#include "GeometryFunc.h"
#include "AdjacenyLUT.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <set>
#include <assert.h>
#include <vector>
#include <cmath>
#include <queue>
#include <string>

//#define TESTMODE

using namespace std;

class Graf
{
    public:

        void generateGraph(int siz);

        void testFunc()
        {
            cout << GeometryFunc::curve(nodes[0], nodes[1], nodes[2]) << '\n';
            cout << GeometryFunc::curve(nodes[0], nodes[1], nodes[2]) << '\n';
            cout << GeometryFunc::curve(nodes[1], nodes[0], nodes[2]) << '\n';
        }

        int n, m; //nr noduri si nr muchii

        void afisare();

        void embed();

        void embed(int nodeIndice, int left, int right);

        int compareEdges(Edge &e1, Edge &e2);


        void liptonTarjanColoring();

        int maxColorGroup();

        Graf();
        virtual ~Graf();
        friend std::istream& operator >> (std::istream& f, Graf& G);

    protected:
    private:
    void addEdge(int a, int b, int i);
    void addEdge(int a, int b);
    void connectToNodesInFace(int startingNode, int startingEdgeIndice);
    int computeInsideCost(int v, int w, int commonParent);
    int partOfCostInside(int &startingNode, vector<int> &parent, vector<int> &inCycle, vector<int> &descendantCost, int &commonAncestor, int fromParentToChild);
    int calculateCycleCost(int &node1, int node1CycleNumber, int &node2, int node2CycleNumber, vector<int> &parent, vector<int> &inCycle, vector<int> &descendantCost, int &commonAncestor, int fromParentToChild);
    int picky(int &v, int &w);
    int nextEdge(int &vertex, int &edgeIndice);
    int previousEdge(int &vertex, int &edgeIndice);
    vector<int> computeDescendantCosts(int root);

    bool liptonTarjanSeparator(vector<int> &C, double &epsilon);

    std::vector<Node> nodes;
    std::vector<Edge> edges[NMax];
};

#endif // GRAF_H
