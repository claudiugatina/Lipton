#include "Graf.h"

Graf::Graf()
{
    //ctor
}

Graf::~Graf()
{
    //dtor
}

void Graf::generateGraph(int siz)
{
    int edgeNumber = 0;
    n = siz * siz;
    m = 2 * n - 2 * siz;
    for(int i = 0; i < siz; ++i)
        for(int j = 0; j < siz; ++j)
        {
            nodes.push_back(Node(i, j));
        }
    for(int i = 0; i < siz; ++i)
        for(int j = 0; j < siz; ++j)
        {
            if(j < siz - 1)
            {
                addEdge(siz * i + j, siz * i + j + 1, edgeNumber);
                ++edgeNumber;
            }
            if(i < siz - 1)
            {
                addEdge(siz * i + j, siz * i + j + siz, edgeNumber);
                ++edgeNumber;
            }
        }
}


void Graf::embed()
{
    for(int i = 0; i < n; ++i)
        embed(i, 0, edges[i].size() - 1);
}

void Graf::embed(int nodeIndice, int left, int right)
{
    // sort edges clockwise
    if(left == right)
        return;
    int middle = (left + right) / 2;
    embed(nodeIndice, left, middle);
    embed(nodeIndice, middle + 1, right);

    Edge aux[right - left + 1];
    int k = 0;
    int j = middle + 1;

    for(int i = left; i <= middle; ++i)
    {
        while(j <= right && compareEdges(edges[nodeIndice][i], edges[nodeIndice][j]) == -1)
        {
            aux[k++] = edges[nodeIndice][j++];
        }
        aux[k++] = edges[nodeIndice][i];
    }

    while(j <= right)
        aux[k++] = edges[nodeIndice][j++];
    for(int i = left; i <= right; ++i)
    {
        edges[nodeIndice][i] = aux[i - left];
        edges[edges[nodeIndice][i].targetNode][edges[nodeIndice][i].indiceInTargetNode].indiceInTargetNode = i; //ii spunem muchiei din partea celuilalt nod la ce indice apare in nodul curent din partea nodului curent
    }
}

int Graf::compareEdges(Edge &e1, Edge &e2)
{
    //se presupune ca se compara muchii care pornesc de la acelasi nod
    Position p1 = nodes[e1.originNode].getPosition();
    Position p2 = nodes[e1.targetNode].getPosition();
    Position p3 = nodes[e2.targetNode].getPosition();
    double unghi1, unghi2;
    if(p1.x == p2.x)
        unghi1 = PI / 2.0 + (p2.y < p1.y) * PI;
    else
        unghi1 = atan((p2.y - p1.y) / (p2.x - p1.x)) + (p2.x < p1.x) * PI;
    if(p1.x == p3.x)
        unghi2 = PI / 2.0 + (p3.y < p1.y) * PI;
    else
        unghi2 = atan((p3.y - p1.y) / (p3.x - p1.x)) + (p3.x < p1.x) * PI;
    if(unghi1 == unghi2)
        return 0;
    if(unghi1 > unghi2)
        return 1;
    if(unghi1 < unghi2)
        return -1;
}

void Graf::afisare()
{
    for(int i = 0; i < n; ++i)
    {
        cout << i << '\n';
        for(int j = 0; j < edges[i].size(); ++j)
        {
            cout << edges[i][j].targetNode << ' ' << edges[i][j].color << "  ";
        }
        cout << "\n\n";
    }
}

void Graf::addEdge(int a, int b, int i)
{
    int ind1 = edges[a].size();
    int ind2 = edges[b].size();
    Edge e1(a, b, i, ind2), e2(b, a, i, ind1);
    edges[a].push_back(e1);
    edges[b].push_back(e2);
}

void Graf::addEdge(int a, int b)
{
    int ind1 = edges[a].size();
    int ind2 = edges[b].size();
    Edge e1(a, b, ind2), e2(b, a, ind1);
    edges[a].push_back(e1);
    edges[b].push_back(e2);
}

void Graf::connectToNodesInFace(int startingNode, int startingEdgeIndice)
{
    Node firstNode = nodes[startingNode];

    vector<Edge> result;
    enum direction{LEFT = 1, RIGHT = -1};
    direction firstCurve;
    int numberOfEdgesAdded = 0;
    int x = edges[startingNode][startingEdgeIndice].targetNode;
    int edgeIndice = edges[startingNode][startingEdgeIndice].indiceInTargetNode;
    Node secondNode = nodes[x];

    edgeIndice = (edgeIndice + 1) % edges[x].size();
    int futureEdgeIndice = edges[x][edgeIndice].indiceInTargetNode;
    x = edges[x][edgeIndice].targetNode;
    edgeIndice = futureEdgeIndice;
    Node thirdNode = nodes[x];

    firstCurve = (GeometryFunc::curve(firstNode, secondNode, thirdNode) > 0) ? LEFT : RIGHT;

    while(x != startingNode)
    {
        if(edges[x][(edgeIndice + 1) % edges[x].size()].targetNode != startingNode)
        {
            edges[x].insert(edges[x].begin() + edgeIndice + 1, Edge(x, startingNode, startingEdgeIndice + numberOfEdgesAdded + 1));
            ++numberOfEdgesAdded;
            result.push_back(Edge(startingNode, x, edgeIndice + 1));
            for(int edgeIterator = edgeIndice + 2; edgeIterator < edges[x].size(); ++edgeIterator) // not a proper iterator
            {
                int otherNode = edges[x][edgeIterator].targetNode;
                int edgeIndiceInOtherNode = edges[x][edgeIterator].indiceInTargetNode;
                ++edges[otherNode][edgeIndiceInOtherNode].indiceInTargetNode;
            }
        }
        edgeIndice = (edgeIndice + 1) % edges[x].size();
        futureEdgeIndice = edges[x][edgeIndice].indiceInTargetNode;
        x = edges[x][edgeIndice].targetNode;
        edgeIndice = futureEdgeIndice;
    }
    edges[startingNode].insert(edges[startingNode].begin() + startingEdgeIndice + 1, result.begin(), result.end());
    for(int edgeIterator = startingEdgeIndice + numberOfEdgesAdded; edgeIterator < edges[startingNode].size(); ++edgeIterator) // not a proper iterator
    {
        int otherNode = edges[startingNode][edgeIterator].targetNode;
        int edgeIndiceInOtherNode = edges[startingNode][edgeIterator].indiceInTargetNode;
        edges[otherNode][edgeIndiceInOtherNode].indiceInTargetNode += numberOfEdgesAdded;
    }
}

int Graf::partOfCostInside(int &startingNode, vector<int> &parent, vector<int> &inCycle, vector<int> &descendantCost, int &commonAncestor, int fromParentToChild)
{
    int inside = 0;
    int node = startingNode;
    while(node != commonAncestor)
    {
        int side[2];
        side[0] = 0;
        side[1] = 0;
        int whereToAdd = 0;
        bool parentFirstChildSecond = false;
        bool childFirstParentSecond = false;
        bool isLeaf = false;
        bool isCommonAncestor = false;

        for(int i = 0; i < edges[node].size(); ++i)
        {
            int targetNode = edges[node][i].targetNode;
            //cazuri: si-a gasit fiul din ciclu prima data, apoi parintele lui
            //si-a gasit mai intai parintele, apoi fiul
            if(inCycle[targetNode])
            {
                whereToAdd ^= 1;
                if(!parentFirstChildSecond && !childFirstParentSecond)
                {
                    if(parent[node] == targetNode)
                        parentFirstChildSecond = true;
                    else
                        childFirstParentSecond = true;
                }
            }
            else
            {
                if(node == parent[targetNode])
                {
                    side[whereToAdd] += descendantCost[targetNode];
                }
            }

        }
        if(parentFirstChildSecond)
        {
            inside += side[0 ^ fromParentToChild];
        }
        else
            inside += side[1 ^ fromParentToChild];
        node = parent[node];
    }
    return inside;
}

int Graf::calculateCycleCost(int &node1, int node1CycleNumber, int &node2, int node2CycleNumber, vector<int> &parent, vector<int> &inCycle, vector<int> &descendantCost, int &commonAncestor, int fromParentToChild)
{
    int cost = partOfCostInside(node1, parent, inCycle, descendantCost, commonAncestor, fromParentToChild);
    cost += partOfCostInside(node2, parent, inCycle, descendantCost, commonAncestor, 1 ^ fromParentToChild);
    int side[2];
    side[0] = 0;
    side[1] = 0;
    int whereToAdd = 0;
    bool yFirstwSecond = false;
    bool wFirstySecond = false;
    int node = commonAncestor;
    for(int i = 0; i < edges[node].size(); ++i)
    {
        int targetNode = edges[node][i].targetNode;
        //cazuri: si-a gasit fiul din ciclu prima data, apoi parintele lui
        //si-a gasit mai intai parintele, apoi fiul
        if(inCycle[targetNode] && node == parent[targetNode])
        {
            whereToAdd ^= 1;
            if(!yFirstwSecond && !wFirstySecond)
            {
                if(inCycle[targetNode] == node1CycleNumber)
                    yFirstwSecond = true;
                if(inCycle[targetNode] == node2CycleNumber)
                    wFirstySecond = true;
            }
        }
        else
        {
            if(node == parent[targetNode])
            {
                side[whereToAdd] += descendantCost[targetNode];
            }
        }
    }
    if(yFirstwSecond)
    {
        cost += side[1 ^ fromParentToChild];
    }
    else
        cost += side[fromParentToChild];
    return cost;
}

int Graf::picky(int &v, int &w)
{
    int y;
    int edgeIndice = 0;
    for(int i = 0; i < edges[v].size(); ++i)
    {
        if(edges[v][i].targetNode == w)
        {
            edgeIndice = i;
            break;
        }
    }
    edgeIndice = (edgeIndice + 1) % edges[v].size();
    int nextNode = edges[v][edgeIndice].targetNode;
    int nextNodeEdgeIndice = edges[v][edgeIndice].indiceInTargetNode;
    nextNodeEdgeIndice = (nextNodeEdgeIndice + 1) % edges[nextNode].size();
    if(edges[nextNode][nextNodeEdgeIndice].targetNode == w)
        y = nextNode;
    else
    {
        edgeIndice = edges[v][edgeIndice].indiceInTargetNode;
        edgeIndice = edgeIndice - 1 + (edgeIndice == 0) * edges[w].size();
        nextNode = edges[w][edgeIndice].targetNode;
        nextNodeEdgeIndice = edges[w][edgeIndice].indiceInTargetNode;
        nextNodeEdgeIndice = nextNodeEdgeIndice - 1 + (nextNodeEdgeIndice == 0) * edges[nextNode].size();
        assert(edges[nextNode][nextNodeEdgeIndice].targetNode == v);
        y = nextNode;
    }
    return y;
}

int Graf::nextEdge(int &vertex, int &edgeIndice)
{
    return (edgeIndice + 1) % edges[vertex].size();
}

int Graf::previousEdge(int &vertex, int &edgeIndice)
{
    return (edgeIndice - 1) + edges[vertex].size() * (edgeIndice == 0);
}

int Graf::maxColorGroup()
{
    vector<int> colorNumber(n + 1, 0);
    int maximum = 0;
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < edges[i].size(); ++j)
        {
            ++colorNumber[edges[i][j].color];
        }
    }
    for(int i = 0; i < colorNumber.size(); ++i)
    {
        if(colorNumber[i] > maximum)
            maximum = colorNumber[i];
    }
    return maximum / 2;
}

void Graf::liptonTarjanColoring()
{
    double epsilon = pow(double(n), -(1.0/3.0));
    vector<int> C(n, 0);
    while(liptonTarjanSeparator(C, epsilon));
    vector<int> coloring(m, -1);
    int nextColor = 0;
    for(int i = 0; i < n; ++i)
    {
        if(C[i] == 1)
        {
            for(int j = 0; j < edges[i].size(); ++j)
            {
                coloring[edges[i][j].edgeNumber] = nextColor;
                edges[i][j].color = nextColor;
                edges[edges[i][j].targetNode][edges[i][j].indiceInTargetNode].color = nextColor;
            }
        }
    }
    ++nextColor;
    vector<int> v(n, 0);
    for(int i = 0; i < n; ++i)
    {
        if(!v[i] && !C[i])
        {
            queue<int> BFS;
            BFS.push(i);
            while(!BFS.empty())
            {
                int node = BFS.front();
                v[node] = 1;
                for(int j = 0; j < edges[node].size(); ++j)
                {
                    int node2 = edges[node][j].targetNode;
                    if(!v[node2] && !C[node2])
                    {
                        coloring[edges[node][j].edgeNumber] = nextColor;
                        edges[node][j].color = nextColor;
                        edges[node2][edges[node][j].indiceInTargetNode].color = nextColor;
                        assert(edges[node2][edges[node][j].indiceInTargetNode].targetNode == node);
                        BFS.push(node2);
                    }
                }
                BFS.pop();
            }
        }
        ++nextColor;
    }
}

bool Graf::liptonTarjanSeparator(vector<int> &C, double &epsilon)
{
    ///THEROEM 4
    vector<int> v(n, 0);

    for(int i = 0; i < n; ++i)
    {
        ///mai intai cautam o componenta conexa din G \ C cu mai mult de n/epsilon elemente, apoi aplicam teorema
        if(!v[i] && !C[i])
        {
            vector<vector<int> > levels;
            vector<int> parent(n, -1);
            vector<int> currentPartition(n, 0);
            levels.push_back(vector<int>());
            v[i] = 1;
            currentPartition[i] = 1;
            levels.push_back(vector<int>(1, i));
            //vom folosi o parcurgere BFS pentru a gasi componente conexe din G \ C
            int c = 1; //in c se retine numarul elementelor dintr-o componenta conexa
            int lastLevel = levels.size() - 1;
            while(levels[lastLevel].size() != 0)
            {
                vector<int> nextLevel;
                for(int j = 0; j < levels[lastLevel].size(); ++j)
                {
                    int currentNode = levels[lastLevel][j];
                    for(int k = 0; k < edges[currentNode].size(); ++k)
                    {
                        int nextLevelNode = edges[currentNode][k].targetNode;
                        if(!v[nextLevelNode] && !C[nextLevelNode])
                        {
                            v[nextLevelNode] = 1;
                            nextLevel.push_back(nextLevelNode);
                            currentPartition[nextLevelNode] = 1;
                            parent[nextLevelNode] = currentNode;
                            ++c;
                        }
                    }
                }
                levels.push_back(nextLevel);
                ++lastLevel;
            }

            if((double) c > epsilon * n)
            {
                step5:

                vector<int> sumVector;
                sumVector.push_back(0); // folosim un vector de sume pentru a calcula mai repede cate noduri se afla intre diferite niveluri
                for(int i = 1; i < levels.size(); ++i)
                    sumVector.push_back(sumVector[i - 1] + levels[i].size());
                int l1 = 0;
                while(sumVector[l1] < c / 2)
                    ++l1;
                int k = sumVector[l1];
                int l0 = l1;
                while(!(levels[l0].size() + 2 * (l1 - l0) <= 2 * sqrt(k)))
                    --l0;
                int l2 = l1 + 1;
                while(!(levels[l2].size() + 2 * (l2 - l1 - 1) <= 2 * sqrt(n - k)))
                    ++l2;

                step6:
                /*Delete all verticez on level l2 and above.
                Construct a new vertex x to represent all
                vertices on levels 0 through l0.*/
                //nu vreau sa stric graful, deci o sa fac
                //altul care contine ce trebuie

                Graf G;
                int totalCost = sumVector[l2 - 1] - sumVector[l0] + 1;
                // daca nu avem macar trei noduri nu se poate pune problema unui separator
                if(totalCost < 3)
                {
                    for(int i = 0; i < levels[l0].size(); ++i)
                        C[levels[l0][i]] = 1;
                    for(int i = 0; i < levels[l2].size(); ++i)
                        C[levels[l2][i]] = 1;
                    return true;
                }
                int x = levels[l0][0];
                vector<int> froml0(n, 0);
                vector<int> froml2(n, 0);
                for(int i = 0; i < levels[l0].size(); ++i)
                {
                    froml0[levels[l0][i]] = 1;
                }
                for(int i = 0; i < levels[l2].size(); ++i)
                {
                    froml2[levels[l2][i]] = 1;
                }
                G.n = n;

                // copy all edges
                for(int i = l0 + 1; i < l2; ++i)
                {
                    for(int j = 0; j < levels[i].size(); ++j)
                    {
                        int currentNode = levels[i][j];
                        for(int k = 0; k < edges[currentNode].size(); ++k)
                        {
                            int targetNode = edges[currentNode][k].targetNode;
                            if(froml0[targetNode])
                            {
                                G.addEdge(currentNode, x);
                                parent[currentNode] = x;
                            }
                            else
                            {
                                G.edges[currentNode].push_back(edges[currentNode][k]);
                            }
                        }
                    }
                }
                parent[x] = -1;

                // cleanly delete edges to l2
                for(int i = 0; i < levels[l2 - 1].size(); ++i)
                {
                    vector<Edge> aux;
                    int skipped = 0;
                    int currentNode = levels[l2 - 1][i];
                    for(int j = 0; j < edges[currentNode].size(); ++j)
                    {
                        int targetNode = edges[currentNode][j].targetNode;
                        if(!froml2[targetNode])
                        {
                            aux.push_back(edges[currentNode][j]);
                            G.edges[targetNode][edges[currentNode][j].indiceInTargetNode].indiceInTargetNode -= skipped;
                        }
                        else
                        {
                            ++skipped;
                        }
                    }
                    G.edges[currentNode] = aux;
                }

                // cleanly delete duplicate edges of x
                {
                    AdjacenyLUT adj;
                    vector<Edge> aux;
                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode = G.edges[x][i].targetNode;
                        // daca am mai gasit o muchie din asta inseamna ca muchia curenta este duplicat.
                        // in acest caz o stergem din targetNode, si actualizam pentru toti vecinii lui targetNode de dupa aceasta muchie indiceInTargetNode
                        // trebuie sa stergem muchia si din nodul x
                        // facem un vector in care punem toate muchiile bune care pornesc de la x
                        // dupa ce copiem vectorul bun in muchiile lui x ramane un pas
                        // muchiile lui x stiu la ce indice se afla in nodurile adiacente, dar invers nu
                        // trebuie doar sa actualizam in vecinii lui x indicii noi.
                        if(adj.search(x, targetNode))
                        {
                            int indice = G.edges[x][i].indiceInTargetNode;
                            G.edges[targetNode].erase(G.edges[targetNode].begin() + indice);
                            for(int j = indice; j < G.edges[targetNode].size(); ++j)
                            {
                                int targetNode2 = G.edges[targetNode][j].targetNode;
                                int indice2 = G.edges[targetNode][j].indiceInTargetNode;
                                --G.edges[targetNode2][indice2].indiceInTargetNode;
                            }
                        }
                        else
                        {
                            adj.insert(x, targetNode);
                            aux.push_back(G.edges[x][i]);
                        }
                    }
                    G.edges[x] = aux;
                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode = G.edges[x][i].targetNode;
                        int indice = G.edges[x][i].indiceInTargetNode;
                        G.edges[targetNode][indice].indiceInTargetNode = i;
                    }
                }
                #ifdef TESTMODE
                {//testblock
                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode = G.edges[x][i].targetNode;
                        int indice = G.edges[x][i].indiceInTargetNode;
                        assert(G.edges[targetNode][indice].targetNode == x);
                    }

                    for(int i = l0 + 1; i < l2; ++i)
                    {
                        for(int j = 0; j < levels[i].size(); ++j)
                        {
                            int currentNode = levels[i][j];
                            for(int k = 0; k < G.edges[currentNode].size(); ++k)
                            {
                                int targetNode = G.edges[currentNode][k].targetNode;
                                int indice = G.edges[currentNode][k].indiceInTargetNode;
                                assert(G.edges[targetNode][indice].targetNode == currentNode);
                            }
                        }
                    }
                }
                #endif

                //Make all faces of the new graph into triangles by scanning the boundary of each face and adding edges when necessary
                //Though in reality we're just going to add an edge wherever we find two adjacent edges (v, a) and (v, b) where a and b aren't connected
                // making vab a triangle.
                {
                    AdjacenyLUT adj;
                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode = G.edges[x][i].targetNode;
                        adj.insert(x, targetNode);
                    }

                    for(int i = l0 + 1; i < l2; ++i)
                    {
                        for(int j = 0; j < levels[i].size(); ++j)
                        {
                            int currentNode = levels[i][j];
                            for(int k = 0; k < G.edges[currentNode].size(); ++k)
                            {
                                int targetNode = G.edges[currentNode][k].targetNode;
                                adj.insert(currentNode, targetNode);
                            }
                        }
                    }

                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode1 = G.edges[x][i].targetNode;
                        int indice2 = G.nextEdge(x, i);
                        int targetNode2 = G.edges[x][indice2].targetNode;
                        if(!adj.search(targetNode1, targetNode2))
                        {
                            int indice1 = G.edges[x][i].indiceInTargetNode;
                            indice2 = G.edges[x][indice2].indiceInTargetNode;
                            indice2 = G.nextEdge(targetNode2, indice2);
                            G.edges[targetNode1].insert(G.edges[targetNode1].begin() + indice1, Edge(targetNode1, targetNode2, indice2));
                            G.edges[targetNode2].insert(G.edges[targetNode2].begin() + indice2, Edge(targetNode2, targetNode1, indice1));
                            // and now clean up the mess (make indices right again)
                            ++indice1;
                            ++indice2;
                            while(indice1 < G.edges[targetNode1].size())
                            {
                                int aflictedNode = G.edges[targetNode1][indice1].targetNode;
                                int aflictedEdge = G.edges[targetNode1][indice1].indiceInTargetNode;
                                ++G.edges[aflictedNode][aflictedEdge].indiceInTargetNode;
                                ++indice1;
                            }
                            while(indice2 < G.edges[targetNode2].size())
                            {
                                int aflictedNode = G.edges[targetNode2][indice2].targetNode;
                                int aflictedEdge = G.edges[targetNode2][indice2].indiceInTargetNode;
                                ++G.edges[aflictedNode][aflictedEdge].indiceInTargetNode;
                                ++indice2;
                            }
                            adj.insert(targetNode1, targetNode2);
                        }
                    }

                    for(int i = l0 + 1; i < l2; ++i)
                    {
                        for(int j = 0; j < levels[i].size(); ++j)
                        {
                            int currentNode = levels[i][j];
                            for(int k = 0; k < G.edges[currentNode].size(); ++k)
                            {
                                int targetNode1 = G.edges[currentNode][k].targetNode;
                                int indice2 = G.nextEdge(currentNode, k);
                                int targetNode2 = G.edges[currentNode][indice2].targetNode;
                                if(!adj.search(targetNode1, targetNode2))
                                {
                                    int indice1 = G.edges[currentNode][k].indiceInTargetNode;
                                    indice2 = G.edges[currentNode][indice2].indiceInTargetNode;
                                    indice2 = G.nextEdge(targetNode2, indice2);
                                    G.edges[targetNode1].insert(G.edges[targetNode1].begin() + indice1, Edge(targetNode1, targetNode2, indice2));
                                    G.edges[targetNode2].insert(G.edges[targetNode2].begin() + indice2, Edge(targetNode2, targetNode1, indice1));
                                    // and now clean up the mess (make indices right again)
                                    ++indice1;
                                    ++indice2;
                                    while(indice1 < G.edges[targetNode1].size())
                                    {
                                        int aflictedNode = G.edges[targetNode1][indice1].targetNode;
                                        int aflictedEdge = G.edges[targetNode1][indice1].indiceInTargetNode;
                                        ++G.edges[aflictedNode][aflictedEdge].indiceInTargetNode;
                                        ++indice1;
                                    }
                                    while(indice2 < G.edges[targetNode2].size())
                                    {
                                        int aflictedNode = G.edges[targetNode2][indice2].targetNode;
                                        int aflictedEdge = G.edges[targetNode2][indice2].indiceInTargetNode;
                                        ++G.edges[aflictedNode][aflictedEdge].indiceInTargetNode;
                                        ++indice2;
                                    }
                                    adj.insert(targetNode1, targetNode2);
                                }
                            }
                        }
                    }
                }

                #ifdef TESTMODE
                {//testblock
                    for(int i = 0; i < G.edges[x].size(); ++i)
                    {
                        int targetNode = G.edges[x][i].targetNode;
                        int indice = G.edges[x][i].indiceInTargetNode;
                        assert(G.edges[targetNode][indice].targetNode == x);
                    }

                    for(int i = l0 + 1; i < l2; ++i)
                    {
                        for(int j = 0; j < levels[i].size(); ++j)
                        {
                            int currentNode = levels[i][j];
                            for(int k = 0; k < G.edges[currentNode].size(); ++k)
                            {
                                int targetNode = G.edges[currentNode][k].targetNode;
                                int indice = G.edges[currentNode][k].indiceInTargetNode;
                                assert(G.edges[targetNode][indice].targetNode == currentNode);
                            }
                        }
                    }
                }
                #endif

                #ifdef TESTMODE
                {

                    for(int i = l0 + 1; i < l2; ++i)
                    {
                        for(int j = 0; j < levels[i].size(); ++j)
                        {
                            int currentNode = levels[i][j];
                            AdjacenyLUT adj;
                            for(int k = 0; k < G.edges[currentNode].size(); ++k)
                            {
                                int targetNode = G.edges[currentNode][k].targetNode;
                                if(adj.search(currentNode, targetNode))
                                    cout << "Duplicate: " << currentNode << ' ' << targetNode << '\n';
                                assert(!adj.search(currentNode, targetNode));
                                adj.insert(currentNode, targetNode);
                            }
                        }
                    }
                }
                #endif

                //choose any non-tree edge (v, w)
                Edge nonTreeEdge;
                int v, w;
                for(int i = l0 + 1; i < l2; ++i)
                {
                    for(int j = 0; j < levels[i].size(); ++j)
                    {
                        for(int k = 0; k < G.edges[levels[i][j]].size(); ++k)
                        {
                            nonTreeEdge = G.edges[levels[i][j]][k];
                            v = nonTreeEdge.originNode;
                            w = nonTreeEdge.targetNode;
                            if(!(v == parent[w] || w == parent[v]))
                                goto restOfStep8;
                        }
                    }
                }

                restOfStep8:
                int numberOfNodesInCycle = 0;
                vector<int> inCycle(n, 0);
                inCycle[v] = 1;
                int node = v;
                int commonAncestor;
                while(node != -1)
                {
                    inCycle[node] = 1;
                    node = parent[node];
                    ++numberOfNodesInCycle;
                }
                node = w;
                while(!inCycle[node])
                {
                    inCycle[node] = 2;
                    node = parent[node];
                    ++numberOfNodesInCycle;
                }
                commonAncestor = node;
                while(parent[node] != -1)
                {
                    node = parent[node];
                    inCycle[node] = 0;
                    --numberOfNodesInCycle;
                }
                vector<int> descendantCost(n, 1);
                //compute descending costs
                for(int i = l2 - 1; i > l0; --i)
                {
                    for(int j = 0; j < levels[i].size(); ++j)
                    {
                        descendantCost[parent[levels[i][j]]] += descendantCost[levels[i][j]];
                    }
                }

                //compute cost inside and outside of cycle
                //parcurgem separat muchiile lui commonAncestor
                int inside = 0;
                int insideIsOutside = 0;
                int side[2];
                side[0] = 0;
                side[1] = 0;
                int whereToAdd = 0;
                int order = 0;
                bool vFirstwSecond = false;
                bool wFirstvSecond = false;
                for(int i = 0; i < G.edges[commonAncestor].size(); ++i)
                {
                    int targetNode = G.edges[commonAncestor][i].targetNode;
                    // inside is outside daca gaseste parintele intre ramura 1 si ramura 2
                    // ramurile acestea fiind cele care duc la v si w
                    if(inCycle[targetNode])
                    {
                        whereToAdd ^= 1;
                        order = order * 10 + inCycle[targetNode]; // 1 - v, 2 - w, 3 - parent

                    }
                    else
                    {
                        if(commonAncestor == parent[targetNode])
                        {
                            side[whereToAdd] += descendantCost[targetNode];
                        }
                        else
                            if(targetNode == parent[commonAncestor])
                            {
                                order = order * 10 + 3;
                                side[whereToAdd] += totalCost - descendantCost[commonAncestor];
                            }
                    }
                }
                switch(order)
                {
                    case 132:
                        insideIsOutside = 1;
                        break;
                    case 213:
                        insideIsOutside = 1;
                        break;
                    case 321:
                        insideIsOutside = 1;
                        break;
                    default:
                        break;
                }
                if((order % 10 == 1 || (order / 10) % 10 == 1) && ((order / 10) % 10 == 2 || order / 100 == 2))
                {
                    inside += side[1 ^ insideIsOutside];
                }
                else
                    inside += side[insideIsOutside];

                inside += partOfCostInside(v, parent, inCycle, descendantCost, commonAncestor, insideIsOutside);
                inside += partOfCostInside(w, parent, inCycle, descendantCost, commonAncestor, 1 ^ insideIsOutside);

                int outside = totalCost - inside - numberOfNodesInCycle;
                if(outside > inside)
                {
                    inside = outside;
                    insideIsOutside ^= 1;
                }

                step9:
                if(inside <= 2 * totalCost / 3)
                {
                    //zicem ca toate nodurile din l0, l2 si din ciclul gasit fac parte din C (separatorul final)
                    while(v != commonAncestor)
                    {
                        C[v] = 1;
                        v = parent[v];
                    }
                    while(w != commonAncestor)
                    {
                        C[w] = 1;
                        w = parent[w];
                    }
                    C[commonAncestor] = 1;
                    for(int i = 0; i < levels[l0].size(); ++i)
                        C[levels[l0][i]] = 1;
                    for(int i = 0; i < levels[l2].size(); ++i)
                        C[levels[l2][i]] = 1;
                    return true;
                }

                if(!insideIsOutside)
                {
                    int y = G.picky(v, w);
                    // ne bazam ca y nu e parintele amandurora
                    if(y == parent[w])
                    {
                        inCycle[w] = 0;
                        w = y;
                        goto step9;
                    }
                    if(y == parent[v])
                    {
                        inCycle[v] = 0;
                        v = y;
                        goto step9;
                    }

                    // ramura 'if neither .. is a tree edge':
                    int z = y;
                    int costFromYToZ = 0;
                    while(!inCycle[z])
                    {
                        inCycle[z] = 3;
                        z = parent[z];
                        ++costFromYToZ;
                    }

                    int costInsideYW = partOfCostInside(y, parent, inCycle, descendantCost, z, 0);
                    costInsideYW += partOfCostInside(w, parent, inCycle, descendantCost, z, 1);
                    {
                        // there should be a function for this. Though it seems a lot of times a block like this is used, there are slight differences
                        // in relevant nodes and information gathered
                        int side[2];
                        side[0] = 0;
                        side[1] = 0;
                        int whereToAdd = 0;
                        bool yFirstwSecond = false;
                        bool wFirstySecond = false;
                        node = z;
                        for(int i = 0; i < G.edges[node].size(); ++i)
                        {
                            int targetNode = G.edges[node][i].targetNode;
                            //cazuri: si-a gasit fiul din ciclu prima data, apoi parintele lui
                            //si-a gasit mai intai parintele, apoi fiul
                            if(inCycle[targetNode] && node == parent[targetNode])
                            {
                                whereToAdd ^= 1;
                                if(!yFirstwSecond && !wFirstySecond)
                                {
                                    if(inCycle[targetNode] == 3)
                                        yFirstwSecond = true;
                                    if(inCycle[targetNode] == 2)
                                        wFirstySecond = true;
                                }
                            }
                            else
                            {
                                if(node == parent[targetNode])
                                {
                                    side[whereToAdd] += descendantCost[targetNode];
                                }
                            }
                        }
                        if(yFirstwSecond)
                        {
                            costInsideYW += side[1];
                        }
                        else
                            costInsideYW += side[0];
                    }

                    int costInsideYV = inside - costFromYToZ - costInsideYW;
                    if(costInsideYW > costInsideYV)
                    {
                        if(inCycle[z] == 1)
                        {
                            while(v != z)
                            {
                                inCycle[v] = 0;
                                v = parent[v];
                            }
                            v = y;
                            int cleaner = y;
                            while(inCycle[cleaner] == 3)
                            {
                                inCycle[cleaner] == 1;
                                cleaner = parent[cleaner];
                            }
                        }
                        else
                        {
                            commonAncestor = z;
                            int cleaner = parent[commonAncestor];
                            while(inCycle[cleaner])
                            {
                                inCycle[cleaner] = 0;
                                cleaner = parent[cleaner];
                            }
                            cleaner = v;
                            while(inCycle[cleaner])
                            {
                                inCycle[cleaner] = 0;
                                cleaner = parent[cleaner];
                            }
                            v = y;
                            cleaner = y;
                            while(inCycle[cleaner] == 3)
                            {
                                inCycle[cleaner] == 1;
                                cleaner = parent[cleaner];
                            }
                        }
                    }
                    else
                    {
                        if(inCycle[z] == 1)
                        {
                            commonAncestor = z;
                            int cleaner = parent[commonAncestor];
                            while(cleaner != -1 && inCycle[cleaner])
                            {
                                inCycle[cleaner] = 0;
                                cleaner = parent[cleaner];
                            }
                            cleaner = y;
                            while(cleaner != commonAncestor)
                            {
                                inCycle[cleaner] = 2;
                                cleaner = parent[cleaner];
                            }
                            cleaner = w;
                            while(cleaner != -1 && inCycle[cleaner])
                            {
                                inCycle[cleaner] = 0;
                                cleaner = parent[cleaner];
                            }
                            w = y;
                            cleaner = w;
                            while(cleaner != -1 && inCycle[cleaner] == 3)
                            {
                                inCycle[cleaner] == 2;
                                cleaner = parent[cleaner];
                            }
                        }
                        else
                        {
                            int cleaner = w;
                            while(cleaner != z)
                            {
                                inCycle[cleaner] = 0;
                                cleaner = parent[cleaner];
                            }
                            w = y;
                            cleaner = w;
                            while(inCycle[cleaner] == 3)
                            {
                                inCycle[cleaner] == 2;
                                cleaner = parent[cleaner];
                            }
                        }
                    }
                    goto step9;
                }
                else
                {
                    int y;
                    int edgeIndice;
                    for(int i = 0; i < G.edges[v].size(); ++i)
                    {
                        if(G.edges[v][i].targetNode == w)
                        {
                            edgeIndice = i;
                            break;
                        }
                    }
                    edgeIndice = G.previousEdge(v, edgeIndice);
                    int nextNode = G.edges[v][edgeIndice].targetNode;
                    int nextNodeEdgeIndice = G.edges[v][edgeIndice].indiceInTargetNode;
                    nextNodeEdgeIndice = G.previousEdge(nextNode, nextNodeEdgeIndice);
                    if(G.edges[nextNode][nextNodeEdgeIndice].targetNode == w)
                        y = nextNode;
                    else
                    {
                        // daca intra pe ramura asta inseamna ca ceva nu a mers bine
                        edgeIndice = G.nextEdge(v, edgeIndice);
                        edgeIndice = G.edges[v][edgeIndice].indiceInTargetNode;
                        edgeIndice = G.nextEdge(w, edgeIndice);
                        nextNode = G.edges[w][edgeIndice].targetNode;
                        nextNodeEdgeIndice = G.edges[w][edgeIndice].indiceInTargetNode;
                        nextNodeEdgeIndice = G.previousEdge(nextNode, nextNodeEdgeIndice);
                        assert(G.edges[nextNode][nextNodeEdgeIndice].targetNode == v);
                        y = nextNode;
                    }

                    if(v == parent[y])
                    {
                        inside -= 1;
                        v = y;
                        goto step9;
                    }
                    if(w == parent[y])
                    {
                        inside -= 1;
                        w = y;
                        goto step9;
                    }

                    // ramura 'if neither .. is a tree edge'
                    // gasim un nou 'common ancestor' urmarind parintii de la y si de la common ancestor.
                    // ciclul nou va fi format din y si nodul din {v, w} pentru care nodul din {v, w} care nu a fost ales se afla in interiorul ciclului
                    // interior insemnand interior pe bune, nu ca ii zicem noi interior la exterior doar pentru ca e mai mare

                    int z = y;
                    while(!inCycle[z] && z != -1)
                    {
                        inCycle[z] = 3;
                        z = parent[z];
                    }

                    if(z != commonAncestor)
                    {
                        // cazul in care z a ajuns pana la radacina fara sa dea de ciclu
                        if(z == -1)
                        {
                            int commonAncestorLastValue = commonAncestor;
                            while(!inCycle[commonAncestor] && parent[commonAncestor] != -1)
                            {
                                inCycle[commonAncestor] = 1;
                                commonAncestor = parent[commonAncestor];

                            }
                            // curatam vectorul inCycle, ca sa nu contina 1 si pentru noduri care de fapt nu fac parte din ciclu
                            z = commonAncestor;
                            while(z != -1 && inCycle[z])
                            {
                                z = parent[z];
                                inCycle[z] = 0;
                            }
                            // parcurgem muchiile lui commonAncestor pentru a vedea in ce ordine apar muchiile relevante
                            // cazuri:
                            // v si w apar intre parinte si fiul din ciclu => il luam pe cel care contine ciclul precedent
                            // v si w sunt separati de parinte si de fiu => nu putem lua nici un ciclu care sa il contina pe cel curent
                            // deci vom calcula costurile pentru (y, v) si (y, w), si il vom alege pe cel cu costul mai mare.
                            int order[4];
                            int indice = 0;
                            bool oneAfterAnother = false;
                            for(int i = 0; i < G.edges[commonAncestor].size(); ++i)
                            {
                                int targetNode = G.edges[commonAncestor][i].targetNode;
                                if(targetNode == v)
                                    order[indice++] = 1;
                                if(targetNode == w)
                                    order[indice++] = 2;
                                if(targetNode == parent[commonAncestor])
                                    order[indice++] = 3;
                                if(parent[targetNode] == commonAncestor && inCycle[targetNode])
                                    order[indice++] = 4;
                            }
                            // luam toate cazurile in care v si w sunt una dupa alta
                            // cazuri posibile:
//                            1234 1243 1324 1342 1423 1432
//                            2134 2143 2314 2341 2413 2431
//                            3124 3142 3214 3241 3412 3421
//                            4123 4132 4213 4231 4312 4321
                            // conditia ca v si w sa se afle pe pozitiile i si j este order[i] + order[j] == 3
                            // verificam daca se afla pe pozitii consecutive
                            for(indice = 0; indice < 4; ++i)
                                if(order[indice] + order[(indice + 1) % 4] == 3)
                                {
                                    oneAfterAnother = true;
                                    break;
                                }
                            if(oneAfterAnother)
                            {
                                if(order[(indice - 1 + 4 * (indice == 0))] == 3)
                                {
                                    if(order[indice] == v)
                                    {
                                        // luam ciclul (v, y)
                                        // recalculam costul ciclului
                                        inside = partOfCostInside(v, parent, inCycle, descendantCost, parent[commonAncestor], 0);
                                        w = y;
                                        goto step9;
                                    }
                                    else
                                    {
                                        // luam ciclul (w, y)
                                        inside = partOfCostInside(w, parent, inCycle, descendantCost, parent[commonAncestor], 1);
                                        v = y;
                                        goto step9;
                                    }
                                }
                                else
                                {
                                    if(order[indice] == v)
                                    {
                                        // luam ciclul (w, y)
                                        inside = partOfCostInside(w, parent, inCycle, descendantCost, parent[commonAncestor], 1);
                                        v = y;
                                        goto step9;
                                    }
                                    else
                                    {
                                        inside = partOfCostInside(v, parent, inCycle, descendantCost, parent[commonAncestor], 0);
                                        w = y;
                                        goto step9;
                                    }
                                }
                            }
                            else
                            {
                                // cazul in care v si w sunt de o parte si de alta
                                int inside2 = partOfCostInside(w, parent, inCycle, descendantCost, parent[commonAncestor], 1);
                                inside = partOfCostInside(v, parent, inCycle, descendantCost, parent[commonAncestor], 0);
                                if(inside2 > inside)
                                {
                                    v = y;
                                    goto step9;
                                }
                                else
                                {
                                    w = y;
                                    goto step9;
                                }
                            }
                        }
                        if(inCycle[z] == 1)
                        {
                            inside -= partOfCostInside(y, parent, inCycle, descendantCost, z, 0);
                            v = y;
                            while(inCycle[y] == 3)
                            {
                                inCycle[y] = 1;
                                y = parent[y];
                            }
                            goto step9;
                        }
                        if(inCycle[z] == 2)
                        {
                            inside -= partOfCostInside(y, parent, inCycle, descendantCost, z, 1);
                            w = y;
                            while(inCycle[y] == 3)
                            {
                                inCycle[y] = 2;
                                y = parent[y];
                            }
                            goto step9;
                        }


                    }
                    else
                    {

                    }
                }
            }
            else
            {
                levels.erase(levels.begin(), levels.begin() + levels.size());
                c = 0;
            }
        }
    }
    return false;
}


std::istream& operator >> (std::istream& f, Graf& G)
{
    f >> G.n >> G.m;
    for(int i = 0; i < G.n; ++i)
    {
        double x, y;
        f >> x >> y;
        G.nodes.push_back(Node(x, y));
    }
    for(int i = 0; i < G.m; ++i)
    {
        int a, b;

        f >> a >> b;
        if(a >= G.n || b >= G.n || b < 0 || a < 0)
            cerr << "Input format is incorrect";
        G.addEdge(a, b, i);
    }
}

