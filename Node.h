#ifndef NODE_H
#define NODE_H
#include <iostream>
#include "Position.h"
#define DMAX 50
#include "Edge.h"

using namespace std;


class Node
{
    public:
        Position getPosition();

        int getNumber();

        void setPosition(int x, int y);

        void afisare();

        Node();
        Node(double x, double y);
        virtual ~Node();
    protected:
    private:
        Position position;
};

#endif // NODE_H
