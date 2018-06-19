#include "Node.h"

Node::Node()
{
    //ctor
}

Node::Node(double x, double y)
{
    position = Position(x, y);
}

Node::~Node()
{
    //dtor
}


Position Node::getPosition()
{
    return position;
}

void Node::setPosition(int a, int b)
{
    position.x = a;
    position.y = b;
}
