#include "GeometryFunc.h"

GeometryFunc::GeometryFunc()
{
    //ctor
}

GeometryFunc::~GeometryFunc()
{
    //dtor
}

double GeometryFunc::curve(Node &a, Node &b, Node &c)
{
    Position p1 = a.getPosition();
    Position p2 = b.getPosition();
    Position p3 = c.getPosition();
    double result = p1.x * p2.y + p2.x * p3.y + p3.x * p1.y - p2.y * p3.x - p3.y * p1.x - p1.y * p2.x;
    return result;
}
