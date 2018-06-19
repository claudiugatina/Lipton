#ifndef GEOMETRYFUNC_H
#define GEOMETRYFUNC_H
#include "Node.h"


class GeometryFunc
{
    public:
        /// @InOutCorrelation
        /// a curve to the right correlates to a negative return. A curve to the left correlates to a positive return.
        static double curve(Node &a, Node &b, Node &c);
        GeometryFunc();
        virtual ~GeometryFunc();
    protected:
    private:
};

#endif // GEOMETRYFUNC_H
