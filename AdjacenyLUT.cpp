#include "AdjacenyLUT.h"

AdjacenyLUT::AdjacenyLUT()
{
    //ctor
}

AdjacenyLUT::~AdjacenyLUT()
{
    //dtor
}

bool AdjacenyLUT::search(int a, int b)
{
    if(a > b)
        return search(b, a);
    else
    {
        if(adjacencySets[a].count(b) == 0)
            return false;
        else
            return true;
    }
}

void AdjacenyLUT::insert(int a, int b)
{
    if(a > b)
    {
        insert(b, a);
    }
    else
    {
        adjacencySets[a].insert(b);
    }
}
