#ifndef ADJACENYLUT_H
#define ADJACENYLUT_H
#include <set>
#define NMAX 600

class AdjacenyLUT
{
    public:
        bool search(int a, int b);
        void insert(int a, int b);
        AdjacenyLUT();
        virtual ~AdjacenyLUT();
    protected:
    private:
    std::set<int> adjacencySets[NMAX];
};

#endif // ADJACENYLUT_H
