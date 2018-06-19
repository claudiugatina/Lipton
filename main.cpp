#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#define NMax 100000
#define MMax 600000
#define PI 3.1415
#include "Graf.h"
using namespace std;

int main()
{
//    freopen("test3.in", "r", stdin);
    freopen("liptontarjan.out", "w", stdout);
//
//    Graf G;
//    cin >> G;
//  //  G.testFunc();
//    G.embed();
//    G.afisare();
//    cout << "\n\n\n\n";
//    G.liptonTarjanColoring();
//    G.afisare();
//    cout << '\n' << G.maxColorGroup();
//    Graf G;
//    G.generateTree(100, 5);
//    G.afisare();
//    G.liptonTarjanColoring();
//    G.afisare();
//    cout << '\n' << G.maxColorGroup();
    {
        for(int n = 31; n <= 45; ++n)
        {

                Graf G;
                G.generateGraph(n);
              //  G.afisare();
                G.embed();
                clock_t before = clock();
                G.liptonTarjanColoring();
                clock_t after = clock();
                double duration = (after - before) / (double)CLOCKS_PER_SEC;
                int maxColorGroup = G.maxColorGroup();
                cout << n << " & " << duration << " & " << maxColorGroup << " \\\\ \n";
        }
    }
    return 0;
}
