

#include <iostream>
#include <math.h>
#include "Topology1d.h"
#include "TopologyQuad.h"

using std::cin;
using std::cout;
using std::endl;

int main()
{
    std::cout << "\n\n";
    Topology1d reta;
    std::cout << "Linear element\n"
              << "Number of sides = " << reta.nSides << endl
              << "Number of corners = "<< reta.nCorners << endl
              << "Dimension = "<< reta.Dimension<< endl <<endl;
    
    TopologyQuad quadrilatero;
    std::cout << "Quadrilateral element\n"
              << "Number of sides = " << quadrilatero.nSides << endl
              << "Number of corners = "<< quadrilatero.nCorners << endl
              << "Dimension = "<< quadrilatero.Dimension<< endl <<endl;
    
    return 0;
}