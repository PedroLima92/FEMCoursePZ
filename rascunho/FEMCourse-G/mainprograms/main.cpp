#include <iostream>
#include <fstream>
#include <math.h>
#include <TMatrix.h>

#include "GeoNode.h"
#include "Geom1d.h"

int main () {

    VecDouble coordinates(3);

    coordinates[0] = 0.0;
    coordinates[1] = 0.0;
    coordinates[2] = 0.0;
    GeoNode node1(coordinates);

    coordinates[0] = 1.0;
    coordinates[1] = 1.0;
    coordinates[2] = 1.0;
    GeoNode node2(coordinates);

    node1.Print(std::cout);
    node2.Print(std::cout);

    Geom1d line_test();

    return 0;
}

