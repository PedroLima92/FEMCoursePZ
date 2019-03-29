/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyTetrahedron.h"
#include "tpanic.h"

int TopologyTetrahedron::NSideNodes(int side) {
    if (side >= nSides) 
    {
        std::cout<<"TopologyTetrahedron::NSideNodes bad parameter\n";
        DebugStop();
    }
    static int nsidenodes[15] = {1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4};
    return nsidenodes[side];
}

int TopologyTetrahedron::SideNodeIndex(int side, int node) {
    int SideNodes[6][2] =
        {
            {0, 1},         //4
            {1, 2},         //5
            {2, 0},         //6
            {0, 3},         //7
            {1, 3},         //8
            {2, 3}          //9
        };
    int FaceNodes[4][3] = 
        {
            {0, 1, 2},      //10
            {0, 1, 3},      //11
            {1, 2, 3},      //12
            {0, 2, 3}       //13
        };

    if (side < 4 && node == 0) return side;
    if (side >= 4 && side < 10 && node < 2) return SideNodes[side - 4][node];
    if (side >= 10 && side < 14 && node < 3) return FaceNodes[side - 10][node];
    if (side == 14 && node < 4) return node;
    std::cout << "TPZTetrahedron::SideNodeIndex bad parameter\n";
    DebugStop();
    return -1;
}

ElementType TopologyTetrahedron::Type() {
    return ETetraedro;
}
