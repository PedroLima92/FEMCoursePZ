/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 Pedro Lima */

#include "Topology1d.h"
#include "tpanic.h"

int Topology1d::NSideNodes(int side) 
{
    if (side >= nSides) 
    {
        std::cout<<"Topology1d::NSideNodes bad parameter\n";
        DebugStop();
    }
    static int nsidenodes[3] = {1, 1, 2};
    return nsidenodes[side];
}

int Topology1d::SideNodeIndex(int side, int node) {
    if (side < 2 && node == 0) return side;
    if (side == 2 && node < 2) return node;
    std::cout << "Topology1d::SideNodeIndex inconsistent side or node (s" << side
            << " n" << node << ")" << std::endl;
    return -1;
}

ElementType Topology1d::Type() {
    return EOned;
}
