/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyQuad.h"
#include "tpanic.h"

int TopologyQuad::NSideNodes(int side)
{
    if (side >= nSides) 
    {
        std::cout<<"TopologyQuad::NSideNodes bad parameter\n";
        DebugStop();
    }
    static int nsidenodes[9] = {1, 1, 1, 1, 2, 2, 2, 2, 4};
    return nsidenodes[side];
}

int TopologyQuad::SideNodeIndex(int side, int node)
{
    if (side < 4 && node == 0)
        return side;
    if (side >= 4 && side < 8 && node < 2)
        return ((side + node) % 4); 
    if (side == 8)
        return node;
    std::cout << "TopologyQuad::SideNodeIndex bad parameter\n";
    DebugStop();
    return -1;
}

ElementType TopologyQuad::Type()
{
    return EQuadrilateral;
}
