//
//  Topology1d.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include "Topology1d.h"
#include "tpanic.h"

int Topology1d::NSideNodes(int side)
{
    if (side>2) {
        std::cout << "TopologyQuad::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[3] = {1,1,2};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int Topology1d::SideNodeIndex(int side, int node)
{
    if(side <2 && node == 0)
        return side;
    if(side == 2 && node <2)
        return node;
    std::cout << "Topology1d::SideNodeIndex inconsistent side or node" << std::endl;
    return -1;
}

// return the enumerated element type
ElementType Topology1d::Type(){
    return EOned;
}
