//
//  TopologyQuad.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "TopologyQuad.h"
#include "tpanic.h"

int TopologyQuad::NSideNodes(int side)
{
    if (side>8) {
        std::cout << "TopologyQuad::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[9] = {1,1,1,1,2,2,2,2,4};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int TopologyQuad::SideNodeIndex(int side, int node)
{
    if(side<4 && node==0)
        return side;
    if(side>=4 && side<8 && node<2)
        return (side+node)%4;
    if(side==8 && node <4)
        return node;
    
    std::cout << "TopologyQuad::SideNodeIndex inconsistent side or node" << std::endl;
    return EXIT_FAILURE;
}

// return the enumerated element type
ElementType TopologyQuad::Type(){
    return EQuadrilateral;
}
