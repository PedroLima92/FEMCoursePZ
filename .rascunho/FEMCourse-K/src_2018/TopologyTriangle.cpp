//
//  TopologyTriangle.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "TopologyTriangle.h"
#include "tpanic.h"

int TopologyTriangle::NSideNodes(int side)
{
    if (side>6) {
        std::cout << "TopologyTriangle::NSideNodes: Bad parameter side" << std::endl;
        DebugStop();
        return EXIT_FAILURE;
    }
    
    int nsidenodes[7] = {1,1,1,2,2,2,3};
    return nsidenodes[side];
}

// local node index of a node associated with a side
int TopologyTriangle::SideNodeIndex(int side, int node)
{
    if(side<3 && node==0)
        return side;
    if(side>=3 && side<6 && node<2)
        return (side+node)%3;
    if(side==6 && node <3)
        return node;
    
    std::cout << "TopologyTriangle::SideNodeIndex inconsistent side or node" << std::endl;
    return EXIT_FAILURE;
}

// return the enumerated element type
ElementType TopologyTriangle::Type(){
    return ETriangle;
}
