/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 * 
 * Pedro Lima 2019
 */

#include "Topology1d.h"
#include "tpanic.h"

int Topology1d::NSideNodes(int side)
{
    if (side>2) {
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
    return -1;
}

// return the enumerated element type
ElementType Topology1d::Type(){
    return EOned;
}
