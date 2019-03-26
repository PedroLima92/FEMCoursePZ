/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "TopologyTetrahedron.h"
#include "tpanic.h"

int TopologyTetrahedron::NSideNodes(int side) {
    DebugStop();
}

int TopologyTetrahedron::SideNodeIndex(int side, int node) {
    DebugStop();
}

ElementType TopologyTetrahedron::Type() {
    ETetraedro;
}
