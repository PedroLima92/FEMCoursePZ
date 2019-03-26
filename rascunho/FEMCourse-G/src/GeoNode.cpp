/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoNode.h"
#include "tpanic.h"
#include <iomanip>

void GeoNode::Print(std::ostream &out) {

    out << ":: Node Info" << std::endl
        << " Coordinates: [" << xco[0] << ", " << xco[1] << ", " << xco[2]
        << "]" << std::endl;
}
