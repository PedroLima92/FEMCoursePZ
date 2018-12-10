//
//  GeoNode.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "DataTypes.h"
#include "GeoNode.h"


void GeoNode::Print(std::ostream &out){
    out << "    Coordinates";
    for(int i=0;i<3;i++) out << "\t" << xco[i];
    out << "\n";
}
