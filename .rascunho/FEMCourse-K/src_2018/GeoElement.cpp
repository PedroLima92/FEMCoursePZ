//
//  GeoElement.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "CompElementTemplate.h"
#include "GeoMesh.h"

#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"

#include "tpanic.h"

GeoElement::GeoElement(){
    GMesh = 0;
    MaterialId = 0;
    Reference = NULL;
    Index = 0;
}

GeoElement::GeoElement(int materialid, GeoMesh *mesh, int index)
{
    GMesh = mesh;
    MaterialId = materialid;
    Index = index;
    if (index+1>=GetMesh()->NumElements()) {
        GetMesh()->SetNumElements(index+1);
    }
    GMesh->SetElement(index, this);
    this->Reference = NULL;
}

GeoElement::~GeoElement(){
    
}

GeoElement::GeoElement(const GeoElement &copy){
    GMesh = copy.GMesh;
    MaterialId = copy.MaterialId;
    Reference = copy.Reference;
    Index = copy.Index;
}

CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index){
    
    switch (this->Type()) {
            
        case EOned:
            return new CompElementTemplate<Shape1d>(index, mesh, this);
            break;
            
        case ETriangle:
            return new CompElementTemplate<ShapeTriangle>(index, mesh, this);
            break;
            
        case EQuadrilateral:
            return new CompElementTemplate<ShapeQuad>(index, mesh, this);
            break;
            
        case ETetraedro:
            return new CompElementTemplate<ShapeTetrahedron>(index, mesh, this);
            break;
            
        default:
            std::cout << "GeoElement::CreateCompEl : Incorrect GeoElement Type" << std::endl;
            DebugStop();
            break;
    }
}

void GeoElement::Print(std::ostream &out){
    
    out << "Element index    " << GetIndex() << std::endl;
    out << "Material index    " << Material() << std::endl;
    out << "Number of nodes    " << NNodes() << std::endl;
    out << "Corner nodes       " << NCornerNodes() << std::endl;
    out << "Nodes indexes          ";
    int i;
    for (i = 0;i < NNodes();i++) out << NodeIndex(i) << " ";
    out << "\nNumber of sides    " << NSides() << std::endl;
    out << std::endl;
    for (i = 0;i < NSides();i++) {
        out << "Neighbours for side   " << i << " : ";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this,i);
        if (!(neighbour.Element() != 0 && neighbour.Side() > -1)) {
            out << "No neighbour\n";
        }
        else {
            while ((neighbour == thisside) == false) {
                out << neighbour.Element()->GetIndex() <<  "/" << neighbour.Side() << ' ';
                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
    }
}
