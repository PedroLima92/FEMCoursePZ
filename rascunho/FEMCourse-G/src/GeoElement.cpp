/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElement.h"
#include "tpanic.h"

GeoElement::GeoElement() {
    DebugStop();

}

GeoElement::GeoElement(int materialid, GeoMesh *mesh, int index) {
    DebugStop();
}

GeoElement::GeoElement(const GeoElement &copy) {
    DebugStop();
}

GeoElement::~GeoElement() {
}

CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index) {
    DebugStop();
}

void GeoElement::Print(std::ostream &out) {
    DebugStop();
}
