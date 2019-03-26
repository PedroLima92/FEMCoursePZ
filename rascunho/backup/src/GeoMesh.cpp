/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoMesh.h"
#include "tpanic.h"

GeoMesh::GeoMesh() {
    DebugStop();
}

GeoMesh::GeoMesh(const GeoMesh &copy) {
    DebugStop();
}

GeoMesh &GeoMesh::operator=(const GeoMesh &copy) {
    DebugStop();
}

void GeoMesh::SetNumNodes(int nnodes) {
    DebugStop();
}

void GeoMesh::SetNumElements(int numelements) {
    DebugStop();
}

int GeoMesh::NumNodes() {
    DebugStop();
}

int GeoMesh::NumElements() {
    DebugStop();
}

GeoNode &GeoMesh::Node(int node) {
    DebugStop();
}

void GeoMesh::SetElement(int elindex, GeoElement *gel) {
    DebugStop();
}

GeoElement *GeoMesh::Element(int elindex) {
    DebugStop();
}

void GeoMesh::BuildConnectivity() {
    DebugStop();
}

void GeoMesh::Print(std::ostream &out) {
    DebugStop();
}
