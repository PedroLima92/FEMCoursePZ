/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTetrahedron.h"
#include "tpanic.h"

GeomTetrahedron::GeomTetrahedron() {
    DebugStop();
}

GeomTetrahedron::~GeomTetrahedron() {
}

GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy) {
    DebugStop();
}

GeomTetrahedron& GeomTetrahedron::operator=(const GeomTetrahedron& copy){
    DebugStop();
}

void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    DebugStop();
}

void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    DebugStop();
}


void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    DebugStop();
}


void GeomTetrahedron::SetNodes(const VecInt &nodes) {
    DebugStop();
}

void GeomTetrahedron::GetNodes(VecInt &nodes) {
    DebugStop();
}

int GeomTetrahedron::NodeIndex(int node) {
    DebugStop();
}

int GeomTetrahedron::NumNodes() {
    DebugStop();

}

GeoElementSide GeomTetrahedron::Neighbour(int side) {
    DebugStop();
}

void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour) {
    DebugStop();
}
