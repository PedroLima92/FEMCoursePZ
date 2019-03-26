/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTetrahedron.h"
#include "tpanic.h"

GeomTetrahedron::GeomTetrahedron() {
    fNodeIndices = {-1, -1, -1, -1};
}

GeomTetrahedron::~GeomTetrahedron() {
}

GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy) {
    fNodeIndices = copy.fNodeIndices;

    for (auto i = 0; i < nSides; i++) {
        fNeighbours[i] = copy.fNeighbours[i];
    }
}

GeomTetrahedron& GeomTetrahedron::operator=(const GeomTetrahedron& copy){
    fNodeIndices = copy.fNodeIndices;

    for (auto i = 0; i < nSides; i++) {
        fNeighbours[i] = copy.fNeighbours[i];
    }
    return *this;
}

void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    phi[0] = 1. - xi[0] - xi[1] - xi[2];
    phi[1] = xi[0];
    phi[2] = xi[1];
    phi[3] = xi[2];

    dphi(0, 0) = -1.;
    dphi(0, 1) = 1.;
    dphi(0, 2) = 0.;
    dphi(0, 3) = 0.;

    dphi(1, 0) = -1.;
    dphi(1, 1) = 0.;
    dphi(1, 2) = 1.;
    dphi(1, 3) = 0.;

    dphi(2, 0) = -1.;
    dphi(2, 1) = 0.;
    dphi(2, 2) = 0.;
    dphi(2, 3) = 1.;
}

void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    VecDouble phi(nCorners, 0.);
    Matrix dphi(nCorners, nCorners, 0.);

    Shape(xi, phi, dphi);
    int dimensions = NodeCo.Rows();

    for(int i = 0; i < dimensions; i++) {
        x[i] = 0.0;
        for(int j = 0; j < nCorners; j++) {
            x[i] += phi[j] * NodeCo.GetVal(i, j);
        }
    }
}

void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    int dimensions = NodeCo.Rows();
    int ncol = NodeCo.Cols();

    gradx.Resize(dimensions, 1);
    gradx.Zero();

    VecDouble phi(2);
    Matrix dphi(2, 2);
    Shape(xi, phi, dphi);

    for (int i = 0; i < dimensions; i++) {
        for (int j = 0; j < ncol; j++) {
            gradx(i, 0) += NodeCo.GetVal(i, j) * dphi(0, j);
        }
    }
}


void GeomTetrahedron::SetNodes(const VecInt &nodes) {
    if (nodes.size() != nCorners) {
        DebugStop();
    }

    fNodeIndices = nodes;
}

void GeomTetrahedron::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomTetrahedron::NodeIndex(int node) {
    if (node < 0 || node > (nCorners - 1)) {
        DebugStop();
    }

    return fNodeIndices[node];
}

int GeomTetrahedron::NumNodes() {
    return nCorners;
}

GeoElementSide GeomTetrahedron::Neighbour(int side) {
    if (side < 0 || (side > nSides - 1)) {
        DebugStop();
    }

    return fNeighbours[side];
}

void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour) {
    if (side < 0 || (side > nSides - 1)) {
        DebugStop();
    }

    fNeighbours[side] = neighbour;
}
