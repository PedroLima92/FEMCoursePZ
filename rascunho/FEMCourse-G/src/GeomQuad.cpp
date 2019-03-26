/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"
#include "tpanic.h"

GeomQuad::GeomQuad() {
    fNodeIndices = {-1, -1, -1, -1};
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;

    for (auto i = 0; i < nSides; i++) {
        fNeighbours[i] = copy.fNeighbours[i];
    }
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;

    for (auto i = 0; i < nSides; i++) {
        fNeighbours[i] = copy.fNeighbours[i];
    }
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    phi[0] = (1. - xi[0]) * (1. - xi[1]) / 4.;
    phi[1] = (1. + xi[0]) * (1. - xi[1]) / 4.;
    phi[2] = (1. + xi[0]) * (1. + xi[1]) / 4.;
    phi[3] = (1. - xi[0]) * (1. + xi[1]) / 4.;

    dphi(0, 0) = (+xi[1] - 1.) / 4.;
    dphi(0, 1) = (-xi[1] + 1.) / 4.;
    dphi(0, 2) = (+xi[1] + 1.) / 4.;
    dphi(0, 3) = (-xi[1] - 1.) / 4.;

    dphi(1, 0) = (+xi[0] - 1.) / 4.;
    dphi(1, 1) = (-xi[0] - 1.) / 4.;
    dphi(1, 2) = (+xi[0] + 1.) / 4.;
    dphi(1, 3) = (-xi[0] + 1.) / 4.;
}

void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
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

void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
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

void GeomQuad::SetNodes(const VecInt &nodes) {
    if (nodes.size() != nCorners) {
        DebugStop();
    }

    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) {
    if (node < 0 || node > (nCorners - 1)) {
        DebugStop();
    }

    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) {
    if (side < 0 || (side > nSides - 1)) {
        DebugStop();
    }

    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    if (side < 0 || (side > nSides - 1)) {
        DebugStop();
    }

    fNeighbours[side] = neighbour;
}
