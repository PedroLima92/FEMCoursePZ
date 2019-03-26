/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTriangle.h"
#include "tpanic.h"

GeomTriangle::GeomTriangle() {
    DebugStop();
}

GeomTriangle::~GeomTriangle() {
}

GeomTriangle::GeomTriangle(const GeomTriangle &copy) {
    DebugStop();

}

GeomTriangle& GeomTriangle::operator=(const GeomTriangle& copy) {
    DebugStop();
}

void GeomTriangle::Shape(const VecDouble& xi, VecDouble& phi, Matrix& dphi) {
    DebugStop();
}

void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    DebugStop();
}

void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    DebugStop();
}

void GeomTriangle::SetNodes(const VecInt &nodes) {
    DebugStop();
}

void GeomTriangle::GetNodes(VecInt &nodes) {
    DebugStop();
}

int GeomTriangle::NodeIndex(int node) {
    DebugStop();
}

int GeomTriangle::NumNodes() {
    DebugStop();
}

GeoElementSide GeomTriangle::Neighbour(int side) {
    DebugStop();
}

void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour) {
    DebugStop();
}
