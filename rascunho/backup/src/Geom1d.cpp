/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"
#include "tpanic.h"

Geom1d::Geom1d() {
    DebugStop();
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    DebugStop();
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    DebugStop();
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    DebugStop();
}

void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    DebugStop();
}

void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    DebugStop();
}

void Geom1d::SetNodes(const VecInt &nodes) {
    DebugStop();
}

void Geom1d::GetNodes(VecInt &nodes) {
    DebugStop();
}

int Geom1d::NodeIndex(int node) {
    DebugStop();
}

int Geom1d::NumNodes() {
    DebugStop();
}

GeoElementSide Geom1d::Neighbour(int side) {
    DebugStop();
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    DebugStop();
}
