/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "DOF.h"
#include "tpanic.h"

DOF::DOF() {
    DebugStop();
}

DOF::DOF(const DOF &copy) {
    DebugStop();
}

DOF &DOF::operator=(const DOF &copy) {
    DebugStop();
}

DOF::~DOF() {
}

int64_t DOF::GetFirstEquation() {
    DebugStop();
}

void DOF::SetFirstEquation(int64_t first) {
    DebugStop();
}

void DOF::SetNShapeStateOrder(int NShape, int NState, int Order) {
    DebugStop();
}

int DOF::GetNShape() const {
    DebugStop();
}

int DOF::GetNState() const {
    DebugStop();
}

int DOF::GetOrder() const {
    DebugStop();
}

void DOF::Print(const CompMesh &mesh, std::ostream & out) {
    DebugStop();
}
