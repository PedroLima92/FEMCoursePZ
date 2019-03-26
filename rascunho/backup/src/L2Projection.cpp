/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "tpanic.h"

L2Projection::L2Projection() {
    DebugStop();
}

L2Projection::L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2) {
    DebugStop();
}

L2Projection::L2Projection(const L2Projection &copy) {
    DebugStop();
}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    DebugStop();
}

L2Projection *L2Projection::Clone() const {
    DebugStop();
}

L2Projection::~L2Projection() {
}

Matrix L2Projection::GetProjectionMatrix() const {
    DebugStop();
}

void L2Projection::SetProjectionMatrix(const Matrix &proj) {
    DebugStop();
}

void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    DebugStop();
}

int L2Projection::NEvalErrors() const {
    DebugStop();
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    DebugStop();
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    DebugStop();
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    DebugStop();
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    DebugStop();
}

void L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const {
    DebugStop();
}
