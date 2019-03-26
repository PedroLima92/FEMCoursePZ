/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include "tpanic.h"

Poisson::Poisson() {
    DebugStop();
}

Poisson::Poisson(int materialid, Matrix &perm) {
    DebugStop();
}

Poisson::Poisson(const Poisson &copy) {
    DebugStop();
}

Poisson &Poisson::operator=(const Poisson &copy) {
    DebugStop();
}

Poisson *Poisson::Clone() const {
    DebugStop();
}

Poisson::~Poisson() {
}

Matrix Poisson::GetPermeability() const {
    DebugStop();
}

void Poisson::SetPermeability(const Matrix &perm) {
    DebugStop();
}

int Poisson::NEvalErrors() const {
    DebugStop();
}

int Poisson::VariableIndex(const PostProcVar var) const {
    DebugStop();
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    DebugStop();
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    DebugStop();
}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    DebugStop();
}

void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    DebugStop();
}

void Poisson::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    DebugStop();
}

double Poisson::Inner(Matrix &S, Matrix & T) const {
    DebugStop();
}
