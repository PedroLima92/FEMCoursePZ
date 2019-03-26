/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "tpanic.h"

Assemble::Assemble() {
    DebugStop();
}

Assemble::Assemble(CompMesh *mesh) {
    DebugStop();
}

Assemble::Assemble(const Assemble &copy) {
    DebugStop();
}

Assemble &Assemble::operator=(const Assemble &copy) {
    DebugStop();
}

void Assemble::SetMesh(CompMesh *mesh) {
    DebugStop();
}

int64_t Assemble::NEquations() {
    DebugStop();
}

void Assemble::OptimizeBandwidth() {    
}

void Assemble::Compute(Matrix &globmat, Matrix &rhs) {
    DebugStop();
}
