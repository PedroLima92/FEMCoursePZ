/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "tpanic.h"
//Karol mandou fazer

using namespace std;

Analysis::Analysis() {
    DebugStop();

}

Analysis::Analysis(const Analysis &cp) {
    DebugStop();
}

Analysis &Analysis::operator=(const Analysis &cp) {
    DebugStop();
}

Analysis::~Analysis() {
}

Analysis::Analysis(CompMesh *mesh) {
    DebugStop();
}

void Analysis::SetMesh(CompMesh *mesh) {
    DebugStop();
}

CompMesh *Analysis::Mesh() const {
    DebugStop();
}

void Analysis::RunSimulation() {
    DebugStop();
}

void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const {
    DebugStop();
}

VecDouble Analysis::PostProcessError(std::ostream &out, PostProcess &defPostProc) const {
    DebugStop();
}

