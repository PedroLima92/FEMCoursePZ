/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompMesh.h"
#include "tpanic.h"

CompMesh::CompMesh() {
    DebugStop();
}

CompMesh::CompMesh(const CompMesh &copy) {
    DebugStop();
}

CompMesh::CompMesh(GeoMesh *gmesh) {
    DebugStop();
}

CompMesh::~CompMesh() {
}

GeoMesh *CompMesh::GetGeoMesh() const {
    DebugStop();
}

void CompMesh::SetGeoMesh(GeoMesh *gmesh) {
    DebugStop();
}

void CompMesh::SetNumberElement(int64_t nelem) {
    DebugStop();
}

void CompMesh::SetNumberDOF(int64_t ndof) {
    DebugStop();
}

void CompMesh::SetNumberMath(int nmath) {
    DebugStop();
}

void CompMesh::SetElement(int64_t elindex, CompElement *cel) {
    DebugStop();
}

void CompMesh::SetDOF(int64_t index, const DOF &dof) {
    DebugStop();
}

void CompMesh::SetMathStatement(int index, MathStatement *math) {
    DebugStop();
}

DOF &CompMesh::GetDOF(int64_t dofindex) {
    DebugStop();
}

CompElement *CompMesh::GetElement(int64_t elindex) const {
    DebugStop();
}

MathStatement *CompMesh::GetMath(int matindex) const {
    DebugStop();
}

std::vector<CompElement *> CompMesh::GetElementVec() const {
    DebugStop();
}

std::vector<DOF> CompMesh::GetDOFVec() const {
    DebugStop();
}

std::vector<MathStatement *> CompMesh::GetMathVec() const {
    DebugStop();
}

void CompMesh::SetElementVec(const std::vector<CompElement *> &vec) {
    DebugStop();
}

void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec) {
    DebugStop();
}

void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec) {
    DebugStop();
}

void CompMesh::AutoBuild() {
    DebugStop();
}

void CompMesh::Resequence() {
    DebugStop();
}

void CompMesh::Resequence(VecInt &DOFindices) {
    DebugStop();
}

std::vector<double> &CompMesh::Solution() {
    DebugStop();
}

void CompMesh::LoadSolution(std::vector<double> &Sol) {
    DebugStop();
}

void CompMesh::Print(std::ostream & out) {
    DebugStop();
}
