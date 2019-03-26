/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompElement.h"
#include "tpanic.h"

using namespace std;

CompElement::CompElement() {
    DebugStop();
}

CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo) {
    DebugStop();
}

CompElement::CompElement(const CompElement &copy) {
    DebugStop();
}

CompElement &CompElement::operator=(const CompElement &copy) {
    DebugStop();
}

CompElement::~CompElement() {
}

CompElement *CompElement::Clone() const {
    DebugStop();
}

MathStatement *CompElement::GetStatement() const {
    DebugStop();
}

void CompElement::SetStatement(MathStatement *statement) {
    DebugStop();
}

IntRule *CompElement::GetIntRule() const {
    DebugStop();
}

void CompElement::SetIntRule(IntRule *irule) {
    DebugStop();
}

void CompElement::SetIndex(int64_t ind) {
    DebugStop();
}

GeoElement *CompElement::GetGeoElement() const {
    DebugStop();
}

void CompElement::SetGeoElement(GeoElement *element) {
    geoel = element;
}

CompMesh *CompElement::GetCompMesh() const {
    DebugStop();
}

void CompElement::SetCompMesh(CompMesh *mesh) {
    DebugStop();
}

void CompElement::InitializeIntPointData(IntPointData &data) const {
    DebugStop();
}

void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const {
    DebugStop();
}

void CompElement::Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const {
    DebugStop();
}

void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const {
    DebugStop();
}

void CompElement::EvaluateError(std::function<void(const VecDouble &loc, VecDouble &val, Matrix &deriv) > fp, VecDouble &errors) const {
    DebugStop();
}

void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol) const {
    DebugStop();
}
