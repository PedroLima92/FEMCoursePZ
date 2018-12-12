/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "CompElementTemplate.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"
#include "CompElementTemplate.h"
#include "tpanic.h"

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() {
    DebugStop();
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo) {
    DebugStop();
}

template<class Shape >
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate & copy) {
    DebugStop();
}

template<class Shape >
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate & copy) {
    DebugStop();
}

template<class Shape >
CompElementTemplate<Shape>::~CompElementTemplate() {
}

template<class Shape >
CompElement * CompElementTemplate<Shape>::Clone() const {
    DebugStop();
}

template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix & dphi) const {
    DebugStop();
}

template<class Shape>
void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble & coefs) const {
    DebugStop();
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const {
    DebugStop();
}

template<class Shape>
void CompElementTemplate<Shape>::SetNDOF(int64_t ndof) {
    DebugStop();
}

template<class Shape>
void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex) {
    DebugStop();
}

template<class Shape>
int64_t CompElementTemplate<Shape>::GetDOFIndex(int i) {
    DebugStop();
}

template<class Shape>
int CompElementTemplate<Shape>::NDOF() const {
    DebugStop();
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const {
    DebugStop();
}

template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order) {
    DebugStop();
}

template<class Shape>
void CompElementTemplate<Shape>::Print(std::ostream &out) {
    DebugStop();
}



template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
