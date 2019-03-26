/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElementTemplate.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "tpanic.h"

using namespace std;

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) {
    DebugStop();
}

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy) {
    DebugStop();
}

template<class TGeom>
GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy) {
    DebugStop();
}

template<class TGeom>
ElementType GeoElementTemplate<TGeom>::Type() {
    DebugStop();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x) {
    DebugStop();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx) {
    DebugStop();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac, Matrix &axes, double &detjac, Matrix &jacinv) {
    DebugStop();
}

template<class TGeom>
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds) {
    DebugStop();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out) {
    DebugStop();
}

template<class TGeom>
int GeoElementTemplate<TGeom>::SideIsUndefined(int side) {
    DebugStop();
}


template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<GeomTetrahedron>;
