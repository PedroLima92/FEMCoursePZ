//
//  GeoElementTemplate.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/28/18.
//
//

#include <stdio.h>
#include "GeoElementTemplate.h"

#include "GeoElement.h"
#include "GeoNode.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeoMesh.h"

#include "tpanic.h"

/// constructor
template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) : GeoElement(materialid, gmesh, index){
    Geom.SetNodes(nodeindices);
    for(int side=0; side<TGeom::nSides; side++)
    {
        Geom.SetNeighbour(side,GeoElementSide(this,side));
    }
}

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy){
    Geom = copy.Geom;
}

template<class TGeom>
GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy){
    Geom = copy.Geom;
    return * this;
}

/// return the enumerated element type
template<class TGeom>
ElementType GeoElementTemplate<TGeom>::Type(){
    return Geom.Type();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x){
    Matrix NodeCoor(3,NNodes());
    
    int i,j;
    for(i=0; i<this->NNodes(); i++) {
        int id = NodeIndex(i);
        for(j=0;j<3;j++) {
            NodeCoor(j,i) = GMesh->Node(id).Coord(j);
        }
    }
    
    Geom.X(xi, NodeCoor, x);
};

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx){
    Matrix NodeCoor(3,this->NNodes(),0.);
    
    int i,j;
    for(i=0; i<this->NNodes(); i++) {
        int id = NodeIndex(i);
        for(j=0;j<3;j++) {
            NodeCoor(j,i) = GMesh->Node(id).Coord(j);
        }
    }
    
    Geom.GradX(xi, NodeCoor, x, gradx);
    
};

template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac, Matrix &axes, double &detjac, Matrix &jacinv){
    
    detjac = 0.0;
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim   = GetReference()->Dimension();
    
    switch (dim) {
        case 0:
        {
            detjac=1.;
            break;
        }
            
        case 1:
        {
            jac.Zero();
            
            VecDouble v_1(3,0.);
            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
            }
            
            double norm_v_1 = 0.;
            for(int i = 0; i < nrows; i++) {
                norm_v_1 += v_1[i]*v_1[i];
            }
            norm_v_1    = sqrt(norm_v_1);
            jac(0,0)    = norm_v_1;
            detjac      = norm_v_1;
            jacinv(0,0) = 1.0/detjac;
            
            detjac = fabs(detjac);
            
            for(int i=0; i < 3; i++) {
                axes(0,i) = v_1[i]/norm_v_1;
            }
            
        }
            break;
        case 2:
        {
            jac.Zero();
            
            VecDouble v_1(3,0.), v_2(3,0.);
            
            /**  Definitions: v_1_til and v_2_til -> asscoiated orthonormal vectors to v_1 and v_2 */
            VecDouble v_1_til(3,0.), v_2_til(3,0.);
            
            for (int i = 0; i < nrows; i++) {
                v_1[i]  = gradx.GetVal(i,0);
                v_2[i]  = gradx.GetVal(i,1);
            }
            
            double norm_v_1_til = 0.0;
            double norm_v_2_til = 0.0;
            double v_1_dot_v_2  = 0.0;
            for(int i = 0; i < 3; i++) {
                norm_v_1_til    += v_1[i]*v_1[i];
                v_1_dot_v_2     += v_1[i]*v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);
            
            for(int i=0 ; i < 3; i++) {
                v_1_til[i]          = v_1[i] / norm_v_1_til; // Normalizing
                v_2_til[i]          = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til   += v_2_til[i]*v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);
            
            
            jac(0,0) = norm_v_1_til;
            jac(0,1) = v_1_dot_v_2/norm_v_1_til;
            jac(1,1) = norm_v_2_til;
            
            detjac = jac(0,0)*jac(1,1)-jac(1,0)*jac(0,1);
            
            jacinv(0,0) = +jac(1,1)/detjac;
            jacinv(1,1) = +jac(0,0)/detjac;
            jacinv(0,1) = -jac(0,1)/detjac;
            jacinv(1,0) = -jac(1,0)/detjac;
            
            detjac = fabs(detjac);
            
            for(int i=0; i < 3; i++) {
                v_2_til[i] /= norm_v_2_til; // Normalizing
                axes(0,i)  = v_1_til[i];
                axes(1,i)  = v_2_til[i];
            }
            
        }
            break;
        case 3:
        {
            jac.Zero();
            
            for (int i = 0; i < ncols; i++) {
                jac(i,0)  = gradx.GetVal(i,0);
                jac(i,1)  = gradx.GetVal(i,1);
                jac(i,2)  = gradx.GetVal(i,2);
            }
            
            detjac -= jac(0,2)*jac(1,1)*jac(2,0);//- a02 a11 a20
            detjac += jac(0,1)*jac(1,2)*jac(2,0);//+ a01 a12 a20
            detjac += jac(0,2)*jac(1,0)*jac(2,1);//+ a02 a10 a21
            detjac -= jac(0,0)*jac(1,2)*jac(2,1);//- a00 a12 a21
            detjac -= jac(0,1)*jac(1,0)*jac(2,2);//- a01 a10 a22
            detjac += jac(0,0)*jac(1,1)*jac(2,2);//+ a00 a11 a22
            
            jacinv(0,0) = (-jac(1,2)*jac(2,1)+jac(1,1)*jac(2,2))/detjac;//-a12 a21 + a11 a22
            jacinv(0,1) = ( jac(0,2)*jac(2,1)-jac(0,1)*jac(2,2))/detjac;//a02 a21 - a01 a22
            jacinv(0,2) = (-jac(0,2)*jac(1,1)+jac(0,1)*jac(1,2))/detjac;//-a02 a11 + a01 a12
            jacinv(1,0) = ( jac(1,2)*jac(2,0)-jac(1,0)*jac(2,2))/detjac;//a12 a20 - a10 a22
            jacinv(1,1) = (-jac(0,2)*jac(2,0)+jac(0,0)*jac(2,2))/detjac;//-a02 a20 + a00 a22
            jacinv(1,2) = ( jac(0,2)*jac(1,0)-jac(0,0)*jac(1,2))/detjac;//a02 a10 - a00 a12
            jacinv(2,0) = (-jac(1,1)*jac(2,0)+jac(1,0)*jac(2,1))/detjac;//-a11 a20 + a10 a21
            jacinv(2,1) = ( jac(0,1)*jac(2,0)-jac(0,0)*jac(2,1))/detjac;//a01 a20 - a00 a21
            jacinv(2,2) = (-jac(0,1)*jac(1,0)+jac(0,0)*jac(1,1))/detjac;//-a01 a10 + a00 a11
            
            detjac = fabs(detjac);
            
            axes.Zero();
            axes(0,0) = 1.0;
            axes(1,1) = 1.0;
            axes(2,2) = 1.0;
            
        }
            break;
            
        default:
        {
            std::cout << " CompElement::Jacobian : Wrong dimension" << std::endl;
            DebugStop();
        }
            break;
    }
    
}

template<class TGeom>
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds){
    std::cout << "GeoElementTemplate::WhichSide : Not implemented yet" << std::endl;
    DebugStop();
    return 0;
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out){
    out << "Element index    " << GetIndex() << std::endl;
    out << "Material index    " << Material() << std::endl;
    out << "Number of nodes    " << NNodes() << std::endl;
    out << "Corner nodes       " << NCornerNodes() << std::endl;
    out << "Nodes indexes          ";
    int i;
    for (i = 0;i < NNodes();i++) out << NodeIndex(i) << " ";
    out << "\nNumber of sides    " << NSides() << std::endl;
    out << std::endl;
    for (i = 0;i < NSides();i++) {
        out << "Neighbours for side   " << i << " : ";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this,i);
        if (!(neighbour.Element() != 0 && neighbour.Side() > -1)) {
            out << "No neighbour\n";
        }
        else {
            while ((neighbour == thisside) == false) {
                out << neighbour.Element()->GetIndex() <<  "/" << neighbour.Side() << ' ';
                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
    }
};

template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<GeomTetrahedron>;

