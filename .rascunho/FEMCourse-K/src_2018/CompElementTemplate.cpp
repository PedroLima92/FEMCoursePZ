//
//  CompElementTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElement.h"
#include "CompElementTemplate.h"
#include "CompMesh.h"
#include "DataTypes.h"

#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"

#include "GeoElement.h"
#include "GeoElementSide.h"

#include "MathStatement.h"

// Default constructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() : dofindexes(0), intrule(0){
    
}

// Constructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo) : CompElement(ind,cmesh,geo){
    
    this->SetIntRule(&intrule);
    intrule.SetOrder(cmesh->GetDefaultOrder()*2);
    
    geo->SetReference(this);
    
    int nsides = geo->NSides();
    std::vector<DOF> dofvec = cmesh->GetDOFVec();
    this->SetNDOF(nsides);
    int64_t dofid = dofvec.size();
    
    for (int isides=0; isides<nsides; isides++) {
        
        GeoElementSide neighbour = geo->Neighbour(isides);
        GeoElementSide thisside(geo, isides);
        
        DOF celdof;
        int order = cmesh->GetDefaultOrder();
        int nstate = this->GetStatement()->NState();
        if (neighbour != thisside) {
            if (neighbour.Element()->GetIndex() < ind) {
                CompElement * cel = neighbour.Element()->GetReference();
                this->SetDOFIndex(isides, cel->GetDOFIndex(neighbour.Side()));
            }
            else{
                this->SetDOFIndex(isides, dofid);
                celdof.SetNShapeStateOrder(Shape::NShapeFunctions(isides,order), nstate, order);
                dofvec.push_back(celdof);
                dofid++;
            }
        }
        else {
            neighbour.Element()->GetIndex();
            this->SetDOFIndex(isides, dofid);
            celdof.SetNShapeStateOrder(Shape::NShapeFunctions(isides,order), nstate, order);
            dofvec.push_back(celdof);
            dofid++;
        }
    }
    cmesh->SetDOFVec(dofvec);
}

// Copy constructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &copy) : CompElement(copy){
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
}

// Operator of copy
template<class Shape>
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &copy){
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
    return *this;
}

// Destructor of CompElementTemplate
template<class Shape>
CompElementTemplate<Shape>::~CompElementTemplate(){
    
}

// Method for creating a copy of the element
template<class Shape>
CompElement *CompElementTemplate<Shape>::Clone() const{
//    CompElement * cel = new CompElementTemplate<Shape>(*this);
//    return cel;
}

// Compute shape functions set at point x
template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const{
    int nsides = this->GetGeoElement()->NSides();
    VecInt orders(nsides, 0.);
    for (int isides = 0; isides<nsides; isides++) {
        int order = GetCompMesh()->GetDOF(dofindexes[isides]).GetOrder();
        orders[isides] = order;
    }
    Shape::Shape(intpoint, orders, phi, dphi);
}

// Get Multiplying Coeficients
template<class Shape>
void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble &coefs) const{
    
    int64_t nstate = GetStatement()->NState();
    
    coefs.resize(0);
    
    CompMesh * cmesh = GetCompMesh();
    VecDouble sol = cmesh->Solution();
    
    for (int64_t idof=0; idof<NDOF(); idof++) {
        int64_t indexDOF = GetDOFIndex(idof);
        int64_t indexSOL = cmesh->GetDOF(indexDOF).GetFirstEquation();
        int64_t nshape = cmesh->GetDOF(indexDOF).GetNShape();
        
        for (int64_t istate=0; istate<nstate*nshape; istate++) {
            if (indexSOL+istate < sol.size()) {
                coefs.push_back(sol[indexSOL+istate]);
            }
        }
    }
    
}

// Return the number of shape functions
template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const{
    int nsides = GetGeoElement()->NSides();
    VecInt orders(nsides, 0);
    for (int idofs=0; idofs<NDOF(); idofs++) {
        DOF dof = GetCompMesh()->GetDOF(GetDOFIndex(idofs));
        orders[idofs] = dof.GetOrder();
    }
    return Shape::NShapeFunctions(orders);
}

template<class Shape>
void CompElementTemplate<Shape>::SetNDOF(int64_t ndof){
    dofindexes.resize(ndof, 0.);
}

// Se DOF index in vector position i
template<class Shape>
void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex){
    dofindexes[i] = dofindex;
}

// Se DOF index in vector position i
template<class Shape>
int64_t CompElementTemplate<Shape>::GetDOFIndex(int i) const { //MUDEI
    return dofindexes[i];
}

// Return the number of degree of freedom
template<class Shape>
int CompElementTemplate<Shape>::NDOF() const{
    return dofindexes.size();
}

// Return the number of shape functions stored in the DOF data structure
template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const{
    CompMesh * compmesh = GetCompMesh();
    return compmesh->GetDOF(doflocindex).GetNShape();
}

// Use the Shape template class to compute the number of shape functions
template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order){
    VecInt orders(1,order);
    return Shape::NShapeFunctions(orders);
}

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;
