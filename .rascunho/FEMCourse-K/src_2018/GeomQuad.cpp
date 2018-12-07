//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "GeomQuad.h"

 /// Constructor
GeomQuad::GeomQuad() : fNodeIndices(0){
    
}

/// destructor
GeomQuad::~GeomQuad(){
    
}

/// copy constructor
GeomQuad::GeomQuad(const GeomQuad &copy){
    this->operator=(copy);
}

/// operator=
GeomQuad &GeomQuad::operator=(const GeomQuad &copy){
    fNodeIndices = copy.fNodeIndices;
    for (int isides=0; isides<nSides; isides++) {
        fNeighbours[isides] = copy.fNeighbours[isides];
    }
    return *this;
}

/// Computes the shape functions associated with the geometric map
void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
    
    phi[0] = 0.25*(1.-xi[0])*(1.-xi[1]);
    phi[1] = 0.25*(1.+xi[0])*(1.-xi[1]);
    phi[2] = 0.25*(1.+xi[0])*(1.+xi[1]);
    phi[3] = 0.25*(1.-xi[0])*(1.+xi[1]);
    
    dphi(0,0) = 0.25*(xi[1]-1.);
    dphi(1,0) = 0.25*(xi[0]-1.);
    
    dphi(0,1) = 0.25*(1.-xi[1]);
    dphi(1,1) =-0.25*(1.+xi[0]);
    
    dphi(0,2) = 0.25*(1.+xi[1]);
    dphi(1,2) = 0.25*(1.+xi[0]);
    
    dphi(0,3) =-0.25*(1.+xi[1]);
    dphi(1,3) = 0.25*(1.-xi[0]);
}

/// Computes the value of x for a given point in parameter space as a function of corner coordinates
void GeomQuad::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
    
    VecDouble phi(nCorners,0.);
    Matrix dphi(2,nCorners);
    Shape(xi,phi,dphi);
    int space = NodeCo.Rows();
    
    for(int i = 0; i < space; i++) {
        x[i] = 0.0;
        for(int j = 0; j < 4; j++) {
            x[i] += phi[j]*NodeCo.GetVal(i,j);
        }
    }
}

/// Computes the value of x and gradx for a given point in parameter space
void GeomQuad::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
    
    int nrow = NodeCo.Rows();
    gradx.Zero();
    
    VecDouble phi(nCorners);
    Matrix dphi(2,nCorners);
    Shape(xi,phi,dphi);
    
    for(int i = 0; i < dphi.Cols(); i++){
        for(int j = 0; j < nrow; j++)
        {
            gradx(0,j) += NodeCo.GetVal(j,i)*dphi(0,i);
            gradx(1,j) += NodeCo.GetVal(j,i)*dphi(1,i);
        }
    }
}

int GeomQuad::NumNodes(){
    return nCorners;
}

/// Set the node indices of the element
void GeomQuad::SetNodes(const VecInt &nodes){
    fNodeIndices = nodes;
}

/// Set the node indices of the element
void GeomQuad::GetNodes(VecInt &nodes){
    nodes = fNodeIndices;
}

/// Return the index of a node
int GeomQuad::NodeIndex(int node){
    return fNodeIndices[node];
}

GeoElementSide GeomQuad::Neighbour(int side){
    return fNeighbours[side];
}

// Initialize the neighbour data structure
void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
