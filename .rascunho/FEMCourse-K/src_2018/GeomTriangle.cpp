//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "GeomTriangle.h"

 /// Constructor
GeomTriangle::GeomTriangle() : fNodeIndices(0){
    
}

/// destructor
GeomTriangle::~GeomTriangle(){
    
}

/// copy constructor
GeomTriangle::GeomTriangle(const GeomTriangle &copy){
    this->operator=(copy);
}

/// operator=
GeomTriangle &GeomTriangle::operator=(const GeomTriangle &copy){
    fNodeIndices = copy.fNodeIndices;
    for (int isides=0; isides<nSides; isides++) {
        fNeighbours[isides] = copy.fNeighbours[isides];
    }
    return *this;
}

/// Computes the shape functions associated with the geometric map
void GeomTriangle::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
    double qsi = xi[0], eta = xi[1];
    phi[0] = 1.-qsi-eta;
    phi[1] = qsi;
    phi[2] = eta;
    dphi(0,0) = dphi(1,0) = -1.;
    dphi(0,1) = dphi(1,2) =  1.;
    dphi(1,1) = dphi(0,2) =  0.;
}

/// Computes the value of x for a given point in parameter space as a function of corner coordinates
void GeomTriangle::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
    
//    int nrow = NodeCo.Rows();
//    for(int i = 0; i < nrow; i++)
//        x[i] = NodeCo.GetVal(i,0)*(1.-xi[i])*0.5+NodeCo.GetVal(i,1)*(1.+xi[i])*0.5;
    
    VecDouble phi(nCorners,0.);
    Matrix dphi(2,nCorners);
    Shape(xi,phi,dphi);
    int space = NodeCo.Rows();
    
    for(int i = 0; i < space; i++) {
        x[i] = 0.0;
        for(int j = 0; j < nCorners; j++) {
            x[i] += phi[j]*NodeCo.GetVal(i,j);
        }
    }
    
}

/// Computes the value of x and gradx for a given point in parameter space
void GeomTriangle::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
    
    int space = gradx.Cols();
    int nrow = gradx.Rows();
    gradx.Zero();
    
    VecDouble phi(nCorners);
    Matrix dphi(2,nCorners);
    Shape(xi,phi,dphi);
    
    for(int i = 0; i < space; i++){
        for(int j = 0; j < nrow; j++)
        {
            gradx(j,0) += NodeCo.GetVal(j,i)*dphi(0,i);
            gradx(j,1) += NodeCo.GetVal(j,i)*dphi(1,i);
        }
    }
}

/// Set the node indices of the element
int GeomTriangle::NumNodes(){
    return nCorners;
}

/// Set the node indices of the element
void GeomTriangle::SetNodes(const VecInt &nodes){
    fNodeIndices = nodes;
}

/// Set the node indices of the element
void GeomTriangle::GetNodes(VecInt &nodes){
    nodes = fNodeIndices;
}

/// Return the index of a node
int GeomTriangle::NodeIndex(int node){
    return fNodeIndices[node];
}

GeoElementSide GeomTriangle::Neighbour(int side){
    return fNeighbours[side];
}

// Initialize the neighbour data structure
void GeomTriangle::SetNeighbour(int side, const GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
