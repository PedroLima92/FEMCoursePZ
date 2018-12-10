//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "GeomTetrahedron.h"

 /// Constructor
GeomTetrahedron::GeomTetrahedron() : fNodeIndices(0){
    
}

/// destructor
GeomTetrahedron::~GeomTetrahedron(){
    
}

/// copy constructor
GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy){
    this->operator=(copy);
}

/// operator=
GeomTetrahedron &GeomTetrahedron::operator=(const GeomTetrahedron &copy){
    fNodeIndices = copy.fNodeIndices;
    for (int isides=0; isides<nSides; isides++) {
        fNeighbours[isides] = copy.fNeighbours[isides];
    }
    return *this;
}

/// Computes the shape functions associated with the geometric map
void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi){
    phi[0] = 1.0-xi[0]-xi[1]-xi[2];
    phi[1] = xi[0];
    phi[2] = xi[1];
    phi[3] = xi[2];
    
    dphi(0,0) = -1.0;
    dphi(1,0) = -1.0;
    dphi(2,0) = -1.0;
    dphi(0,1) =  1.0;
    dphi(1,1) =  0.0;
    dphi(2,1) =  0.0;
    dphi(0,2) =  0.0;
    dphi(1,2) =  1.0;
    dphi(2,2) =  0.0;
    dphi(0,3) =  0.0;
    dphi(1,3) =  0.0;
    dphi(2,3) =  1.0;
}

/// Computes the value of x for a given point in parameter space as a function of corner coordinates
void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x){
    
//    int nrow = NodeCo.Rows();
//    for(int i = 0; i < nrow; i++)
//        x[i] = NodeCo.GetVal(i,0)*(1.-xi[i])*0.5+NodeCo.GetVal(i,1)*(1.+xi[i])*0.5;
    VecDouble phi(4);
    Matrix dphi(3,4);
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
void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx){
    
    int nrow = NodeCo.Rows();
    gradx.Zero();
    
    VecDouble phi(nCorners);
    Matrix dphi(3,nCorners);
    Shape(xi,phi,dphi);
    
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 3; j++)
        {
            gradx(j,0) += NodeCo.GetVal(j,i)*dphi(0,i);
            gradx(j,1) += NodeCo.GetVal(j,i)*dphi(1,i);
            gradx(j,2) += NodeCo.GetVal(j,i)*dphi(2,i);
        }
    }
}

/// Set the node indices of the element
int GeomTetrahedron::NumNodes(){
    return nCorners;
}

/// Set the node indices of the element
void GeomTetrahedron::GetNodes(VecInt &nodes){
    nodes = fNodeIndices;
}

void GeomTetrahedron::SetNodes(const VecInt &nodes){
    fNodeIndices = nodes;
}

/// Return the index of a node
int GeomTetrahedron::NodeIndex(int node){
    return fNodeIndices[node];
}

GeoElementSide GeomTetrahedron::Neighbour(int side){
    return fNeighbours[side];
}

// Initialize the neighbour data structure
void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour){
    fNeighbours[side] = neighbour;
}
