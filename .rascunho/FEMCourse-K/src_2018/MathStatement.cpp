//
//  MathStatement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "MathStatement.h"
#include "DataTypes.h"
#include "IntPointData.h"

// Constructor of MathStatement
MathStatement::MathStatement() : MathDim(0), matid(0){
    
}

// Copy constructor of MathStatement
MathStatement::MathStatement(const MathStatement &copy)
{
    this->operator=(copy);
}

// Operator of copy
MathStatement &MathStatement::operator=(const MathStatement &copy){ // MUDEI
    matid = copy.matid;
    MathDim = copy.MathDim;
    return * this;
}

// Destructor of MathStatement
MathStatement::~MathStatement(){
    
}

void MathStatement::Axes2XYZ(const Matrix &dudaxes, Matrix &dudx, const Matrix &axesv, bool colMajor) const{
    
    Matrix axes(axesv.Rows(),axesv.Cols());
    for (int r=0; r<axes.Rows(); r++) {
        for (int c=0; c<axes.Cols(); c++) {
            axes(r,c) = axesv.GetVal(r,c);
        }
    }
    
    if( colMajor ){
        Matrix axesT;
        axes.Transpose(axesT);
        
        if(dudx.Rows() != axesT.Rows() || dudx.Cols() != dudaxes.Cols())
        {
            dudx.Resize(axesT.Rows(), dudaxes.Cols());
        }
        dudx.Zero();
        axesT.Multiply(dudaxes,dudx,0);
    }
    else{
        dudx.Resize(dudaxes.Rows(), axes.Cols());
        dudx.Zero();
        dudaxes.Multiply(axes,dudx,0);
    }
}

void MathStatement::Print(std::ostream &out){
    out << "Material: " << matid << std::endl;
    out << "Dimension: " << MathDim << std::endl;
    out << "Big number: " << gBigNumber << std::endl << std::endl;
}

double MathStatement::gBigNumber = 1e12;
