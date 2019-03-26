//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//


#include "MathStatement.h"
#include "DataTypes.h"
#include "IntPointData.h"
#include "L2Projection.h"

// Default constructor of L2Projection
L2Projection::L2Projection() : BCType(0), projection(), BCVal1(), BCVal2(), forceFunction(0), SolutionExact(0){
    
}

// Constructor of L2Projection
L2Projection::L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2) : MathStatement()
{
    BCType = bctype;
    SetMatID(materialid);
    projection = proj;
    BCVal1 = Val1;
    BCVal2 = Val2;
    forceFunction = 0;
    SolutionExact = 0;
}

// Copy constructor of L2Projection
L2Projection::L2Projection(const L2Projection &copy) : MathStatement(copy)
{
    this->operator=(copy);
}

// Operator of copy
L2Projection &L2Projection::operator=(const L2Projection &copy){
    BCType = copy.BCType;
    projection = copy.projection;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    return *this;
}

// Method for creating a copy of the element
L2Projection *L2Projection::Clone() const{
//    return new L2Projection(*this);
}

// Default contructor of L2Projection
L2Projection::~L2Projection(){
    
}

// Return the L2 projection matrix
Matrix L2Projection::GetProjectionMatrix() const{
    return projection;
}

// Set the L2 projection matrix
void L2Projection::SetProjectionMatrix(const Matrix &proj){
    projection = proj;
}

int L2Projection::NEvalErrors() const{
    return 3;
};

// Method to implement integral over element's volume
void L2Projection::Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const{
    
    VecDouble  &phi = integrationpointdata.phi;
    Matrix dsol_u = integrationpointdata.dsoldx;
    
    int64_t nphis = phi.size();
    VecDouble v1(3,0.);
    VecDouble v2(3,0.);
    for (int icols=0; icols<Val1().Cols(); icols++) {
        v1[icols]=Val1()(0,icols);
    }
    for (int icols=0; icols<Val2().Cols(); icols++) {
        v2[icols]=Val2()(1,icols);
    }
    
    //	Dirichlet condition for each state variable
    int boundarycond = GetBCType();
    switch (boundarycond) {
        case 0:
            if (SolutionExact) {
                Matrix deriv(2,3);
                SolutionExact(integrationpointdata.x, v2, deriv);
            }
            for(int64_t in = 0 ; in < nphis; in++)
            {
                //	Contribution for load Vector
                for (int istate=0; istate<NState(); istate++) {
                    EF(NState()*in+istate,0) += gBigNumber*phi[in]*weight*v2[istate];
                }
                
                for (int64_t jn = 0 ; jn < nphis; jn++)
                {
                    //	Contribution for Stiffness Matrix
                    for (int istate=0; istate<NState(); istate++) {
                        EK(in*NState()+istate,jn*NState()+istate) += gBigNumber*phi[in]*phi[jn]*weight;
                    }
                }
            }
            break;
            
        case 1:
            if (forceFunction) {
                forceFunction(integrationpointdata.x, v2);
            }
            for(int64_t in = 0; in < nphis; in++)
            {
                for (int istate=0; istate<NState(); istate++) {
                    EF(in*NState()+istate,0) += v2[istate]*phi[in]*weight;
                }
            }
            break;
            
        default:
            break;
    }
    
    
}

int L2Projection::VariableIndex(PostProcVar var) const{
    switch (var) {
        case ESol:
            return 0;
            break;
            
        case EDSol:
            return 1;
            
        default:
            std::cout << "L2Projection::VariableIndex : Please choose one of the PostProcVar" << std::endl;
            DebugStop();
            return -1;
            break;
    }
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string &name){
    if (!strcmp("ESol", name.c_str())) return ESol;
    if (!strcmp("EDSol", name.c_str())) return EDSol;
}

// Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
int L2Projection::NSolutionVariables(const L2Projection::PostProcVar var){\
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->NState();
}

// Method to implement error over element's volume
void L2Projection::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
    std::cout << "L2Projection::ContributeError : Not implemented yet " << std::endl;
    DebugStop();
}


// Prepare and print post processing data
void L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const{
    std::cout << "L2Projection::PostProcessSolution : Not implemented yet " << std::endl;
    DebugStop();
}

