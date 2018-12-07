//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "MathStatement.h"
#include "DataTypes.h"
#include "IntPointData.h"
#include "Poisson.h"

 // Default constructor of Poisson
Poisson::Poisson() : permeability(), forceFunction(0), SolutionExact(0){
    
}

// Constructor of Poisson
Poisson::Poisson(int materialid, Matrix &perm) : MathStatement()
{
    permeability = perm;
    SetMatID(materialid);
    forceFunction = 0;
    SolutionExact = 0;
}

// Copy constructor of Poisson
Poisson::Poisson(const Poisson &copy) : MathStatement(copy)
{
    this->operator=(copy);
}

// Operator of copy
Poisson &Poisson::operator=(const Poisson &copy)
{
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    return *this;
}

// Method for creating a copy of the element
Poisson *Poisson::Clone() const{
//    return new Poisson(*this);
}

// Destructor of Poisson
Poisson::~Poisson(){
    
}

// Return the permeability matrix
Matrix Poisson::GetPermeability() const{
    return permeability;
}

// Set the permeability matrix
void Poisson::SetPermeability(const Matrix &perm){
    permeability = perm;
}

int Poisson::NEvalErrors() const{
    return 3;
}

// Method to implement integral over element's volume
void Poisson::Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const{
    
    VecDouble &phi = integrationpointdata.phi;
    Matrix &dphi = integrationpointdata.dphidx;
    VecDouble &x = integrationpointdata.x;
    
    int phr = phi.size();
    
    Matrix perm = GetPermeability();
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        VecDouble ff(NState(),0.);
        if(forceFunction) {
            forceFunction(x, ff);
        }
        for (int istate=0; istate<NState(); istate++) {
            EF(in*NState()+istate, 0) +=  weight * ff[istate] * phi[in];
        }
        
        for(int jn = 0; jn < phr; jn++ ) {
            for(int kd=0; kd<Dimension(); kd++) {
                for (int istate=0; istate<NState(); istate++) {
                    EK(NState()*in+istate,NState()*jn+istate) += weight*(perm(kd,kd)*(dphi(kd,in)*dphi(kd,jn)));
                }
                
            }
        }
    }
//    EK.Print();
//    EF.Print();
}

// Method to implement error over element's volume
void Poisson::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
    
    errors.resize(NEvalErrors(), 0.);
    
    Matrix perm = GetPermeability();
    VecDouble val(Dimension());
    
    VecDouble u = integrationpointdata.solution;
    
    Matrix dudaxes = integrationpointdata.dsoldx;
    Matrix dudx = integrationpointdata.dsoldx;
    Matrix axesv = integrationpointdata.axes;
    this->Axes2XYZ(dudaxes, dudx, axesv);
    
    ///L2 norm
    errors[1] = (u[0] - u_exact[0])*(u[0] - u_exact[0]);
    
    ///semi norma de H1
    errors[2] = 0.;
    for(int i = 0; i < du_exact.Rows(); i++){
        errors[2] += (dudx(i,0) - du_exact(i,0))*(dudx(i,0) - du_exact(i,0));
    }
    
    // Energy Norm
    errors[0] = errors[1]+errors[2];
    ///H1 norm
    
}

// Prepare and print post processing data
void Poisson::PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const{
    VecDouble solint = integrationpointdata.solution;
    Matrix gradu = integrationpointdata.dsoldx;
    
    int nstate = this->NState();
    
    switch (var) {
        case 0: //None
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
            
        case 1: //ESol
        {
            sol = solint;
        }
            break;
            
        case 2: //EDSol
        {
            sol.resize(gradu.Rows() * gradu.Cols());
            for (int i = 0; i < gradu.Rows(); i++) {
                for (int j = 0; j < gradu.Cols(); j++) {
                    sol[i * gradu.Cols() + j] = gradu(i, j);
                }
            }
            
        }
            break;
        case 3: //EFlux
        {
            
        }
            break;
            
        case 4: //EForce
        {
            sol.resize(nstate);
            VecDouble result(nstate);
            this->forceFunction(integrationpointdata.x, sol);
        }
            break;
            
        case 5: //ESolExact
        {
            sol.resize(nstate);
            VecDouble solint(nstate);
            Matrix dsol(nstate, nstate);
            this->SolutionExact(integrationpointdata.x, sol, dsol);
            
        }
            break;
        case 6: //EDSolExact
        {
            sol.resize(gradu.Rows() * gradu.Cols());
            VecDouble solint(nstate);
            Matrix dsol(nstate, nstate);
            this->SolutionExact(integrationpointdata.x, solint, dsol);
            
            for (int i = 0; i < gradu.Rows(); i++) {
                for (int j = 0; j < gradu.Cols(); j++) {
                    sol[i * gradu.Cols() + j] = dsol(i, j);
                }
            }
        }
            break;
            
            
        default:
        {
            std::cout << "Poisson::PostProcessSolution Wrong input var" << std::endl;
            DebugStop();
        }
    }
}

int Poisson::VariableIndex(const PostProcVar var) const {
    if (var == ENone) return ENone;
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;
    if (var == EFlux) return EFlux;
    if (var == EForce) return EForce;
    if (var == ESolExact) return ESolExact;
    if (var == EDSolExact) return EDSolExact;
    
    std::cout << "Poisson::VariableIndex : wrong input var" << std::endl;
    DebugStop();
    return -1;
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    if (!strcmp("Sol", name.c_str())) return ESol;
    if (!strcmp("DSol", name.c_str())) return EDSol;
    if (!strcmp("Flux", name.c_str())) return EFlux;
    if (!strcmp("Force", name.c_str())) return EForce;
    if (!strcmp("SolExact", name.c_str())) return ESolExact;
    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
    
    std::cout << "Poisson::VariableIndex : wrong input name" << std::endl;
    
    DebugStop();
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->Dimension();
    if (var == EFlux) return this->Dimension();
    if (var == EForce) return this->NState();
    if (var == ESolExact) return this->NState();
    if (var == EDSolExact) return this->Dimension();
    
    std::cout << "Poisson::NSolutionVariables : wrong input var" << std::endl;
    
    DebugStop();
    return -1;
}

double Poisson::Inner(Matrix &S,Matrix &T) const{
    std::cout << "Poisson::Inner : Not implemented yet" << std::endl;
    DebugStop();
    return -1;
};
