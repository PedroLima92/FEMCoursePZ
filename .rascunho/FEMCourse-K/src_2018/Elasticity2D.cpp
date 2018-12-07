//
//  Header.h
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 7/5/18.
//
//

#include "Elasticity2D.h"
    
// Default constructor of Poisson
Elasticity2D::Elasticity2D(){
    fPlaneStress = 1;
};

// Constructor of Poisson
Elasticity2D::Elasticity2D(int materialid, double elast, double nu){
    fE = elast;
    fnu = nu;
    SetMatID(materialid);
    fPlaneStress = 1;
    
};

// Copy constructor of Poisson
Elasticity2D::Elasticity2D(const Elasticity2D &copy){
    fnu = copy.fnu;
    fE = copy.fE;
    ff = copy.ff;
    flambda = copy.flambda;
    fmu = copy.fmu;;
    fPlaneStress = copy.fPlaneStress;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
}

// Operator of copy
Elasticity2D &Elasticity2D::operator=(const Elasticity2D &copy){
    fnu = copy.fnu;
    fE = copy.fE;
    ff = copy.ff;
    flambda = copy.flambda;
    fmu = copy.fmu;;
    fPlaneStress = copy.fPlaneStress;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    return *this;
}

// Method for creating a copy of the element
Elasticity2D *Elasticity2D::Clone() const{
    
}

// Destructor of Poisson
Elasticity2D::~Elasticity2D(){
    
};

// Return the permeability matrix
double Elasticity2D::GetYoungModulus() const{
    return fE;
};

// Set the permeability matrix
void Elasticity2D::SetYoungModulus(const double &elast){
    fE = elast;
};

// Return the permeability matrix
double Elasticity2D::GetPoisson() const{
    return fnu;
};

// Set the permeability matrix
void Elasticity2D::SetPoisson(const double &nu){
    fnu= nu;
};

int Elasticity2D::NEvalErrors() const{
    return 3;
};

// Return the number of state variables
int Elasticity2D::NState() const{
    return 2;
};

int Elasticity2D::VariableIndex(const PostProcVar var) const{
    DebugStop();
};

// Return the variable index associated with the name
Elasticity2D::PostProcVar Elasticity2D::VariableIndex(const std::string &name){
    DebugStop();
};

// Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
int Elasticity2D::NSolutionVariables(const PostProcVar var){
    DebugStop();
};

void Elasticity2D::SetPlaneStressState(int state){
    fPlaneStress = state;
}

// Method to implement integral over element's volume
void Elasticity2D::Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const{
    
    VecDouble &phi = integrationpointdata.phi;
    Matrix &dphi = integrationpointdata.dphidx;
    VecDouble  &x = integrationpointdata.x;
    Matrix  &axes = integrationpointdata.axes;
    
    double LambdaL     = (fE*fnu)/((1+fnu)*(1-2*fnu));
    double MuL         = fE/(2*(1+fnu));
    
    int nphi = phi.size();
    int dim = dphi.Rows();
    
    VecDouble du(dim,0.);
    VecDouble dv(dim,0.);
    
    VecDouble P(3,0.);
    if (forceFunction) {
        forceFunction(x,P);
    }
    
    for( int iphi = 0; iphi < nphi; iphi++ ) {
        
        du[0] = (dphi(0,iphi)*axes(0,0)+dphi(1,iphi)*axes(1,0)); //du/dx
        du[1] = (dphi(0,iphi)*axes(0,1)+dphi(1,iphi)*axes(1,1)); //du/dy
        
        EF(2*iphi,0)     +=    weight*P[0]*phi[iphi];    // direcao x
        EF(2*iphi+1,0)   +=    weight*P[1]*phi[iphi];    // direcao y
        
        for( int jphi = 0; jphi < nphi; jphi++ ) {
            
            dv[0] = (dphi(0,jphi)*axes(0,0)+dphi(1,jphi)*axes(1,0)); // dv/dx
            dv[1] = (dphi(0,jphi)*axes(0,1)+dphi(1,jphi)*axes(1,1)); // dv/dy
            
            if (this->fPlaneStress == 1)
            {
                /* Plain stress state */
                EK(2*iphi, 2*jphi)	     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*dv[0]*du[0]		+ (MuL)*dv[1]*du[1]);
                
                EK(2*iphi, 2*jphi+1)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*dv[0]*du[1]			+ (MuL)*dv[1]*du[0]);
                
                EK(2*iphi+1, 2*jphi)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*dv[1]*du[0]			+ (MuL)*dv[0]*du[1]);
                
                EK(2*iphi+1, 2*jphi+1)     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*dv[1]*du[1]		+ (MuL)*dv[0]*du[0]);
            }
            else
            {
                /* Plain Strain State */
                EK(2*iphi,2*jphi)         += weight*	((LambdaL + 2*MuL)*dv[0]*du[0]	+ (MuL)*dv[1]*du[1]);
                
                EK(2*iphi,2*jphi+1)       += weight*	(LambdaL*dv[0]*du[1]			+ (MuL)*dv[1]*du[0]);
                
                EK(2*iphi+1,2*jphi)       += weight*	(LambdaL*dv[1]*du[0]			+ (MuL)*dv[0]*du[1]);
                
                EK(2*iphi+1,2*jphi+1)     += weight*	((LambdaL + 2*MuL)*dv[1]*du[1]	+ (MuL)*dv[0]*du[0]);
                
            }
        }
    }
    
};

// Method to implement error over element's volume
void Elasticity2D::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
    
    errors.resize(NEvalErrors(), 0.);
    
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
    
    
//    Matrix dudx(2,2), stress(2,2), stressexact(2,2);
//    Axes2XYZ(integrationpointdata.dsoldx, dudx, integrationpointdata.axes);
//    
//    dudx = integrationpointdata.dsoldx;
//    if (fPlaneStress == 1)
//    {
//        stress(0,0) = (4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*dudx.g(0,0)+(2*(fmu)*(flambda)/(flambda+2*fmu))*dudx.g(1,1);
//        stress(0,1) = (fmu)*(dudx.g(1,0)+dudx.g(0,1));
//        stress(1,0) = stress(0,1);
//        stress(1,1) = (2*(fmu)*(flambda)/(flambda+2*fmu))*dudx.g(0,0)+(4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*dudx.g(1,1);
//        /* Plain stress state */
//    }
//    else
//    {
//        stress(0,0) = (flambda + 2*fmu)*dudx.g(0,0)+flambda*dudx.g(1,1);
//        stress(1,0) = fmu*(dudx.g(1,0)+dudx.g(0,1));
//        stress(0,1) = stress(1,0);
//        stress(1,1) = (flambda + 2*fmu)*dudx.g(1,1)+flambda*dudx.g(0,0);
//        /* Plain Strain State */
//    }
//    
//    if (fPlaneStress == 1)
//    {
//        stressexact(0,0) = (4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*du_exact.g(0,0)+(2*(fmu)*(flambda)/(flambda+2*fmu))*du_exact.g(1,1);
//        stressexact(0,1) = (fmu)*(du_exact.g(1,0)+du_exact.g(0,1));
//        stressexact(1,0) = stressexact(0,1);
//        stressexact(1,1) = (2*(fmu)*(flambda)/(flambda+2*fmu))*du_exact.g(0,0)+(4*(fmu)*(flambda+fmu)/(flambda+2*fmu))*du_exact.g(1,1);
//        /* Plain stress state */
//    }
//    else
//    {
//        stressexact(0,0) = (flambda + 2*fmu)*du_exact.g(0,0)+flambda*du_exact.g(1,1);
//        stressexact(1,0) = fmu*(du_exact.g(1,0)+du_exact.g(0,1));
//        stressexact(0,1) = stressexact(1,0);
//        stressexact(1,1) = (flambda + 2*fmu)*du_exact.g(1,1)+flambda*du_exact.g(0,0);
//        /* Plain Strain State */
//    }
//    
//    VecDouble uexact = u_exact;
//    VecDouble sol = integrationpointdata.solution;
//    Matrix duexact = du_exact;
//    
//    double L2 = 0.;
//    L2 = (sol[0]-uexact[0])*(sol[0]-uexact[0])+(sol[1]-uexact[1])*(sol[1]-uexact[1]);
//    
//    double H1 = 0.;
//    double energy = 0.;
//    for (int i=0; i<2; i++) {
//        for (int j=0; j<2; j++) {
//            H1 += (dudx(i,j)-duexact(i,j))*(dudx(i,j)-duexact(i,j));
//            energy += (stress(i,j)-stressexact(i,j))*(dudx(i,j)-duexact(i,j));
//        }
//    }
//    errors[0] = energy;
//    errors[1] = L2;
//    errors[2] = H1;
};

// Prepare and print post processing data
void Elasticity2D::PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const{
    DebugStop();
};

double Elasticity2D::Inner(Matrix &S,Matrix &T) const{
    DebugStop();
};
