//
//  CompElement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElement.h"
#include "CompElementTemplate.h"
#include "GeoElement.h"
#include "IntPointData.h"
#include "MathStatement.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "tpanic.h"
#include "GeoElement.h"

#include "Analysis.h"
#include <functional>

// Default constructor of CompElement
CompElement::CompElement() : compmesh(0), index(-1), geoel(0), intrule(0), mat(0){
    
}

// Constructor of CompElement
CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo){
    compmesh = cmesh;
    index = ind;
    geoel = geo;
    mat = 0;
    for (int64_t imaths=0; imaths<cmesh->GetMathVec().size(); imaths++) {
        if (cmesh->GetMath(imaths)->GetMatID() == geo->Material()) {
            this->SetStatement(cmesh->GetMath(imaths));
            break;
        }
    }
    compmesh->SetElement(index, this);
}

// Copy constructor of CompElement
CompElement::CompElement(const CompElement &copy){
    index = copy.index;
    geoel = copy.geoel;
    compmesh = copy.compmesh;
    intrule = copy.intrule;
    mat = copy.mat;
    compmesh->SetElement(index, this);
}
    
// Operator of copy
CompElement &CompElement::operator=(const CompElement &copy){
    index = copy.index;
    geoel = copy.geoel;
    compmesh = copy.compmesh;
    intrule = copy.intrule;
    mat = copy.mat;
    compmesh->SetElement(index, this);
    return *this;
}
    
// Destructor of CompElement
CompElement::~CompElement(){
    
}

CompElement *CompElement::Clone() const{
//    return new CompElement(*this);
    DebugStop();
}

// Return the material object associated with the element
MathStatement *CompElement::GetStatement() const{
    return mat;
}
    
// Set the material object associated with the element
void CompElement::SetStatement(MathStatement *statement){
    mat = statement;
}
    
// Return integration rule established
IntRule *CompElement::GetIntRule() const{
    return intrule;
}
    
// Set integration rule established
void CompElement::SetIntRule(IntRule *integ){
    intrule = integ;
}
    
// Set element index
void CompElement::SetIndex(int64_t ind){
    index = ind;
}
    
// Return the geometric element associated
GeoElement *CompElement::GetGeoElement() const{
    return geoel;
}
    
// Set the geometric element associated
void CompElement::SetGeoElement(GeoElement *element){
    geoel = element;
}
    
// Return a pointer to the element computational mesh
CompMesh *CompElement::GetCompMesh() const{
    return compmesh;
}
    
// Set a pointer to the element computational mesh
void CompElement::SetCompMesh(CompMesh *mesh){
    compmesh = mesh;
}
    
// Initialize integration points data object
void CompElement::InitializeIntPointData(IntPointData &data) const{
    
    const int nstate = mat->NState();
    const int nshape = this->NShapeFunctions();
    const int dim = Dimension();
    
    data.weight = 0;
    data.detjac = 0;
    
    data.ksi.resize(3, 0.);
    data.x.resize(3, 0.);
    data.phi.resize(nshape, 0.);
    
    data.dphidksi.Resize(dim, nshape);
    data.dphidksi.Zero();
    
    data.gradx.Resize(3,3);
    data.gradx.Zero();
    
    data.axes.Resize(dim, 3);
    data.axes.Zero();
    
    data.dphidx.Resize(dim, nshape);
    data.dphidx.Zero();
    
    data.solution.resize(nshape*nstate, 0.);
    
    data.dsoldksi.Resize(dim, nshape*nstate);
    data.dsoldksi.Zero();
    
    data.dsoldx.Resize(dim, nshape*nstate);
    data.dsoldx.Zero();
}

// Compute and fill integration points data object
void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const{
    
    if (!geoel){
        DebugStop();
        std::cout << "CompElement::ComputeRequiredData : Null GeoElement associated" << std::endl;
    }
    
    geoel->X(intpoint, data.x);
    geoel->GradX(intpoint, data.x, data.gradx);

    int ncols = data.gradx.Cols();
    int dim   = ncols;
    
    Matrix jac(dim,dim);
    Matrix jacinv(dim,dim);

    geoel->Jacobian(data.gradx, jac, data.axes, data.detjac, jacinv);
    this->ShapeFunctions(intpoint, data.phi, data.dphidksi);
    this->Convert2Axes(data.dphidksi, jacinv, data.dphidx);
    
}

void CompElement::Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const
{
    int nshape = this->NShapeFunctions();
    int dim = this->Dimension();
    int ieq;
    switch(dim){
        case 0:
        {
            
        }
            break;
        case 1:
        {
            for(ieq = 0; ieq < nshape; ieq++) {
                dphidx(0,ieq) = dphi.GetVal(0,ieq);
                dphidx(0,ieq) *= jacinv.GetVal(0,0);
            }
        }
            break;
        case 2:
        {
            for(ieq = 0; ieq < nshape; ieq++) {
                dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq);
                dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq);
            }
        }
            break;
        case 3:
        {
            for(ieq = 0; ieq < nshape; ieq++) {
                dphidx(0,ieq) = jacinv.GetVal(0,0)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,0)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,0)*dphi.GetVal(2,ieq);
                dphidx(1,ieq) = jacinv.GetVal(0,1)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,1)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,1)*dphi.GetVal(2,ieq);
                dphidx(2,ieq) = jacinv.GetVal(0,2)*dphi.GetVal(0,ieq) + jacinv.GetVal(1,2)*dphi.GetVal(1,ieq) + jacinv.GetVal(2,2)*dphi.GetVal(2,ieq);
            }
        }
            break;
        default:
        {
            std::cout << "CompElement::Convert2Axes : Not implemented Jacobian" << std::endl;
            DebugStop();
        }
    } //switch
    
}

void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
    
    IntRule * intrule = this->GetIntRule();
    MathStatement * material = this->GetStatement();
    
    if(!material){
        std::cout << "CompElement:: Null MathStatement" << std::endl;
        DebugStop();
    }
    
    IntPointData data;
    this->InitializeIntPointData(data);
    double weight = 0.;
    
    for(int int_ind = 0; int_ind < intrule->NPoints(); int_ind++){
        intrule->Point(int_ind, data.ksi, data.weight);
        weight = data.weight;
        this->ComputeRequiredData(data, data.ksi);
        weight *= fabs(data.detjac);
        material->Contribute(data, weight, ek, ef);
    }
    
}

void CompElement::EvaluateError(std::function<void(const VecDouble &loc,VecDouble &val,Matrix &deriv)> fp,
                                VecDouble &errors) const{
    
    MathStatement * material = this->GetStatement();
    if (!material) {
        std::cout << "CompElement::EvaluateError : no material for this element\n" << std::endl;
        DebugStop();
    }
    
    int NErrors = material->NEvalErrors();
    
    int dim = Dimension();
    if(compmesh->GetGeoMesh()->Dimension() < dim) return;
    
    // Adjust the order of the integration rule
    IntRule * intrule = this->GetIntRule();
    intrule->SetOrder(14);
    
    int nstate = material->NState();
    int var = 0;
    
    VecDouble u_exact(nstate, 0.);
    Matrix du_exact(dim, nstate);
    du_exact.Zero();
    VecDouble intpoint(dim, 0.), values(NErrors, 0.);
    double weight = 0.;
    
    IntPointData data;
    this->InitializeIntPointData(data);
    int nintpoints = intrule->NPoints();
    
    for(int nint = 0; nint < nintpoints; nint++) {
        
        std::fill(values.begin(), values.end(), 0.);
        intrule->Point(nint,intpoint,weight);
        
        this->ComputeRequiredData(data, intpoint);
        weight *= fabs(data.detjac);
        
        if(fp) {
            fp(data.x,u_exact,du_exact); // solução analítica
            this->Solution(intpoint, var, data.solution, data.dsoldx); // solução interpolada
            
            material->ContributeError(data, u_exact, du_exact, values);
            
            for(int ier = 0; ier < NErrors; ier++){
                errors[ier] += weight*values[ier];
            }
        }
        
    }//fim for : integration rule
    
    //Norma sobre o elemento
    for(int ier = 0; ier < NErrors; ier++){
        errors[ier] = sqrt(errors[ier]);
    }//for ier
    
}

// Compute the solution and its gradient at a parametric point
// for dsol the row indicates the direction, the column indicates the state variable
void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol, TMatrix &dsol) const{
    
    IntPointData data;
    
    this->InitializeIntPointData(data);
    this->ComputeRequiredData(data, intpoint);
    
    GetMultiplyingCoeficients(data.coefs);
    data.ComputeSolution();
    
    sol=data.solution;
    dsol=data.dsoldx;
    
}

