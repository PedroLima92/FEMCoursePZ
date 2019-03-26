//
//  CompElement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompElement.h"
#include "GeoElement.h"
#include "CompMesh.h"
#include "CompElementTemplate.h"
#include "DataTypes.h"
#include "MathStatement.h"
#include "PostProcessTemplate.h"
#include "PostProcess.h"

    CompElement::CompElement() : compmesh(0), index(-1), geoel(0), intrule(0), mat(0) {

    }

    CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo){
        index = ind;
        geoel = geo;
        compmesh = cmesh;
    }

    CompElement::CompElement(const CompElement &copy){
        index=copy.index;
        geoel=copy.geoel;
        intrule=copy.intrule;
        mat=copy.mat;
        compmesh=copy.compmesh;
    }

    CompElement &CompElement::operator=(const CompElement &copy){
        index=copy.index;
        geoel=copy.geoel;
        intrule=copy.intrule;
        mat=copy.mat;
        compmesh=copy.compmesh;
        return *this;
    }

    CompElement::~CompElement(){
        
    }

    CompElement *CompElement::Clone() const{
        
    }

    MathStatement *CompElement::GetStatement() const{
        return mat;
    }
    
    void CompElement::SetStatement(MathStatement *statement){
        mat=statement;
    }
    
    IntRule *CompElement::GetIntRule() const{
        return intrule;
    }
    
    void CompElement::SetIntRule(IntRule *intr){
        intrule = intr;
    }

    void CompElement::SetIndex(int64_t ind){
        index = ind;
    }

    GeoElement *CompElement::GetGeoElement() const{
        return geoel;
    }
    
    void CompElement::SetGeoElement(GeoElement *element){
        element=geoel;
    }
    
    CompMesh *CompElement::GetCompMesh() const{
        return compmesh;
    }
    
    void CompElement::SetCompMesh(CompMesh *mesh){
        compmesh = mesh;
    }
    
    void CompElement::InitializeIntPointData(IntPointData &data) const{
        //MathStatement *material = GetStatement();
        int dim = this->Dimension();
        int nshape = this->NShapeFunctions();
        int nstate = this->GetStatement()->NState();
        data.phi.resize(nshape,1);
        data.dphidksi.Resize(dim,nshape);
        data.dphidx.Resize(dim,nshape);
        data.axes.Resize(dim,3);
//        data.gradx.Resize
        data.x.resize(3);
        data.ksi.resize(dim);
        data.solution.resize(nstate);
        data.dsoldksi.Resize(dim,nstate);
        data.dsoldx.Resize(dim,nstate);
    }
    
    void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const{
        
        GeoElement *gel = this->GetGeoElement();
        Matrix gradx,Jac,JacInv;
        
        gel->X(intpoint, data.x);
        
        gel->GradX(intpoint, data.x, data.gradx);
        
        gel->Jacobian(data.gradx, Jac, data.axes, data.detjac, JacInv);
        this->ShapeFunctions(intpoint, data.phi, data.dphidksi);
        
        this->Convert2Axes(data.dphidksi, JacInv, data.dphidx);
     
    }

    void CompElement::Convert2Axes(const Matrix &dphi, const Matrix &jacinv, Matrix &dphidx) const{
        int nshape = dphi.Cols();
        int dim = dphi.Rows();
        dphidx.Resize(dim,nshape);
        int ieq;
        switch(dim){
            case 0:
            {
                
            }
                break;
            case 1:
            {
                dphidx = dphi;
                dphidx = dphidx*jacinv.GetVal(0,0);
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
                std::cout << "Error at " << __PRETTY_FUNCTION__ << " please implement the " << dim << "d Jacobian and inverse\n" <<std::endl;
            }
        }
        
    }


    void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const{
        MathStatement *material = GetStatement();
        if(!material){
            std::cout << " Material == NULL " << std::endl;
            return;
        }
       // this->InitializeElementMatrix(ek,ef);
        IntPointData data;
        this->InitializeIntPointData(data);
        double weight =0.;
        int nintrulepoints = GetIntRule()->NPoints();
        int dim = Dimension();
        VecDouble intpoint(dim,0.);
        ek.Zero();
        ef.Zero();
        int nshape = data.phi.size();
        int nstate = GetStatement()->NState();
        ek.Resize(nshape*nstate,nshape*nstate);
        ef.Resize(nshape*nstate,nstate);
        
        for (int intd_id = 0; intd_id < nintrulepoints; intd_id++) {
            intrule->Point(intd_id, data.ksi, data.weight);
            this->ComputeRequiredData(data, data.ksi);
            weight=data.weight;
            weight *=fabs(data.detjac);
            material->Contribute(data, weight, ek, ef);
        }
        
    }

    // Compute error and exact solution
    void CompElement::EvaluateError(std::function<void(const VecDouble &loc,VecDouble &val,Matrix &deriv)> fp,VecDouble &errors) const{
        
        MathStatement *mat=this->GetStatement();
        IntRule *intruleError = this->GetIntRule();
        int maxorder = 15;
        intruleError->SetOrder(maxorder);
        
        int nerrors = errors.size();
        std::fill(errors.begin(), errors.end(), 0.);
        
        int nstate = mat->NState();
        int dim = Dimension();
        VecDouble u_exact(nstate);
        Matrix du_exact(dim,nstate);
        VecDouble values(nerrors);
        
        IntPointData data;
        this->InitializeIntPointData(data);
        double weight =0.;
        int nintrulepoints = GetIntRule()->NPoints();
        VecDouble intpoint(dim,0.);
        
        for (int intd_id = 0; intd_id < nintrulepoints; intd_id++) {
            
            GetIntRule()->Point(intd_id, data.ksi, data.weight);
            this->ComputeRequiredData(data, data.ksi);
            GetMultiplyingCoeficients(data.coefs);
            data.ComputeSolution();
            weight=data.weight;
            weight *=fabs(data.detjac);
            
            if(fp) {
                fp(data.x,u_exact,du_exact);
            }
            mat->ContributeError(data,u_exact,du_exact,values);
            
            for(int ier = 0; ier < nerrors; ier++)
            {
                errors[ier] += values[ier]*weight;
            }
            
        }

        //Norma sobre o elemento
        for(int ier = 0; ier < nerrors; ier++){
            errors[ier] = sqrt(errors[ier]);
        }
        
//        intrule->SetOrder(prevorder);
        
        
    }

    
    void CompElement::Solution(VecDouble &intpoint, int var, VecDouble &sol) const{
    
        IntPointData data;
        this->InitializeIntPointData(data);
        GetMultiplyingCoeficients(data.coefs);
        MathStatement *mat = this->GetStatement();
        
        ComputeRequiredData(data,intpoint);
        data.ComputeSolution();
        sol.resize(3);
        
        mat->PostProcessSolution(data,var,sol);

        
    }



