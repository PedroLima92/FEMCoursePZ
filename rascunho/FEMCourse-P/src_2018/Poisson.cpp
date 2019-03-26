//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "Poisson.h"
#include "MathStatement.h"
#include "tpanic.h"

    Poisson::Poisson(){
        
    }
    
    Poisson::Poisson(int materialid, Matrix &perm){
        SetMatID(materialid);
        permeability=perm;
    }
    
    Poisson::Poisson(const Poisson &copy){
        permeability=copy.permeability;
        forceFunction=copy.forceFunction;
    }
    
    Poisson &Poisson::operator=(const Poisson &copy){
        permeability=copy.permeability;
        forceFunction=copy.forceFunction;
        return *this;
    }
    
    Poisson *Poisson::Clone() const{
       // return new MathStatement(*this);
    }

    Poisson::~Poisson(){
        
    }
    
    Matrix Poisson::GetPermeability() const{
        return permeability;
    }
    
    void Poisson::SetPermeability(const Matrix &perm){
        permeability=perm;
    }

    int Poisson::NEvalErrors() const{
        return 3;
    }

    int Poisson::VariableIndex(const PostProcVar var) const{
        
        int nvar = 10;
        for (int i=0; i<nvar; i++) {
            if (var==PostProcVar(i)) {
                return i;
            }
        }
        return 0;
    }

    // Return the variable index associated with the name
    Poisson::PostProcVar Poisson::VariableIndex(const std::string &name){
  
        if (!strcmp("Sol", name.c_str()))  return ESol;
        if (!strcmp("Solution", name.c_str()))  return ESol;
        if (!strcmp("DSol", name.c_str()))  return EDSol;
        if (!strcmp("DSolution", name.c_str()))  return EDSol;
        if (!strcmp("Flux", name.c_str()))         return EFlux;
        if (!strcmp("Force", name.c_str()))   return EForce;
        if (!strcmp("Sol_exact", name.c_str()))   return ESolExact;
        if (!strcmp("DSol_exact", name.c_str()))   return EDSolExact;
        
        std::cout  << " Var index not implemented " << std::endl;
        DebugStop();
        return ENone;
    }

    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    int Poisson::NSolutionVariables(const PostProcVar var){
        
        switch(var) {
                
            case ESol:
                return this->NState(); // Solution, Vector
            case EDSol:
                return this->Dimension();// Derivative of solution, Vector
            case EFlux:
                return this->Dimension(); // Flux, Vector
            case EForce:
                return this->NState(); // Force vector, Vector
            case ESolExact:
                return this->NState(); // Sol_exact, Vector
            case EDSolExact:
                return this->Dimension(); // DSol_exact, Vector

            default:
            {
                std::cout  << " Var index not implemented " << std::endl;
                DebugStop();
            }
        }
        return 0;
        
    }

    void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const{
        
        VecDouble &phi=data.phi;
        Matrix &dphi=data.dphidx;
        VecDouble &x = data.x;
        Matrix &axes =data.axes;
        
       // dphi.Print();
      //  data.axes.Print();
        
        
        Matrix perm = GetPermeability();
//        Matrix Kdphi(2,2,0.);
        
  //      dphi.Print();
        
        Matrix dphiU(Dimension(),dphi.Cols());
        Axes2XYZ(dphi, dphiU, data.axes);
  //      data.axes.Print();
  //      dphiU.Print();
        
        // shape index for nstate
        int dim = Dimension();
        int nphi= phi.size();
        int nshape = nphi*NState();
        int index=0, inormal = 0;
        VecDouble shapeindex(nshape,0.), normalindex(nshape,0.);
        for (int i = 0; i<nphi; i++) {
            inormal = 0;
            for (int s = 0; s<NState(); s++) {
                shapeindex[index]=i;
                normalindex[index]=inormal;
                index++;
                inormal++;
            }
        }
        
        Matrix Normalvec(NState(),NState(),0.);
        for (int in =0; in<NState(); in++) {
            Normalvec(in,in)=1.;
        }
        
        for(int i = 0; i < nshape; i++ )
        {
            int iphi = shapeindex[i];
            int ivec = normalindex[i];
            Matrix phiVi(NState(),1,0.),GradVi(NState(),dim,0.);
            for (int e=0; e<NState(); e++) {
                phiVi(e,0) = phi[iphi]*Normalvec(e,ivec);
            
                // Grad V
                for (int f=0; f<dim; f++) {
                    GradVi(e,f) = Normalvec(e,ivec)*dphiU(f,iphi);
                }
            }
            
            // Force vector :
            
            VecDouble f(NState(),0.);
            forceFunction(x,f);

            double phi_dot_f = 0.0;
            for (int e=0; e<NState(); e++) {
                phi_dot_f += phiVi(e,0)*f[e];
            }
            
            EF(i,0) += phi_dot_f * weight;
            
            
            for(int j = 0; j < nshape; j++){
                int jphi = shapeindex[j];
                int jvec = normalindex[j];
                
                Matrix GradVj(NState(),dim,0.), KGradVj(NState(),dim,0.);
                for (int e=0; e<NState(); e++) {
                    for (int f=0; f<dim; f++) {
                        GradVj(e,f) = Normalvec(e,jvec)*dphiU(f,jphi);
                    }
                }
                
//                if (NState()==dim) {
//                    // K * Grad U
//                    for (int ik=0; ik<dim; ik++) {
//                        for (int jk=0; jk<NState(); jk++) {
//                            for (int l=0; l<dim; l++){
//                                KGradVj(ik,jk) += perm(ik,l)*GradVj(l,jk);
//                            }
//                        }
//                    }
//
//                    double val = Inner(GradVi, KGradVj);
//                    EK(i,j) += weight * val;
//                }else{
                
                    GradVj.Transpose();
                    for (int ik=0; ik<NState(); ik++) {
                        for (int jk=0; jk<dim; jk++) {
                            for (int l=0; l<dim; l++){
                                KGradVj(ik,jk) += perm(jk,l)*GradVj(l,ik);
                            }
                        }
                    }
                    
                    double val = Inner(GradVi, KGradVj);
                    EK(i,j) += weight * val;
                    
//                }
                
                
                
            }
            
        }
        

//        EK.Print();
//        EF.Print();
//        std::cout<<std::endl;

    }

    // Method to implement error over element's volume
    void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const{
        
        du_exact.Resize(Dimension(), NState());
        u_exact.resize(NState());
        errors.resize(NEvalErrors());
        std::fill(errors.begin(), errors.end(),0.);
        VecDouble Sol, DSol;
       
        this->PostProcessSolution(data,ESol,Sol);
      //  DSol = this->PostProcessSolution(data,EDSol);
        
        Matrix &dsol = data.dsoldx;
        Matrix dsolxy;
        Axes2XYZ(dsol, dsolxy, data.axes);
        
        
        //values[0] : erro norma L2
        double diff;
        errors[0] = 0.;
        for(int i=0; i<NState(); i++) {
            diff = Sol[i] - u_exact[i];
            errors[0]  += diff*diff;
        }
        
        //values[1] : erro em semi norma H1
        errors[1] = 0.;
        double diff2;
        for(int i=0; i<Dimension(); i++) {
            for(int j=0; j<NState(); j++) {
      //          diff2 = DSol[j+i*NState()] - du_exact(i,j);
                diff2 = dsolxy(i,j)-du_exact(i,j);
                errors[1]  += diff2*diff2;
            }
        }
        
        //values[2] : erro em norma H1 <=> norma Energia
        errors[2]  = errors[1]+errors[2];
        
        
    }

    void Poisson::PostProcessSolution(const IntPointData &data, const int varindex, VecDouble &Solout) const{
        
        PostProcVar var = PostProcVar(varindex);
        
        VecDouble u_h = data.solution;
        int usize = u_h.size();
        Matrix du_h = data.dsoldx;
        int duRows = du_h.Rows();
        int duCols = du_h.Cols();
        
        Matrix perm = GetPermeability();
        int permRows = perm.Rows();
        int permCols = perm.Cols();


        switch(var) {
                
            case ESol: //ux, uy, uz
            {
                Solout.resize(NState());
                for (int isol = 0; isol< NState(); isol++) {
                    Solout[isol]=u_h[isol];
                }
                
            }
                break;
                
            case EDSol: //Grad u
            {
                Solout.resize(duRows*duCols,0.);
                for (int i = 0; i<duRows ; i++) {
                    for (int j = 0; j<duCols; j++) {
                        Solout[j+i*duCols]=du_h(i,j);
                    }
                }
            }
                break;
                
            case EFlux: //Flux = Perm x Grad u
            {
                Solout.resize(permRows*duCols,0.);
                Matrix flux(permRows,duCols,0.);
                
                for (int ik = 0; ik<permRows; ik++) {
                    for (int jk = 0; jk<duCols; jk++) {
                        for (int lk = 0; lk<permCols; lk++) {
                            flux(ik,jk)+=perm(ik,lk)*du_h(lk,jk);
                        }
                    }
                }
                
                Solout.resize(permRows*duCols,0.);
                for (int i = 0; i<duRows ; i++) {
                    for (int j = 0; j<duCols; j++) {
                        Solout[j+i*duCols]=flux(i,j);
                    }
                }
            }
                break;
                
            case EForce: //f
            {
                Solout.resize(NState());
                VecDouble f(NState(),0.0);
                forceFunction(data.x,f);
                
                for (int isol = 0; isol< NState(); isol++) {
                    Solout[isol]=f[isol];
                }
                
                
            }
                break;
                
            case ESolExact: //u_exact
            {
                Solout.resize(NState());
                VecDouble sol(NState(),0.0);
                Matrix dsol(NState(),1,0.0);
                if(SolutionExact){
                    SolutionExact(data.x,sol,dsol);
                }
                for (int isol = 0; isol< NState(); isol++) {
                    Solout[isol]=sol[isol];
                }
                
//                std::cout << std::endl;
//                for (int isol = 0; isol< NState(); isol++) {
//                    Solout[isol]=sol[isol];
//                    std::cout << "x[i] = " << data.x[isol]<<std::endl;
//                    std::cout << "v[i] = " << sol[isol]<<std::endl;
//                }
//                std::cout << std::endl;
                
            }
                break;
//
//            case ESolExact: //du_exact
//            {
//                TPZVec<STATE> sol(3,0.0);
//                if(this->HasForcingFunctionExact()){
//                    this->fForcingFunctionExact->Execute(datavec[pindex].x, sol, gradu); // @omar::check it!
//                }
//                Solout[0] = sol[2]; // px
//
//            }
//                break;
                
                
            default:
            {
                std::cout  << " Var index not implemented " << std::endl;
                DebugStop();
            }
        }
        
        
    }

double Poisson::Inner(Matrix &S, Matrix &T) const{
    
    double Val = 0.;
    
    for(int i = 0; i < S.Rows(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



//    // Prepare and print post processing data
//    void Poisson::EvaluateSolution(const IntPointData &integrationpointdata, PostProcess &defPostProc) const{
//
//    }

