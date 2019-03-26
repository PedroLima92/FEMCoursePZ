//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//
#include "Assemble.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "DOF.h"
#include "CompMesh.h"

#include "Analysis.h"
#include "MathStatement.h"

using namespace arma;
    
Assemble::Assemble(): cmesh(0){
};

Assemble::Assemble(CompMesh *mesh){
    cmesh = mesh;
};

Assemble::Assemble(const Assemble &copy){
    cmesh = copy.cmesh;
};

Assemble &Assemble::operator=(const Assemble &copy){
    cmesh = copy.cmesh;
    return *this;
};

void Assemble::SetMesh(CompMesh *mesh){
    cmesh = mesh;
};

/// Compute the total number of equations
int64_t Assemble::NEquations(){
    int64_t neq = 0;
    std::vector<DOF> dofs;
    dofs = cmesh->GetDOFVec();
    int64_t ndofs = dofs.size();
    
    for (int64_t idofs=0; idofs< ndofs; idofs++) {
        neq += dofs[idofs].GetNShape()*dofs[idofs].GetNState();
    }
    
    return neq;
};

/// Optimize the bandwidth of the global system of equations
void Assemble::OptimizeBandwidth(){
    std::cout << "Assemble::OptimizeBandwidth : Not implemented yet" << std::endl;
    DebugStop();
};

/// Compute the global stiffness matrix and right hand side
//void Assemble::Compute(sp_mat &globmat, mat &rhs){
void Assemble::Compute(sp_mat &globmat, mat &rhs){
    int nelem = cmesh->GetElementVec().size();
    
    for (int iel=0; iel<nelem; iel++) {
        CompElement * cel = cmesh->GetElement(iel);
        
        if (!cel) {
            continue;
        }
        
        Matrix ek, ef;
        int nshape = cel->NShapeFunctions();
        int nstate = cel->GetStatement()->NState();
        ek.Resize(nshape*nstate, nshape*nstate);
        ef.Resize(nshape*nstate, 1);
        ek.Zero();
        ef.Zero();
        cel->CalcStiff(ek, ef);
        
        int nsides = cel->GetGeoElement()->NSides();
        VecInt DOFfirsteq(0);
        
        for (int64_t idofs=0; idofs<nsides; idofs++) {
            int64_t dofindex = cel->GetDOFIndex(idofs);
            DOF dof = cmesh->GetDOF(dofindex);
            if (dof.GetNShape() > 0) {
                DOFfirsteq.push_back(dof.GetFirstEquation());
            }
        }
        
        for (int64_t idofs=0; idofs<DOFfirsteq.size(); idofs++) {
            int64_t dofindexi = cel->GetDOFIndex(idofs);
            int64_t fei = DOFfirsteq[idofs];
            DOF dofi = cmesh->GetDOF(dofindexi);
            
            for (int64_t jdofs=0; jdofs<DOFfirsteq.size(); jdofs++) {
                int64_t fej = DOFfirsteq[jdofs];
                int64_t dofindexj = cel->GetDOFIndex(jdofs);
                DOF dofj = cmesh->GetDOF(dofindexj);
                
                for (int64_t idim=0; idim<dofi.GetNState()*dofi.GetNShape(); idim++) {
                    for (int64_t jdim=0; jdim<dofj.GetNState()*dofj.GetNShape(); jdim++) {
                        globmat(fei+idim, fej+jdim) += ek(idofs*dofj.GetNState()+idim,jdofs*dofj.GetNState()+jdim);
                    }
                }
            }
            
            for (int64_t idim = 0; idim<dofi.GetNState()*dofi.GetNShape(); idim++) {
                rhs(fei+idim,0) += ef(idofs*dofi.GetNState()+idim,0);
            }
            
        }
    }
    
    std::cout << "Assemble: Compute done!" << std::endl;
};
