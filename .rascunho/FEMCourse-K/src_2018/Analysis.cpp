
//
//  Analysis.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#include "Analysis.h"
#include "Assemble.h"
#include "MathStatement.h"
#include "CompElement.h"
#include "CompElementTemplate.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "VTKGeoMesh.h"
#define ARMA_USE_SUPERLU
#include <armadillo>
#include <iostream>

using namespace arma;

Analysis::Analysis() : cmesh(0), Solution(), GlobalSystem(), RightHandSide(){
    
};

Analysis::Analysis(const Analysis &cp){
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
};

Analysis &Analysis::operator=(const Analysis &cp){
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
    return *this;
};

Analysis::Analysis(CompMesh *compmesh){
    cmesh = compmesh;
};

Analysis::~Analysis(){
    
};

void Analysis::SetMesh(CompMesh *compmesh){
    cmesh = compmesh;
};

CompMesh *Analysis::Mesh() const{
    return cmesh;
};

void Analysis::RunSimulation(){
    
    Assemble assem(cmesh);
    int neq = assem.NEquations();
    
//    Matrix globmat(neq,neq,0.);
//    Matrix rhs(neq,1,0.);
    sp_mat globmat(neq,neq);
    mat rhs(neq,1);
    
    assem.Compute(globmat, rhs);
    
//    globmat.Solve_LU(rhs);
    mat solarma = spsolve(globmat, rhs);
    
    VecDouble sol(neq,0.);
    for (int64_t irows=0; irows<neq; irows++) {
        sol[irows] = solarma(irows,0);
//        sol[irows] = rhs(irows,0);
    }
    
    cmesh->LoadSolution(sol);
    
};

void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const{
    VTKGeoMesh::PrintSolVTK(cmesh, defPostProc, filename);
};

VecDouble Analysis::PostProcessError(std::ostream &out, PostProcess &defPostProc) const{
    
    VecDouble sol = cmesh->Solution();
    int64_t neq = sol.size();
    VecDouble ux(neq,0.);
    VecDouble errors(3,0.);
    VecDouble values(3,0.);
    
    int64_t i, nel = cmesh->GetElementVec().size();
    
    for(i=0; i<nel; i++) {
        CompElement *cel = cmesh->GetElement(i);
        
        if(cel) {
            MathStatement *mat = cel->GetStatement();
            if(mat->GetMatID() > 0)
            {
                std::fill(errors.begin(), errors.end(), 0.);
                cel->EvaluateError(defPostProc.GetExact(), errors);
                for(int ier = 0; ier < 3; ier++)
                {
                    values[ier] += errors[ier] * errors[ier];
                }
            }
        }
    }
    
    int nerrors = errors.size();
    VecDouble ervec(nerrors, -10.);
    
    if (nerrors < 3) {
        std::cout << "Analysis::PostProcessError - At least 3 norms are expected" << std::endl;
        DebugStop();
    }
    else{
        out << "############" << std::endl;
        out <<"Norma H1 or L2 -> p = "  << sqrt(values[0]) << std::endl;
        out <<"Norma L2 or L2 -> u = "    << sqrt(values[1]) << std::endl;
        out << "Semi-norma H1 or L2 -> div = "    << sqrt(values[2])  << std::endl;
        for(int ier = 3; ier < nerrors; ier++)
            out << "other norms = " << sqrt(values[ier]) << std::endl;
    }
    // Returns the calculated errors.
    for(i=0;i<nerrors;i++){
        ervec[i] = sqrt(values[i]);
    }
    return ervec;
};
