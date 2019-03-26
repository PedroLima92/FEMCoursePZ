//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompMesh.h"
#include "GeoMesh.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "MathStatement.h"
#include "DOF.h"
#include <map>

    CompMesh::CompMesh():geomesh(0),compelements(0),dofs(0),mathstatements(0){
        
    }

    CompMesh::CompMesh(GeoMesh *gmesh) : geomesh(gmesh){
        
    }

    CompMesh::CompMesh(const CompMesh &copy){
        compelements = copy.compelements;
        dofs = copy.dofs;
        mathstatements = copy.mathstatements;
    }
    
    CompMesh::~CompMesh(){
        
    }

    GeoMesh *CompMesh::GetGeoMesh() const{
        return geomesh;
    }

    void CompMesh::SetGeoMesh(GeoMesh *gmesh){
        geomesh=gmesh;
    }

    void CompMesh::SetNumberElement(int64_t nelem){
        compelements.resize(nelem);
    }
    
    void CompMesh::SetNumberDOF(int64_t ndof){
        dofs.resize(ndof);
    }
    
    void CompMesh::SetNumberMath(int nmath){
        mathstatements.resize(nmath);
    }
    
    void CompMesh::SetElement(int64_t elindex, CompElement *cel){
        compelements[elindex]=cel;
    }
    
    void CompMesh::SetDOF(int64_t index, const DOF &dof){
        dofs[index]=dof;
    }
    
    void CompMesh::SetMathStatement(int index, MathStatement *math){
        mathstatements[index]=math;
    }
    
    DOF &CompMesh::GetDOF(int64_t dofindex){
        return dofs[dofindex];
    }
    
    CompElement *CompMesh::GetElement(int64_t elindex) const{
        return compelements[elindex];
    }
    
    MathStatement *CompMesh::GetMath(int matindex) const{
        return mathstatements[matindex];
    }

    std::vector<CompElement *> CompMesh::GetElementVec() const{
        return compelements;
    }
    
    std::vector<DOF> CompMesh::GetDOFVec() const{
        return dofs;
    }
    
    std::vector<MathStatement *> CompMesh::GetMathVec() const{
        return mathstatements;
    }
    
    void CompMesh::SetElementVec(const std::vector<CompElement *> &vec){
        compelements=vec;
    }
    
    void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec){
        dofs=dofvec;
    }
    
    void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec){
        mathstatements=mathvec;
    }

    void CompMesh::AutoBuild(){
        int nel = GetGeoMesh()->NumElements();
        SetNumberElement(nel);
        for (int iel = 0; iel< nel; iel++) {
            GeoElement *gel = GetGeoMesh()->Element(iel);
            CompElement *cel = gel->CreateCompEl(this, iel);
            SetElement(iel, cel);
            
        }
        
        this->Resequence();
    }

// Initialize the datastructure FirstEquation of the DOF objects
    void CompMesh::Resequence(){
        
        int64_t ncon = GetNumberDOF();
        int64_t firstEq = 0;
        
        for (int idof=0; idof<ncon; idof++) {
            this->GetDOF(idof).SetFirstEquation(firstEq);
            int dofsize = this->GetDOF(idof).GetNShape()*this->GetDOF(idof).GetNState();
            firstEq += dofsize;
        }
        
    }

    // Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
    void CompMesh::Resequence(VecInt &DOFindices){
        
    }

    std::vector<double> &CompMesh::Solution(){
        return solution;
    }

    void CompMesh::LoadSolution(std::vector<double> &Sol){
        solution=Sol;
    }

    void CompMesh::Print (std::ostream & out) {
        
        out << "\n\t\t COMPUTABLE GRID INFORMATIONS:\n\n";
        
        out << "number of DOFs            = " << GetNumberDOF() << std::endl;
        out << "number of elements            = " << GetElementVec().size() << std::endl;
        out << "number of materials           = " << GetMathVec().size() << std::endl;
        out << "dimension of the mesh         = " <<  GetMathVec()[0]->Dimension() << std::endl;
        
        out << "\n\t Connect Information:\n\n";
        int64_t i, nelem = GetNumberDOF();
        for(i=0; i<nelem; i++) {
            out << " Index " << i << ' ';
            GetDOFVec()[i].Print(*this,out);
        }
        out << "\n\t Computable Element Information:\n\n";
        nelem = GetElementVec().size();
        for(i=0; i<nelem; i++) {
            if(!GetElement(i)) continue;
            CompElement *el = GetElement(i);
            out << "\n Index " << i << ' ';
            el->Print(out);
          }
        out << "\n\t Material Information:\n\n";
      //  std::map<int, MathStatement * >::const_iterator mit;
        nelem = GetMathVec().size();
        
        for(int mit=0; mit< nelem; mit++) {
            MathStatement *mat = GetMathVec()[mit];
            if (!mat) {
                DebugStop();
            }
            mat->Print(out);
        }
    }




