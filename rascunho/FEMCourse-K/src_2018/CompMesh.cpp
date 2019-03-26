//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#include "CompMesh.h"
#include "CompElement.h"
#include "GeoMesh.h"
#include "MathStatement.h"

CompMesh::CompMesh() : geomesh(0), compelements(0), dofs(0), mathstatements(0), solution(0) {
    DefaultOrder = 1;
}

CompMesh::CompMesh(GeoMesh *gmesh) : geomesh(gmesh), compelements(0), dofs(0), mathstatements(0), solution(0)
{
    geomesh = gmesh;
    SetNumberElement(gmesh->NumElements());
    gmesh->SetReference(this);
    DefaultOrder = -1;
    
};

// Copy constructor of CompMesh
CompMesh::CompMesh(const CompMesh &copy){
    geomesh = copy.geomesh;
    compelements = copy.compelements;
    dofs = copy.dofs;
    mathstatements = copy.mathstatements;
    solution = copy.solution;
    DefaultOrder = copy.DefaultOrder;
}

// Destructor of CompMesh
CompMesh::~CompMesh(){
    
}

GeoMesh *CompMesh::GetGeoMesh() const{
    return geomesh;
}

void CompMesh::SetGeoMesh(GeoMesh *gmesh){
    geomesh = gmesh;
}
    
// Set the number of computational elements on the grid
void CompMesh::SetNumberElement(int64_t nelem){
    compelements.resize(nelem);
}
    
// Set the number of degrees of freedom
void CompMesh::SetNumberDOF(int64_t ndof){
    dofs.resize(ndof);
}
    
// Set the number of math statements
void CompMesh::SetNumberMath(int nmath){
    mathstatements.resize(nmath);
}

// Set the computational element associated to an index
void CompMesh::SetElement(int64_t elindex, CompElement *cel){
    compelements[elindex] = cel;
}
    
// Set the degree of freedom associated to an index
void CompMesh::SetDOF(int64_t index, const DOF &dof){
    dofs[index] = dof;
}
    
// Set the math statement object associated to an index
void CompMesh::SetMathStatement(int index, MathStatement *math){
    mathstatements[index] = math;
}
    
// Return the degree of freedom index
DOF &CompMesh::GetDOF(int64_t dofindex){
    return dofs[dofindex];
}
    
// Return the computational element associated to an index
CompElement *CompMesh::GetElement(int64_t elindex) const{
    return compelements[elindex];
}
    
// Return the math statement object associated to an index
MathStatement *CompMesh::GetMath(int matindex) const{
    return mathstatements[matindex];
}

// Return the vector with computational elements
std::vector<CompElement *> CompMesh::GetElementVec() const{
    return compelements;
}
    
// Return the vector with degrees of freedom
std::vector<DOF> CompMesh::GetDOFVec() const{
    return dofs;
}
    
// Return the vector with math statement objects
std::vector<MathStatement *> CompMesh::GetMathVec() const{
    return mathstatements;
}

// Set the vector with computational elements
void CompMesh::SetElementVec(const std::vector<CompElement *> &vec){
    compelements = vec;
}
    
// Set the vector with degrees of freedom
void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec){
    dofs = dofvec;
}
    
// Set the vector with math statement objects
void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec){
    mathstatements = mathvec;
}

// will create the computational elements
void CompMesh::AutoBuild(){
    int nelem = geomesh->NumElements();
    
    int id = 0;
    for (int64_t icompel=0; icompel<nelem; icompel++) {
        GeoElement * gel = geomesh->Element(icompel);
        if (!gel) {
            continue;
        }
        gel->CreateCompEl(this, id);
        id++;
        
    }
    this->Resequence();
    
    std::cout << "Computational Mesh: AutoBuild done!" << std::endl;
}

// Initialize the datastructure FirstEquation of the DOF objects
void CompMesh::Resequence(){
    
    int64_t nelem = compelements.size();
    int64_t first = 0;
    int64_t pos = 0;
    
    for (int64_t iel=0; iel<nelem; iel++) {
        CompElement * cel = this->GetElement(iel);
        if (!cel || !(cel->GetGeoElement())) {
            continue;
        }
        
        int nsides = cel->GetGeoElement()->NSides();

        for (int64_t isides=0; isides<nsides; isides++) {
            int64_t idof = cel->GetDOFIndex(isides);
            
            if (idof >= pos) {
                if (dofs[pos].GetNShape() > 0) {
                    dofs[pos].SetFirstEquation(first);
                    first += dofs[idof].GetNShape()*dofs[idof].GetNState();
                }
                pos ++;
            }
        }
    }
}

// Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
void CompMesh::Resequence(VecInt &DOFindices){
    std::cout << "CompMesh::Resequence :  Not implemented yet" << std::endl;
    
    int nelem = compelements.size();
    int64_t ndofs = DOFindices.size();
    int64_t first = 0;
    
    for (int i = 0; i < nelem; i++) {
        for (int j = 0; j < ndofs; j++) {
            int idof = DOFindices[j];
            dofs[idof].SetFirstEquation(first);
            int result = GetDOF(idof).GetNShape() * GetDOF(idof).GetNState();
            first += result;
        }
    }
}

std::vector<double> &CompMesh::Solution(){
    return solution;
}

void CompMesh::LoadSolution(std::vector<double> &Sol){
    solution = Sol;
}

