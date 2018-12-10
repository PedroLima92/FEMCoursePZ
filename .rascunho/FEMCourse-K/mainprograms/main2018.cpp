#include <iostream>
#include "tpanic.h"

// Geometria e malha computacional
#include "GeoMesh.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "CompElementTemplate.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "Geom1d.h"

// Material
#include "Poisson.h"
#include "Elasticity2D.h"
#include <functional>

// Assemblagem
#include "Assemble.h"

// Condições de contorno
#include "L2Projection.h"

// Solução da malha em elementos finitos:
#include "Analysis.h"
#include "PostProcess.h"
#include "PostProcessTemplate.h"

// Pós-processamento
#include "VTKGeoMesh.h"

using std::cout;
using std::endl;
using std::cin;

int bcleft = -1;
int bcright = -2;
int bcbottom = -3;
int bctop = -4;
int bcfront = -5;
int bcback = -6;

void SolutionTest();
void solexact(const VecDouble &co, VecDouble &result, Matrix &deriv);
void solexactelast(const VecDouble &co, VecDouble &result, Matrix &deriv);
void forcefunction(const VecDouble &co, VecDouble &result);
void forcefunctionelast(const VecDouble &co, VecDouble &result);

GeoMesh * GenerateGeoMeshTrian(int64_t nelem_i);
GeoMesh * GenerateGeoMeshQuad(int64_t nelem_i);
GeoMesh * GenerateGeoMeshTetrahedron(int64_t nelem_i);
CompMesh * GenerateCompMesh(GeoMesh *gmesh, int pOrder);
GeoMesh *TetrahedronGeoMesh(int nnodes_x, int nnodes_y, int nnodes_z, double l);
CompMesh * GenerateCompMeshElasticity(GeoMesh *gmesh, int pOrder);
GeoMesh * GenerateGeoMeshQuadElasticity(int64_t nelem_i);

int main ()
{
    int pOrder = 2;
    VecDouble error0(3,0.);
    int nelem_0 = 0;
    
    bool poisson = true;
    
    if(poisson){
        for (int64_t iel = 0; iel<6; iel++){
            std::cout << iel << std::endl;
            int nelem_i = pow(2, iel);
//            GeoMesh * gmesh = GenerateGeoMeshTrian(nelem_i);
            GeoMesh * gmesh = GenerateGeoMeshQuad(nelem_i);
//            GeoMesh * gmesh = GenerateGeoMeshTetrahedron(nelem_i);
//            gmesh->Print(std::cout);
            
            VTKGeoMesh * testevtk = new VTKGeoMesh();
            testevtk->PrintGMeshVTK(gmesh, "teste.vtk");
            
            CompMesh * cmesh = GenerateCompMesh(gmesh, pOrder);

            Analysis * an = new Analysis(cmesh);
            an->RunSimulation();
            std::string filename("solution.vtk");
            testevtk->PrintCMeshVTK(cmesh, 2, filename);
            
            PostProcessTemplate<Poisson> defPostProc = PostProcessTemplate<Poisson>(an);
            defPostProc.SetExact(&solexact);
            defPostProc.AppendVariable(Poisson::ESol);
            defPostProc.AppendVariable(Poisson::EDSol);
            VecDouble error = an->PostProcessError(std::cout, defPostProc);
            
            double len_x = 1.0;
            std::cout << (log(error[0])-log(error0[0])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
            std::cout << (log(error[1])-log(error0[1])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
            std::cout << (log(error[2])-log(error0[2])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
            error0 = error;
            nelem_0 = nelem_i;
            
            delete testevtk;
            delete an;
            delete gmesh;
            delete cmesh;
        }
    }
    else{
//        for (int64_t iel = 1; iel<10; iel++){
//            std::cout << iel << std::endl;
////            int nelem_i = pow(2, iel);
//            int nelem_i = iel;
////            GeoMesh * gmesh = GenerateGeoMeshTrian(nelem_i);
//            GeoMesh * gmesh = GenerateGeoMeshQuadElasticity(nelem_i);
////            gmesh->Print(std::cout);
//            
//            VTKGeoMesh * testevtk = new VTKGeoMesh();
//            testevtk->PrintGMeshVTK(gmesh, "teste.vtk");
//            
//            CompMesh * cmesh = GenerateCompMeshElasticity(gmesh, pOrder);
//            
//            Analysis * an = new Analysis(cmesh);
//            an->RunSimulation();
//            
//            VecDouble sol = cmesh->Solution();
//
//            std::string filename("solution.vtk");
//            testevtk->PrintCMeshVTK(cmesh, 2, filename);
//            
//            PostProcessTemplate<Elasticity2D> defPostProc = PostProcessTemplate<Elasticity2D>(an);
//            defPostProc.SetExact(&solexactelast);
//            VecDouble error = an->PostProcessError(std::cout, defPostProc);
//            
//            double len_x = 1.0;
//            std::cout << (log(error[0])-log(error0[0])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
//            std::cout << (log(error[1])-log(error0[1])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
//            std::cout << (log(error[2])-log(error0[2])) / (log(len_x/nelem_i)-log(len_x/(nelem_0))) << std::endl;
//            error0 = error;
//            nelem_0 = nelem_i;
//            
//            delete testevtk;
//            delete an;
//            delete gmesh;
//            delete cmesh;
//        }
    }
    
    return EXIT_SUCCESS;
}


GeoMesh * GenerateGeoMeshTrian(int64_t nelem_i){
    GeoMesh *gmesh = new GeoMesh();
    gmesh->SetDimension(2);
    
    double len_x = 1.;
    double len_y = 1.;
    
    VecDouble error0(3,0.);
    VecDouble error(3,0.);
    
    double len_ix = len_x/(nelem_i);
    double len_iy = len_y/(nelem_i);
    int64_t nnodes = nelem_i+1;
    
    int64_t nelem = nelem_i*nelem_i;
    gmesh->SetNumElements(nelem*2+nelem_i*4);
    gmesh->SetNumNodes(nnodes*nnodes);
    
    // Enumeration: vertical order - from the below to the top, and from the left to the right
    VecDouble coord(3,0.);
    int64_t id;
    
    for(int i = 0; i < nnodes; i++){
        for(int j = 0; j < nnodes; j++){
            id = i*nnodes + j;
            coord[0] = (j)*len_ix;
            coord[1] = (i)*len_iy;
            gmesh->Node(id).SetCo(coord);
        }
    }
    
    id = 0;
    int materialid = 1;
    for(int i = 0; i < nelem_i; i++){
        for(int j = 0; j < nelem_i; j++){
            
            VecInt connect(3,0); // sentido anti-horario
            connect[0] = j*nnodes + i;
            connect[1] = connect[0] + 1;
            connect[2] = connect[0] + nnodes;
            
            GeoElement *gel0 = new GeoElementTemplate<GeomTriangle>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel0);
            id++;
            
            connect[0] = j*nnodes + i + 1;
            connect[1] = connect[0] + nnodes;
            connect[2] = j*nnodes + i + nnodes;
            GeoElement *gel1 = new GeoElementTemplate<GeomTriangle>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel1);
            id++;
        }
    }
    
    VecInt connect(2,0);
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcleft, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = nnodes-1;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcright, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = 0;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcbottom, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = nnodes*nnodes-nnodes;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bctop, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

GeoMesh * GenerateGeoMeshQuad(int64_t nelem_i){
    GeoMesh *gmesh = new GeoMesh();
    gmesh->SetDimension(2);
    
    double len_x = 1.;
    double len_y = 1.;
    
    VecDouble error0(3,0.);
    VecDouble error(3,0.);
    
    double len_ix = len_x/(nelem_i);
    double len_iy = len_y/(nelem_i);
    int64_t nnodes = nelem_i+1;
    
    int64_t nelem = nelem_i*nelem_i;
//    gmesh->SetNumElements(nelem+nelem_i*4);
    gmesh->SetNumNodes(nnodes*nnodes);
    
    // Enumeration: vertical order - from the below to the top, and from the left to the right
    VecDouble coord(3,0.);
    int64_t id;
    
    for(int i = 0; i < nnodes; i++){
        for(int j = 0; j < nnodes; j++){
            id = i*nnodes + j;
            coord[0] = (j)*len_ix;
            coord[1] = (i)*len_iy;
            gmesh->Node(id).SetCo(coord);
        }
    }
    
    int materialid = 1;
    for(int i = 0; i < nelem_i; i++){
        for(int j = 0; j < nelem_i; j++){
            int64_t id = i*nelem_i+j;
            
            VecInt connect(4,0); // sentido anti-horario
            connect[0] = j*nnodes + i;
            connect[1] = connect[0]+1;
            connect[3] = connect[0]+nnodes;
            connect[2] = connect[3]+1;
            
            GeoElement *gel = new GeoElementTemplate<GeomQuad>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel);
        }
    }
    
    int64_t elempos = nelem;
    VecInt connect(2,0);
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcleft, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    connect[1] = nnodes-1;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcright, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;    }
    
    connect[1] = 0;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcbottom, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    connect[1] = nnodes*nnodes-nnodes;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bctop, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

GeoMesh * GenerateGeoMeshTetrahedron(int64_t nelem_i){
    GeoMesh * gmesh = new GeoMesh();
    gmesh->SetDimension(3);
    
    double len_x = 1.;
    double len_y = 1.;
    double len_z = 1.;
    
    VecDouble error0(3,0.);
    VecDouble error(3,0.);
    
    double len_ix = len_x/(nelem_i);
    double len_iy = len_y/(nelem_i);
    double len_iz = len_z/(nelem_i);
    int64_t nnodes = nelem_i+1;
    
    int64_t nelem = nelem_i*nelem_i*nelem_i;
    gmesh->SetNumElements(nelem*6 + 6*nelem_i*nelem_i*2);
    gmesh->SetNumNodes(pow(nnodes, 3));
    
    // Enumeration: vertical order - from the below to the top, and from the left to the right
    VecDouble coord(3,0.);
    int64_t id = 0;
    
    for(int i = 0; i < nnodes; i++){
        for(int j = 0; j < nnodes; j++){
            for (int k = 0; k < nnodes; k++) {
                coord[0] = (j)*len_ix;
                coord[1] = (i)*len_iy;
                coord[2] = (k)*len_iz;
                gmesh->Node(id).SetCo(coord);
                id++;
            }
        }
    }
    
    int materialid = 1;
    id = 0;
    for(int i = 0; i < nelem_i; i++){
        for (int j = 0; j < nelem_i; j++) {
            for (int k = 0; k < nelem_i; k++) {
                int64_t index = i*nnodes*nnodes+j*nnodes+k;
                
                VecInt connect(4,0); // sentido anti-horario
                connect[0] = index;
                connect[1] = connect[0]+nnodes;
                connect[2] = connect[0]+1;
                connect[3] = connect[0]+nnodes*nnodes+1;
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
                
                connect[0] = connect[0];
                connect[1] = connect[0]+nnodes;
                connect[2] = connect[0]+nnodes*nnodes+1;
                connect[3] = connect[0]+nnodes*nnodes;
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
                
                connect[0] = nnodes+index+1;
                connect[1] = connect[0]-nnodes;
                connect[2] = connect[0]-1;
                connect[3] = connect[0]+nnodes*nnodes-nnodes;
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
                
                connect[0] = nnodes+index+1;
                connect[1] = connect[0]+nnodes*nnodes-nnodes;
                connect[2] = connect[0]-1;
                connect[3] = connect[0]+nnodes*nnodes;
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
                
                connect[0] = nnodes*nnodes+nnodes+index;
                connect[1] = connect[0]-nnodes;
                connect[2] = connect[1]+1;
                connect[3] = connect[0]-nnodes*nnodes;
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
                
                connect[0] = nnodes*nnodes+nnodes+index;
                connect[1] = connect[0]-nnodes+1;
                connect[2] = connect[0]+1;
                connect[3] = connect[0]-nnodes*nnodes;
                
                new GeoElementTemplate<GeomTetrahedron>(connect, materialid, gmesh, id);
                id++;
            }
            
        }
        
    }
    
    int64_t elempos = nelem*6;
    VecInt connect(3,0);
    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = iel*nnodes+jel;
            connect[1] = connect[0]+1;
            connect[2] = connect[0]+nnodes;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bcleft, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[1];
            connect[1] = connect[0]+nnodes;
            connect[2] = connect[2];
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bcleft, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }
    
    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = nnodes*nnodes*(nnodes-1)+iel*nnodes+jel;
            connect[1] = connect[0]+1;
            connect[2] = connect[0]+nnodes;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bcright, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[1];
            connect[1] = connect[0]+nnodes;
            connect[2] = connect[2];
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bcright, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }

    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = iel*nnodes*nnodes+jel;
            connect[2] = connect[0]+1;
            connect[1] = connect[0]+nnodes*nnodes+1;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bcbottom, gmesh, elempos); // na vdd é bcback
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[0];
            connect[2] = connect[1];
            connect[1] = connect[2]-1;
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bcbottom, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }
    
    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = iel*nnodes*nnodes+jel*nnodes;
            connect[1] = connect[0]+nnodes;
            connect[2] = connect[0]+nnodes*nnodes;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bctop, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[1];
            connect[2] = connect[2];
            connect[1] = connect[2]+nnodes;
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bctop, gmesh, elempos); // na vdd é bcbottom
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }
    
    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = nnodes*(nnodes-1) + nnodes*nnodes*iel+jel;
            connect[1] = connect[0]+nnodes*nnodes+1;
            connect[2] = connect[0]+1;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bcfront, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[0];
            connect[2] = connect[1];
            connect[1] = connect[2]-1;
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bcfront, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }
    
    for (int64_t iel=0; iel<nelem_i; iel++) {
        for (int64_t jel=0; jel<nelem_i; jel++) {
            connect[0] = nnodes + nnodes*nnodes*iel + jel*nnodes-1;
            connect[1] = connect[0]+nnodes;
            connect[2] = connect[0]+nnodes*nnodes;
            GeoElement *gelbc0 = new GeoElementTemplate<GeomTriangle>(connect, bcback, gmesh, elempos);
            gmesh->SetElement(elempos, gelbc0);
            elempos++;
            
            connect[0] = connect[1];
            connect[2] = connect[2];
            connect[1] = connect[2]+nnodes;
            GeoElement *gelbc1 = new GeoElementTemplate<GeomTriangle>(connect, bcback, gmesh, elempos); // na vdd é bctop
            gmesh->SetElement(elempos, gelbc1);
            elempos++;
        }
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

CompMesh * GenerateCompMesh(GeoMesh *gmesh, int pOrder){
    // Material:
    int dim = gmesh->Dimension();
    Matrix perm(dim, dim, 0.);
    for (int64_t idim=0; idim<dim; idim++) {
        perm(idim,idim) = 1.;
    }
    int materialid = 1;
    Poisson * pos = new Poisson(materialid, perm);
    pos->SetDimension(gmesh->Dimension());
    
    pos->SetForceFunction(&forcefunction);
    
    // Condição de contorno:
    Matrix projection(2,2,0.);
    int bctype = 0;
    Matrix Val1(2,1,0.);
    Matrix Val2(2,1,0.);
    
    L2Projection * bc0 = new L2Projection(bctype, bcleft, projection, Val1, Val2);
    bc0->SetExactSolution(&solexact);
    bc0->SetDimension(gmesh->Dimension()-1);
    
    L2Projection * bc1 = new L2Projection(bctype, bcright, projection, Val1, Val2);
    bc1->SetExactSolution(&solexact);
    bc1->SetDimension(gmesh->Dimension()-1);
    
    L2Projection * bc2 = new L2Projection(bctype, bcbottom, projection, Val1, Val2);
    bc2->SetExactSolution(&solexact);
    bc2->SetDimension(gmesh->Dimension()-1);
    
    L2Projection * bc3 = new L2Projection(bctype, bctop, projection, Val1, Val2);
    bc3->SetExactSolution(&solexact);
    bc3->SetDimension(gmesh->Dimension()-1);
    
    L2Projection * bc4 = new L2Projection(bctype, bcfront, projection, Val1, Val2);
    bc4->SetExactSolution(&solexact);
    bc4->SetDimension(gmesh->Dimension()-1);
    
    L2Projection * bc5 = new L2Projection(bctype, bcback, projection, Val1, Val2);
    bc5->SetExactSolution(&solexact);
    bc5->SetDimension(gmesh->Dimension()-1);
    
    // Criando malha computacional:
    CompMesh * cmesh = new CompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetNumberMath(7);
    //    cmesh->SetNumberMath(5);
    cmesh->SetMathStatement(0, pos);
    cmesh->SetMathStatement(1, bc0);
    cmesh->SetMathStatement(2, bc1);
    cmesh->SetMathStatement(3, bc2);
    cmesh->SetMathStatement(4, bc3);
    cmesh->SetMathStatement(5, bc4);
    cmesh->SetMathStatement(6, bc5);
    cmesh->AutoBuild();
    
    return cmesh;
}

GeoMesh * GenerateGeoMeshTrianElasticity(int64_t nelem_i){
    GeoMesh *gmesh = new GeoMesh();
    gmesh->SetDimension(2);
    
    double len_x = 1.;
    double len_y = 1.;
    
    VecDouble error0(3,0.);
    VecDouble error(3,0.);
    
    double len_ix = len_x/(nelem_i);
    double len_iy = len_y/(nelem_i);
    int64_t nnodes = nelem_i+1;
    
    int64_t nelem = nelem_i*nelem_i;
    gmesh->SetNumElements(nelem*2+nelem_i*4);
    gmesh->SetNumNodes(nnodes*nnodes);
    
    // Enumeration: vertical order - from the below to the top, and from the left to the right
    VecDouble coord(3,0.);
    int64_t id;
    
    for(int i = 0; i < nnodes; i++){
        for(int j = 0; j < nnodes; j++){
            id = i*nnodes + j;
            coord[0] = -0.5 + (j)*len_ix;
            coord[1] = -0.5 + (i)*len_iy;
            gmesh->Node(id).SetCo(coord);
        }
    }
    
    id = 0;
    int materialid = 1;
    for(int i = 0; i < nelem_i; i++){
        for(int j = 0; j < nelem_i; j++){
            
            VecInt connect(3,0); // sentido anti-horario
            connect[0] = j*nnodes + i;
            connect[1] = connect[0] + 1;
            connect[2] = connect[0] + nnodes;
            
            GeoElement *gel0 = new GeoElementTemplate<GeomTriangle>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel0);
            id++;
            
            connect[0] = j*nnodes + i + 1;
            connect[1] = connect[0] + nnodes;
            connect[2] = j*nnodes + i + nnodes;
            GeoElement *gel1 = new GeoElementTemplate<GeomTriangle>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel1);
            id++;
        }
    }
    
    VecInt connect(2,0);
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcleft, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = nnodes-1;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcright, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = 0;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcbottom, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    connect[1] = nnodes*nnodes-nnodes;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bctop, gmesh, id);
        gmesh->SetElement(id, gelbc0);
        id++;
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

GeoMesh * GenerateGeoMeshQuadElasticity(int64_t nelem_i){
    GeoMesh *gmesh = new GeoMesh();
    gmesh->SetDimension(2);
    
    double len_x = 1.;
    double len_y = 1.;
    
    VecDouble error0(3,0.);
    VecDouble error(3,0.);
    
    double len_ix = len_x/(nelem_i);
    double len_iy = len_y/(nelem_i);
    int64_t nnodes = nelem_i+1;
    
    int64_t nelem = nelem_i*nelem_i;
    gmesh->SetNumElements(nelem+nelem_i*4);
    gmesh->SetNumNodes(nnodes*nnodes);
    
    // Enumeration: vertical order - from the below to the top, and from the left to the right
    VecDouble coord(3,0.);
    int64_t id;
    
    for(int i = 0; i < nnodes; i++){
        for(int j = 0; j < nnodes; j++){
            id = i*nnodes + j;
            coord[0] = (j)*len_ix;
            coord[1] = (i)*len_iy;
            gmesh->Node(id).SetCo(coord);
        }
    }
    
    int materialid = 1;
    for(int i = 0; i < nelem_i; i++){
        for(int j = 0; j < nelem_i; j++){
            int64_t id = i*nelem_i+j;
            
            VecInt connect(4,0); // sentido anti-horario
            connect[0] = j*nnodes + i;
            connect[1] = connect[0]+1;
            connect[3] = connect[0]+nnodes;
            connect[2] = connect[3]+1;
            
            GeoElement *gel = new GeoElementTemplate<GeomQuad>(connect, materialid, gmesh, id);
            gmesh->SetElement(id, gel);
        }
    }
    
    int64_t elempos = nelem;
    VecInt connect(2,0);
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcleft, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    connect[1] = nnodes-1;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+nnodes;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcright, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    connect[1] = 0;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcbottom, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    connect[1] = nnodes*nnodes-nnodes;
    for (int64_t iel=0; iel<nelem_i; iel++) {
        connect[0] = connect[1];
        connect[1] = connect[0]+1;
        GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bctop, gmesh, elempos);
        gmesh->SetElement(elempos, gelbc0);
        elempos++;
    }
    
    gmesh->BuildConnectivity();
    
    return gmesh;
    
}

CompMesh * GenerateCompMeshElasticity(GeoMesh *gmesh, int pOrder){
    // Material:
    Matrix perm(2, 2, 1.);
    perm(1,0) = 0.;
    perm(0,1) = 0.;
    int materialid = 1;
    
    Elasticity2D * elast = new Elasticity2D(materialid,1.,0.3);
    elast->SetForceFunction(&forcefunctionelast);
    elast->SetPlaneStressState(1);
    elast->SetDimension(2);
    
    // Condição de contorno:
    Matrix projection(2,2,0.);
    int bctype = 0;
    Matrix Val1(2,1,0.);
    Matrix Val2(2,1,0.);
    
    L2Projection * bc0 = new L2Projection(bctype, bcleft, projection, Val1, Val2);
    bc0->SetExactSolution(&solexactelast);
    bc0->SetDimension(1);
    
    L2Projection * bc1 = new L2Projection(bctype, bcright, projection, Val1, Val2);
    bc1->SetExactSolution(&solexactelast);
//    bc1->SetForceFunction(&forcefunctionelast);
    bc1->SetDimension(1);
    
    L2Projection * bc2 = new L2Projection(bctype, bcbottom, projection, Val1, Val2);
    bc1->SetExactSolution(&solexactelast);
//    bc2->SetForceFunction(&forcefunctionelast);
    bc2->SetDimension(1);
    
    L2Projection * bc3 = new L2Projection(bctype, bctop, projection, Val1, Val2);
    bc1->SetExactSolution(&solexactelast);
//    bc3->SetForceFunction(&forcefunctionelast);
    bc3->SetDimension(1);
    
    // Criando malha computacional:
    CompMesh * cmesh = new CompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetNumberMath(5);
    cmesh->SetMathStatement(0, elast);
    cmesh->SetMathStatement(1, bc0);
    cmesh->SetMathStatement(2, bc1);
    cmesh->SetMathStatement(3, bc2);
    cmesh->SetMathStatement(4, bc3);
    cmesh->AutoBuild();
    
    return cmesh;
}

void forcefunction(const VecDouble &co, VecDouble &result)
{
    result.resize(2);
    
    double xv = co[0];
    double yv = co[1];
    
    double Pi = M_PI;
    result[0] = 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
    result[1] = -8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
}

void solexact(const VecDouble &co, VecDouble &result, Matrix &deriv)
{
    result.resize(2);
    deriv.Resize(2,2);
    
    double xv = co[0];
    double yv = co[1];
    
    double Pi = M_PI;
    
    double v_x = cos(2*Pi*yv)*sin(2*Pi*xv);
    
    result[0] = v_x;
    result[1] = -v_x;
    
    deriv(0,0)=  2*Pi*cos(2*Pi*yv)*cos(2*Pi*xv);
    deriv(1,0)= -2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
    deriv(0,1)= -2*Pi*cos(2*Pi*yv)*cos(2*Pi*xv);
    deriv(1,1)=  2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
}

void solexactelast(const VecDouble &co, VecDouble &result, Matrix &deriv)
{
    double xv = co[0];
    double yv = co[1];
    
//    result[0] = 3*xv;
//    result[1] = -10*yv;
//
//    deriv.Zero();
//    deriv(0,0) = 3.;
//    deriv(0,1) = 0;
//    deriv(1,0) = 0;
//    deriv(1,1) = -10.;
    
    result[0] = 1.*xv*cos((M_PI*yv)/4);
    result[1] = -0.381972*sin((M_PI*yv)/4);
    
    deriv.Resize(2,3);
    deriv.Zero();
    deriv(0,0) = cos(M_PI*yv/4);
    deriv(1,0) = -0.785398*xv*sin((M_PI*yv)/4);
    deriv(0,1) = 0.;
    deriv(1,1) = -0.3*cos((M_PI*yv)/4);
}

void forcefunctionelast(const VecDouble &co, VecDouble &result)
{
    double yv = co[1];
    
    result[0] = cos(M_PI*yv/4);
    result[1] = 0.;
    
}
