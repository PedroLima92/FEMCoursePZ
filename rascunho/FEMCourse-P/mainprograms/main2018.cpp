

#include <iostream>
#include <fstream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "VTKGeoMesh.h"
#include "ReadGmsh.h"
#include "tpanic.h"
#include "tpanic.h"
#include "GeomQuad.h"
#include "Geom1d.h"
#include "GeomTetrahedron.h"
#include "GeomTriangle.h"
#include "GeoElement.h"
#include "GeoElementTemplate.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "Assemble.h"
#include "Analysis.h"
#include "PostProcess.h"
#include "PostProcessTemplate.h"


using std::cout;
using std::endl;
using std::cin;

void TestIntegrate();
void TestGmsh();
void TestPoisson2DQuad(int order);
void TestPoisson2DTri(int order);
void TestPoisson3DTetra(int order);

void UXi(VecDouble &coord, VecDouble &uXi, VecDouble &gradu);
VecDouble X(VecDouble &coordXi);
void Jacobian(VecDouble &Coord, TMatrix &jacobian, double &detjac);
double InnerVec(VecDouble &S , VecDouble &T);

GeoMesh *CreateGMesh(int nx, int ny, double hx, double hy, ElementType type);
GeoMesh *CreateGMesh3D(int nx, int ny, int nz, double hx, double hy, double hz, ElementType type);
CompMesh *CMesh(GeoMesh *gmesh, int pOrder, int dim);

void F_source(const VecDouble &x, VecDouble &f);
void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol);
void ComputeRateOfConvergence(std::ostream &out, VecDouble &FirstError, VecDouble &erro, double nel);

int matId = 1;
int bc0 = -1; //define id para um material(cond contorno esq)
int bc1 = -2; //define id para um material(cond contorno dir)
int bc2 = -3; //define id para um material(cond contorno inf)
int bc3 = -4; //define id para um material(cond contorno sup)

int bc4 = -5; //define id para um material(cond contorno Frente)
int bc5 = -6; //define id para um material(cond contorno Atrás)

double hx = 1., hy=1.,hz=1.;
VecDouble FirstError(3,0.); // Vetor de erros
const double Pi=M_PI;

int main ()
{
//  1 - Verifica integração numérica
    
//  TestIntegrate();
    
//  2 - Geração de malha Gmsh, impressão VTK
    
//  TestGmsh()
 
//  3 - Poisson 2D, quadrilátero (linear e quadrático), obtenção de erros e taxas de convergência
    
//    TestPoisson2DQuad(1);
//    TestPoisson2DQuad(2);
    
//  4 - Poisson 2D, triângulo (linear e quadrático), obtenção de erros e taxas de convergência
    
//    TestPoisson2DTri(1);
//    TestPoisson2DTri(2);

//  5 - Poisson 3D, obtenção de erros e taxas de convergência
    
//    TestPoisson3DTetra(1);
    TestPoisson3DTetra(2);
    
    return 0;
}

void TestPoisson3DTetra(int pOrder){
    
    int ndiv = 5;
    
    for (int it=1; it<ndiv; it++) {
        
        double div = pow(2,it);
        // Malha geométrica :
        GeoMesh *geotest = CreateGMesh3D(div+1, div+1, div+1, hx, hy, hz, ETetraedro);
   //     geotest->Print(std::cout);
   //     VTKGeoMesh::PrintGMeshVTK(geotest, "MalhaTeste.vtk");
        
        // Malha computacional :
        CompMesh *cmesh = CMesh(geotest, pOrder,3);
        
        std::ofstream MalhaC("malahComp.txt");
  //      cmesh->Print(MalhaC);
        
        // Análise numérica :
        Analysis an(cmesh);
        an.RunSimulation();
        VecDouble Sol = cmesh->Solution();
        
        // Pós-processamento :
        PostProcess *solpos = new PostProcessTemplate<Poisson>(&an);
        solpos->SetExact(Sol_exact);
        solpos->AppendVariable("Sol");
        solpos->AppendVariable("Sol_exact");
        solpos->AppendVariable("Force");
        an.PostProcessSolution("SolutionPost.vtk", *solpos);
        
        // Calculo do erro :
        std::cout << "Comuting Error " << std::endl;
        VecDouble errovec(6,0.);
        std::ofstream ErroOut("Errors&Rates.txt", std::ofstream::app);
        errovec = an.PostProcessError(ErroOut, *solpos);
        
        // Taxa de convergência :
        
        ComputeRateOfConvergence(ErroOut, FirstError, errovec, div);
    }
    
}


void TestPoisson2DTri(int pOrder){
    
    int ndiv = 5;
    
    for (int it=1; it<ndiv; it++) {
        
        double div = pow(2,it);
        
        // Malha geométrica :
        GeoMesh *geotest = CreateGMesh(div+1, div+1, hx, hy, ETriangle);
        geotest->Print(std::cout);
        VTKGeoMesh::PrintGMeshVTK(geotest, "MalhaTeste.vtk");
        
        // Malha computacional :
        CompMesh *cmesh = CMesh(geotest, pOrder,2);
        
        // Análise numérica :
        Analysis an(cmesh);
        an.RunSimulation();
        VecDouble Sol = cmesh->Solution();
        
        // Pós-processamento :
        PostProcess *solpos = new PostProcessTemplate<Poisson>(&an);
        solpos->SetExact(Sol_exact);
        solpos->AppendVariable("Sol");
        solpos->AppendVariable("Sol_exact");
        solpos->AppendVariable("Force");
        an.PostProcessSolution("SolutionPost.vtk", *solpos);
        
        // Calculo do erro :
        std::cout << "Comuting Error " << std::endl;
        VecDouble errovec(6,0.);
        std::ofstream ErroOut("Errors&Rates.txt", std::ofstream::app);
        errovec = an.PostProcessError(ErroOut, *solpos);
        
        // Taxa de convergência :
        
        ComputeRateOfConvergence(ErroOut, FirstError, errovec, div);
    }
    
}


void TestPoisson2DQuad(int pOrder){
    
    int ndiv = 5;
    
    for (int it=1; it<ndiv; it++) {
        
        double div = pow(2,it);
        
        // Malha geométrica :
        GeoMesh *geotest = CreateGMesh(div+1, div+1, hx, hy, EQuadrilateral);
        geotest->Print(std::cout);
        VTKGeoMesh::PrintGMeshVTK(geotest, "MalhaTeste.vtk");
        
        // Malha computacional :
        CompMesh *cmesh = CMesh(geotest, pOrder,2);
        
        // Análise numérica :
        Analysis an(cmesh);
        an.RunSimulation();
        VecDouble Sol = cmesh->Solution();
        
        // Pós-processamento :
        PostProcess *solpos = new PostProcessTemplate<Poisson>(&an);
        solpos->SetExact(Sol_exact);
        solpos->AppendVariable("Sol");
        solpos->AppendVariable("Sol_exact");
        solpos->AppendVariable("Force");
        an.PostProcessSolution("SolutionPost.vtk", *solpos);
        
        // Calculo do erro :
        std::cout << "Comuting Error " << std::endl;
        VecDouble errovec(6,0.);
        std::ofstream ErroOut("Errors&Rates.txt", std::ofstream::app);
        errovec = an.PostProcessError(ErroOut, *solpos);
        
        // Taxa de convergência :
        
        ComputeRateOfConvergence(ErroOut, FirstError, errovec, div);
        
    }
    
}



void ComputeRateOfConvergence(std::ostream &out, VecDouble &FirstError, VecDouble &erro, double div){

    double h0 = 2.*hx/div;
    double h1 = hx/div;
    
    if (div==2) {
        for (int i = 0; i< erro.size(); i++) {
            FirstError[i] = erro[i];
        }
        return;
    }
    
    out << "----- Taxa de convergência : -----" << std::endl;
    out << "Norma L2 -> u = "<<(log(FirstError[0])-log(erro[0]))/(log(h0)-log(h1)) <<std::endl;
    out << "Norma L2 -> Grad u = "<<(log(FirstError[1])-log(erro[1]))/(log(h0)-log(h1)) <<std::endl;
    out << "Norma H1 -> u = "<<(log(FirstError[2])-log(erro[2]))/(log(h0)-log(h1)) <<std::endl;
    out << "---------------------------------- " << std::endl;
    
    for (int i = 0; i< erro.size(); i++) {
        FirstError[i] = erro[i];
    }
    
}

void F_source(const VecDouble &x, VecDouble &f){
    
    f.resize(2);
    
    double xv = x[0];
    double yv = x[1];
    // double zv = x[2];
    
    double f_x = + 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
    double f_y = - 8.0*Pi*Pi*cos(2.0*Pi*xv)*sin(2.0*Pi*yv);

    
    f[0] = f_x; // x direction
    f[1] = f_y; // y direction


}


void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol){
    
    dsol.Resize(3,3);
    sol.resize(3);
    dsol.Zero();
    
    double xv = x[0];
    double yv = x[1];
  //  double zv = x[2];
    
    double v_x =  cos(2*Pi*yv)*sin(2*Pi*xv);
    double v_y =  -(cos(2*Pi*xv)*sin(2*Pi*yv));
    
    sol[0]=v_x;
    sol[1]=v_y;
    
    // vx direction
    dsol(0,0)= 2*Pi*cos(2*Pi*xv)*cos(2*Pi*yv);
    dsol(0,1)= 2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);

    
    // vy direction
    dsol(1,0)= -2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
    dsol(1,1)= -2*Pi*cos(2*Pi*xv)*cos(2*Pi*yv);
    
}



CompMesh *CMesh(GeoMesh *gmesh, int pOrder, int dim){
 
    Matrix perm(dim,dim,0.), proj(dim,dim,0.);

    for (int i =0 ; i<dim; i++) {
        perm(i,i)=1.;
        proj(i,i)=1.;
    }

    CompMesh * cmesh = new CompMesh(gmesh);

    int nel= gmesh->NumElements();
    for (int iel = 0; iel<nel; iel++) {
        int geoMatID = gmesh->Element(iel)->Material();
        if(geoMatID==matId){
            // Materiais internos (Poisson)
            cmesh->SetNumberMath(iel+1);
            Poisson *material = new Poisson(geoMatID,perm);
            material->SetDimension(dim);
            material->SetForceFunction(F_source);
            material->SetExactSolution(Sol_exact);
            cmesh->SetMathStatement(iel, material);
            
        }else{
            // Condições de contorno (L2Projection)
            cmesh->SetNumberMath(iel+1);
            Matrix val1(2,1,0.), val2(2,1,0.);
            L2Projection *bcmat0 = new L2Projection(0,geoMatID,proj,val1,val2);
            bcmat0->SetExactSolution(Sol_exact);
            cmesh->SetMathStatement(iel, bcmat0);
        }
    
    }
    
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();
    
    return cmesh;
}


GeoMesh *CreateGMesh(int nx, int ny, double hx, double hy, ElementType eltype){
    
    GeoMesh *gmesh = new GeoMesh;

    
    int id, index;
    VecDouble coord(3,0.);
    int nnodes=nx*ny;
    gmesh->SetNumNodes(nnodes);
    
    for(int i = 0; i < ny; i++){
        for(int j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = (i)*hy/(ny - 1);
            gmesh->Node(id).SetCo(coord);
        }
    }
    //(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index)
    VecInt nodeind(nnodes);
    VecInt nodeindBC(2,0.);
    
    
    switch (eltype) {
        case EQuadrilateral:
        {
            for(int iq = 0; iq < (ny - 1); iq++){
                for(int jq = 0; jq < (nx - 1); jq++){
                    index = iq*(nx - 1)+ jq;
                    nodeind[0] = iq*ny + jq;
                    nodeind[1] = nodeind[0]+1;
                    nodeind[2] = nodeind[1]+(nx);
                    nodeind[3] = nodeind[0]+(nx);
                    GeoElement *gel = new GeoElementTemplate<GeomQuad>(nodeind,matId,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gel);
                }
            }
        }
            break;
        case ETriangle:
        {
            
            VecInt nodeindD(3,0);
            VecInt nodeindU(3,0);
  
            index =0;
            for(int iq = 0; iq < (ny - 1); iq++){
                for(int jq = 0; jq < (nx - 1); jq++){
                 //   index = (iq)*(nx - 1)+ (jq);
                    nodeindD[0] = (iq)*ny + (jq);
                    nodeindD[1] = nodeindD[0]+1;
                    nodeindD[2] = nodeindD[0]+nx;
                    GeoElement *gelD = new GeoElementTemplate<GeomTriangle>(nodeindD,matId,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelD);
                    
                    index++;
                    
                    nodeindU[0] = nodeindD[1];
                    nodeindU[1] = nodeindD[1]+nx;
                    nodeindU[2] = nodeindD[0]+nx;
                    GeoElement *gelU = new GeoElementTemplate<GeomTriangle>(nodeindU,matId,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelU);
                    
                    index++;
                }
            }
            index--;
        }
            break;
            
        default:
        {
            std::cout << "Element type is not valid" << std::endl;
            DebugStop();
        }
            break;
    }
    
    

    
    for(int iq = 0; iq < (ny - 1); iq++){
        for(int jq = 0; jq < (nx - 1); jq++){
            nodeind[0] = iq*ny + jq;
            nodeind[1] = nodeind[0]+1;
            nodeind[2] = nodeind[1]+(nx);
            nodeind[3] = nodeind[0]+(nx);
            //Condição esquerda:
            if(jq==0){
                index ++;
                nodeindBC[0]=nodeind[0];
                nodeindBC[1]=nodeind[3];
                GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(nodeindBC,bc0,gmesh,index);
                gmesh->SetNumElements(index+1);
                gmesh->SetElement(index, gelbc0);
            }
            //Condição direita:
            if(jq==nx-2){
                index ++;
                nodeindBC[0]=nodeind[1];
                nodeindBC[1]=nodeind[2];
                GeoElement *gelbc1 = new GeoElementTemplate<Geom1d>(nodeindBC,bc1,gmesh,index);
                gmesh->SetNumElements(index+1);
                gmesh->SetElement(index, gelbc1);
            }
            //Condição esquerda:
            if(iq==0){
                index ++;
                nodeindBC[0]=nodeind[0];
                nodeindBC[1]=nodeind[1];
                GeoElement *gelbc2 = new GeoElementTemplate<Geom1d>(nodeindBC,bc2,gmesh,index);
                gmesh->SetNumElements(index+1);
                gmesh->SetElement(index, gelbc2);
            }
            //Condição esquerda:
            if(iq==ny-2){
                index ++;
                nodeindBC[0]=nodeind[2];
                nodeindBC[1]=nodeind[3];
                GeoElement *gelbc3 = new GeoElementTemplate<Geom1d>(nodeindBC,bc3,gmesh,index);
                gmesh->SetNumElements(index+1);
                gmesh->SetElement(index, gelbc3);
            }
        }
    }
    
    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    return gmesh;
}

GeoMesh *CreateGMesh3D(int nx, int ny, int nz, double hx, double hy, double hz, ElementType eltype){
    
    GeoMesh *gmesh = new GeoMesh;
    
    int id, index;
    VecDouble coord(3,0.);
    int nnodes=nx*ny*nz;
    gmesh->SetNumNodes(nnodes);
    for (int k=0; k<nz; k++) {
        for(int i = 0; i < ny; i++){
            for(int j = 0; j < nx; j++){
                id = i*nx + j+ k*nx*ny;
                coord[0] = (j)*hx/(nx - 1);
                coord[1] = (i)*hy/(ny - 1);
                coord[2] = (k)*hz/(nz - 1);
                gmesh->Node(id).SetCo(coord);
            }
        }
    }

    
    VecInt nodeind(8);
    VecInt nodeindBC(3,0.);
    
    
    switch (eltype) {
        case ETetraedro:
        {
            
            VecInt nodeindD1(4,0);
            VecInt nodeindD2(4,0);
            VecInt nodeindU1(4,0);
            VecInt nodeindU2(4,0);
            VecInt nodeindL1(4,0);
            VecInt nodeindL2(4,0);
            
            index =0;
            for (int kq=0; kq<(nz-1); kq++) {
                for(int iq = 0; iq < (ny - 1); iq++){
                    for(int jq = 0; jq < (nx - 1); jq++){
                        
                        // Plano xy
                        
                        nodeindD1[0] = (iq)*ny + (jq) + kq*nx*ny;
                        nodeindD1[1] = nodeindD1[0]+1;
                        nodeindD1[2] = nodeindD1[0]+nx;
                        nodeindD1[3] = nodeindD1[1] + (1)*nx*ny;
                        GeoElement *gelD1 = new GeoElementTemplate<GeomTetrahedron>(nodeindD1,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelD1);
                    
                        index++;
                    
                        nodeindD2[0] = nodeindD1[1];
                        nodeindD2[1] = nodeindD1[2];
                        nodeindD2[2] = nodeindD1[1]+nx;
                        nodeindD2[3] = nodeindD1[1] + (1)*nx*ny;
                        GeoElement *gelD2 = new GeoElementTemplate<GeomTetrahedron>(nodeindD2,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelD2);

                        index++;

                        nodeindU1[0] = nodeindD1[0] + (1)*nx*ny;
                        nodeindU1[1] = nodeindU1[0]+1;
                        nodeindU1[2] = nodeindU1[0]+nx;
                        nodeindU1[3] = nodeindD1[2];
                        GeoElement *gelU1 = new GeoElementTemplate<GeomTetrahedron>(nodeindU1,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelU1);

                        index++;

                        nodeindU2[0] = nodeindU1[1];
                        nodeindU2[1] = nodeindU1[2];
                        nodeindU2[2] = nodeindU1[1]+nx;
                        nodeindU2[3] = nodeindD1[2];
                        GeoElement *gelU2 = new GeoElementTemplate<GeomTetrahedron>(nodeindU2,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelU2);

                        index++;

                        // Plano xz
                        
                        nodeindL1[0] = nodeindD1[0];
                        nodeindL1[1] = nodeindD1[2];
                        nodeindL1[2] = nodeindL1[0]+nx*ny;
                        nodeindL1[3] = nodeindU1[1];
                        GeoElement *gelL1 = new GeoElementTemplate<GeomTetrahedron>(nodeindL1,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelL1);
                        
                        index++;

                        
                        nodeindL2[0] = nodeindD2[2];
                        nodeindL2[1] = nodeindD1[1]+nx*ny;
                        nodeindL2[2] = nodeindU2[2];
                        nodeindL2[3] = nodeindU2[3];
                        GeoElement *gelL2 = new GeoElementTemplate<GeomTetrahedron>(nodeindL2,matId,gmesh,index);
                        gmesh->SetNumElements(index+1);
                        gmesh->SetElement(index, gelL2);

                        index++;
                        
                        
                    }
                }
            }
            index--;
        }
            break;
            
        default:
        {
            std::cout << "Element type is not valid" << std::endl;
            DebugStop();
        }
            break;
    }
    
    

    for (int kq = 0; kq < (nz-1); kq++) {
        for(int iq = 0; iq < (ny - 1); iq++){
            for(int jq = 0; jq < (nx - 1); jq++){
                nodeind[0] = (iq)*ny + (jq) + kq*nx*ny;
                nodeind[1] = nodeind[0]+1;
                nodeind[2] = nodeind[0]+(nx);
                nodeind[3] = nodeind[0]+(nx+1);

                nodeind[4] = nodeind[0]+nx*ny;
                nodeind[5] = nodeind[4]+1;
                nodeind[6] = nodeind[4]+(nx);
                nodeind[7] = nodeind[4]+(nx+1);
                

                //Condição esquerda:
                if(jq==0){
                    index ++;
                    nodeindBC[0]=nodeind[2];
                    nodeindBC[1]=nodeind[0];
                    nodeindBC[2]=nodeind[4];
                    GeoElement *gelbc00 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc0,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc00);
                    
                    index ++;
                    nodeindBC[0]=nodeind[2];
                    nodeindBC[1]=nodeind[6];
                    nodeindBC[2]=nodeind[4];
                    GeoElement *gelbc01 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc0,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc01);
                    
                }
                //Condição direita:
                if(jq==nx-2){
                    index ++;
                    nodeindBC[0]=nodeind[1];
                    nodeindBC[1]=nodeind[3];
                    nodeindBC[2]=nodeind[5];
                    GeoElement *gelbc10 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc1,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc10);
                    
                    index ++;
                    nodeindBC[0]=nodeind[3];
                    nodeindBC[1]=nodeind[7];
                    nodeindBC[2]=nodeind[5];
                    GeoElement *gelbc11 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc1,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc11);
                    
                }
                //Condição inf:
                if(kq==0){
                    index ++;
                    nodeindBC[0]=nodeind[2];
                    nodeindBC[1]=nodeind[3];
                    nodeindBC[2]=nodeind[1];
                    GeoElement *gelbc20 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc2,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc20);
                    
                    index ++;
                    nodeindBC[0]=nodeind[2];
                    nodeindBC[1]=nodeind[1];
                    nodeindBC[2]=nodeind[0];
                    GeoElement *gelbc21 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc2,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc21);
                    
                }
                //Condição sup:
                if(kq==ny-2){
                    index ++;
                    nodeindBC[0]=nodeind[4];
                    nodeindBC[1]=nodeind[5];
                    nodeindBC[2]=nodeind[6];
                    GeoElement *gelbc30 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc3,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc30);
                    
                    index ++;
                    nodeindBC[0]=nodeind[5];
                    nodeindBC[1]=nodeind[7];
                    nodeindBC[2]=nodeind[6];
                    GeoElement *gelbc31 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc3,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc31);
                }
                //Condição Frente:
                if(iq==0){
                    index ++;
                    nodeindBC[0]=nodeind[0];
                    nodeindBC[1]=nodeind[1];
                    nodeindBC[2]=nodeind[5];
                    GeoElement *gelbc40 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc4,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc40);
 
                    index ++;
                    nodeindBC[0]=nodeind[0];
                    nodeindBC[1]=nodeind[5];
                    nodeindBC[2]=nodeind[4];
                    GeoElement *gelbc41 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc4,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc41);
                    
                }
                //Condição Atrás:
                if(iq==ny-2){
                    index ++;
                    nodeindBC[0]=nodeind[3];
                    nodeindBC[1]=nodeind[2];
                    nodeindBC[2]=nodeind[7];
                    GeoElement *gelbc50 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc5,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc50);
                    
                    index ++;
                    nodeindBC[0]=nodeind[2];
                    nodeindBC[1]=nodeind[6];
                    nodeindBC[2]=nodeind[7];
                    GeoElement *gelbc51 = new GeoElementTemplate<GeomTriangle>(nodeindBC,bc5,gmesh,index);
                    gmesh->SetNumElements(index+1);
                    gmesh->SetElement(index, gelbc51);
                
                }
            }
        }
    }

    //Generate neighborhod information
    gmesh->BuildConnectivity();
    
    return gmesh;
}



void TestGmsh(){
    
    VecDouble vec1;
    ReadGmsh read;
    GeoMesh geomesh;
    read.Read(geomesh, "GeometryBench001.msh");
    VTKGeoMesh::PrintGMeshVTK(&geomesh, "MalhaTeste.vtk");
    geomesh.Print(std::cout);
    
}

void TestIntegrate()
{

    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"----------------- Regra de integração - Tabela de valores :  --------------------"<< endl;
    
    double weight;
    VecDouble CoordXi(2), uXi(1), gradu(2);
    
    //Teste 1 - Regra de integração - Tabela de valores;
    
    int order =1;
    double val1 = 0.;
    int NPoint;
    
    while(fabs(val1-12) >= 1.e-5) {
        IntRuleQuad TestInt1(order);
        NPoint = TestInt1.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt1.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val1 = val1 + weight*uXi[0];
        }
        order++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1  << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    //Teste 2 - Regra de integração - Tabela de valores;
    
    order =17;
    double val2 = 0.;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestInt2(order);
        NPoint = TestInt2.NPoints();
        for (int i=0; i<NPoint; i++) {
            TestInt2.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            val2 += weight*InnerVec(gradu, gradu);
        }
        val2=sqrt(val2);
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    GeoMesh mesh;
    CoordXi.clear();
    TMatrix jacobian, jacinv;
    double detjac;
    
    //    Teste 3 - Regra de integração - Tabela de valores;
    
    order =1;
    weight=0.;
    double valA=0;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    while(fabs(valA-80) >= 1.e-5) {
        IntRuleQuad TestInt3(order);
        
        VecDouble coord;
        
        NPoint = TestInt3.NPoints();
        VecInt nodes(NPoint);
        
        
        for (int i=0; i<NPoint; i++) {
            TestInt3.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            valA = valA + detjac*weight;
        }
    }
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    //Teste 4 - Regra de integração - Tabela de valores;
    
    order =26;
    weight=0.;
    double val3 = 0.;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    while(fabs(val3-239.49661609) >= 1.e-5) {
        IntRuleQuad TestInt4(order);
        NPoint = TestInt4.NPoints();
        VecInt nodes4(NPoint);
        for (int i=0; i<NPoint; i++) {
            TestInt4.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            val3 = val3 + detjac*weight*uXi[0];
        }
        order++;
    }
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< NPoint << " (Obs: para ordem > 19 o método NR Gauleg é chamado)"<< endl;
    
    //Teste 5 - Regra de integração - Tabela de valores;
    
    order =18;
    weight=0.;
    CoordXi[0]=CoordXi[1]=uXi[0]=gradu[0]=gradu[1]=0.;
    
    double val4 = 0.;
    
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestInt5(order);
        NPoint = TestInt5.NPoints();
        
        for (int i=0; i<NPoint; i++) {
            TestInt5.Point(i, CoordXi, weight);
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            val4 += detjac*weight*InnerVec(gradu, gradu);
        }
        
        val4=sqrt(val4);
    }
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< NPoint << endl;
    
    
    
    
    cout<<"---------------------------------------------------------------------------------"<< endl;
    cout<<"---------------- Função Gauleg retirada do Numerical Recipes : ------------------"<< endl;
    
    
    //Teste 1 - Função Gauleg retirada do Numerical Recipes;
    
    int order1=1;
    val1=0.;
    double nPoints = 1;
    
    while(fabs(val1-12.) >= 1.e-5) {
        IntRuleQuad TestIntGauss(order1);
        VecDouble x(nPoints);
        VecDouble wvec(nPoints);
        
        TestIntGauss.gaulegQuad(-1, 1, x, wvec);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x[i];
            CoordXi[1]=x[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            
            val1 += wvec[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(xi,eta):  "<< val1 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 2 - Função Gauleg retirada do Numerical Recipes;
    val2 = 0.;
    nPoints = 9;
    
    while(fabs(val2-5.03941698) >= 1.e-5) {
        IntRuleQuad TestIntGauss2(order1);
        VecDouble x2(nPoints);
        VecDouble wvec2(nPoints);
        
        TestIntGauss2.gaulegQuad(-1, 1, x2, wvec2);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x2[i];
            CoordXi[1]=x2[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            val2 += wvec2[i]*InnerVec(gradu, gradu);
        }
        
        val2=sqrt(val2);
        nPoints++;
    }
    cout<<"Resultado Integral ||Gradu(xi,eta)||:  "<< val2 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    
    //Teste 3 - Função Gauleg retirada do Numerical Recipes;
    
    valA=0.;
    nPoints=1;
    while(fabs(valA-80.) >= 1.e-5) {
        IntRuleQuad TestIntGauss3(order1);
        
        VecDouble x3(nPoints);
        VecDouble wvec3(nPoints);
        
        TestIntGauss3.gaulegQuad(-1, 1, x3, wvec3);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x3[i];
            CoordXi[1]=x3[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            valA += detjac*wvec3[i];
        }
        nPoints++;
    }
    
    cout<<"Resultado Área quad. mapeado:  "<< valA << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 4 - Função Gauleg retirada do Numerical Recipes;
    val3=0.;
    nPoints = 14;
    while(fabs(val3-239.49661609) >= 1.e-5) {
        IntRuleQuad TestIntGauss4(order1);
        VecDouble x4(nPoints);
        VecDouble wvec4(nPoints);
        
        TestIntGauss4.gaulegQuad(-1, 1, x4, wvec4);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x4[i];
            CoordXi[1]=x4[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val3 += detjac*wvec4[i]*uXi[0];
        }
        nPoints++;
    }
    
    cout<<"Resultado Integral u(x,y):  "<< val3 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    
    //Teste 5 - Função Gauleg retirada do Numerical Recipes;
    
    val4=0.;
    nPoints = 10;
    while(fabs(val4-22.53695786) >= 1.e-5) {
        IntRuleQuad TestIntGauss5(order1);
        VecDouble x5(nPoints);
        VecDouble wvec5(nPoints);
        
        TestIntGauss5.gaulegQuad(-1, 1, x5, wvec5);
        for (int i=0; i<nPoints*nPoints; i++) {
            CoordXi[0]=x5[i];
            CoordXi[1]=x5[i+nPoints*nPoints];
            UXi(CoordXi,uXi,gradu);
            //coord = X(CoordXi);
            Jacobian(CoordXi, jacobian, detjac);
            
            val4 += detjac*wvec5[i]*InnerVec(gradu, gradu);
        }
        val4=sqrt(val4);
        nPoints++;
    }
    
    cout<<"Resultado Integral ||Gradu(x,y)||:  "<< val4 << "  -> Número de pontos de integração = "<< nPoints*nPoints << endl;
    cout<<"-----------------------------------"<< endl;
    
}

void UXi(VecDouble &coord, VecDouble &uxi, VecDouble &gradu)
{
    uxi[0] = 3.+cos(3.*coord[1])*sin(4.*coord[0]);
    gradu[0] = 4.*cos(3.*coord[1])* cos(4.*coord[0]);
    gradu[1] = -3.*sin(3.*coord[1])* sin(4.*coord[0]);
}


double InnerVec(VecDouble &S , VecDouble &T){
    
    double Val = 0;
    if(S.size()!=T.size()){
        DebugStop();
    }
    for(int i = 0; i < S.size(); i++){
        Val += S[i]*T[i];
    }
    return Val;
    
}


VecDouble X(VecDouble &CoordXi)
{
    
    VecDouble coord(2);
    
    coord[0] = 5.*CoordXi[0]+0.5*sin(3.*CoordXi[1]);
    coord[1] = 4.*CoordXi[1]+0.3*cos(10.*CoordXi[0]);
    
    return coord;
    
}

void Jacobian(VecDouble &Coord, TMatrix &jacobian, double &detjac)
{
    
    jacobian.Resize(2, 2);
    jacobian.Zero();
    
    jacobian(0,0)=5.;
    jacobian(0,1)=1.5*cos(3.*Coord[1]);
    jacobian(1,0)=-3.*sin(10.*Coord[0]);
    jacobian(1,1)=4.;
    
    detjac = fabs(jacobian(0, 0)*jacobian(1, 1) - jacobian(1, 0)*jacobian(0, 1));
    
}
