//#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include <stdio.h>
//#include <functional>
//#include <math.h>
//#include "tpanic.h"
//
//// Integration rules
//#include "IntRule.h"
//#include "IntRule1d.h"
//#include "IntRuleQuad.h"
//#include "IntRuleTriangle.h"
//#include "IntRuleTetrahedron.h"
//
//// Leitura da malha e pós processamento com VTK
//#include "ReadGmsh.h"
//#include "VTKGeoMesh.h"
//
//// Geometria e malha computacional
//#include "GeoMesh.h"
//#include "GeoElement.h"
//#include "GeoElementTemplate.h"
//#include "CompElementTemplate.h"
//#include "GeomQuad.h"
//#include "Geom1d.h"
//
//// Material
//#include "Poisson.h"
//#include <functional>
//
//// Assemblagem
//#include "Assemble.h"
//
//// Condições de contorno
//#include "L2Projection.h"
//
//// Solução da malha em elementos finitos:
//#include "Analysis.h"
//#include "PostProcess.h"
//#include "PostProcessTemplate.h"
//
//#include <pzlog.h>
//
//using std::cout;
//using std::endl;
//using std::cin;
//
//// Exercício 3:
//void Exercicio3 (int maxpOrder, double tol);
//double Funcao2D (VecDouble &pt);
//double Funcao2D_NormadaDerivada (VecDouble &pt);
//double DetJacobian (VecDouble &pt);
//
//// Testes de integração:
//void IntegrationTests();
//void test1D();
//void testQuad();
//void testTriangle();
//
//// Testes gmesh e cmesh:
//void GmeshCmeshTests();
//
//// Teste matriz de rigidez para 1 elemento:
//void StiffnessMatrix1Test();
//void forcefunction(const VecDouble &co, VecDouble &result);
//
//// Teste da assemblagem:
//void AssembleTest();
//
//// Teste de aplicação das condições de contorno, solução da malha computacional e cálculo do erro
//void SolutionTest();
//void solexact(const VecDouble &co, VecDouble &result, Matrix &deriv);
//
//int main ()
//{
//    //    // Resolução do Exercício 3
//    //    int maxpOrder = 40;
//    //    double tol = 1e-6;
//    //    Exercicio3(maxpOrder, tol);
//    //
//    //    // Testes das regras de integração:
//    //    IntegrationTests();
//    //
//    //    // Teste de leitura da malha, build connectivity e Paraview
//    //    GmeshCmeshTests();
//    //
//    //    // Teste da matriz de rigidez
//    //    StiffnessMatrix1Test();
//    //
//    //    // Teste do Assemble:
//    //    AssembleTest();
//    //
//    // Aplicação das condições de contorno:
//    SolutionTest();
//    
//    return EXIT_SUCCESS;
//}
//
//// ----------------------- EXERCÍCIO 3 -----------------------
//
//void Exercicio3(int maxpOrder, double tol){
//    
//    std::string namefile = "Results";
//    std::ofstream outfile(namefile+".txt");
//    
//    double integralxieta = 0.;
//    double integralxieta_norm = 0.;
//    double integralxy = 0.;
//    double integralxy_norm = 0.;
//    double integralarea = 0.;
//    
//    int pOrder = 1;
//    
//    while ( fabs(integralxieta - 12.0)>tol || fabs(integralxieta_norm - 5.0394169797223975)>tol || fabs(integralxy - 239.49661608858898)>tol || fabs(integralxy_norm - 22.53695786372607)>tol || fabs(integralarea - 80.0)>tol)
//    {
//        integralxieta = 0.;
//        integralxieta_norm = 0.;
//        integralxy = 0.;
//        integralxy_norm = 0.;
//        integralarea = 0.;
//        IntRuleQuad * integ = new IntRuleQuad(pOrder);
//        
//        int64_t npts = integ->NPoints();
//        VecDouble point(3,0.);
//        
//        for (int64_t ipts=0; ipts<npts; ipts++) {
//            
//            double weight;
//            VecDouble weights;
//            
//            integ->Point(ipts, point, weight);
//            
//            integralxieta += weight * Funcao2D(point);
//            integralxieta_norm += weight * Funcao2D_NormadaDerivada(point);
//            integralxy += weight * Funcao2D(point)*DetJacobian(point);
//            integralxy_norm += weight * Funcao2D_NormadaDerivada(point)*DetJacobian(point);
//            integralarea += weight * DetJacobian(point);
//        }
//        std::string namefile_int = "Integr_Order"+std::to_string(pOrder);
//        std::ofstream outfile_int(namefile_int+".txt");
//        integ->Print(outfile_int);
//        
//        outfile << "\nOrdem de integracao: " << pOrder << "; Npts: " << npts << std::endl;
//        outfile << "Integral xi, eta:" << integralxieta << "; Norma xi,etx:"<< sqrt(integralxieta_norm) << std::endl;
//        outfile << "Integral x,y:" << integralxy << "; Norma x,y:"<< sqrt(integralxy_norm) << std::endl;
//        outfile << "Area:" << integralarea << std::endl;
//        
//        std::cout << "\nOrdem de integracao: " << pOrder << "; Npts: " << npts << std::endl;
//        std::cout << "Integral xi, eta:" << integralxieta << "; Norma xi,etx:"<< sqrt(integralxieta_norm) << std::endl;
//        std::cout << "Integral x,y:" << integralxy << "; Norma x,y:"<< sqrt(integralxy_norm) << std::endl;
//        std::cout << "Area:" << integralarea << std::endl;
//        
//        pOrder++;
//        
//        if (pOrder > maxpOrder){
//            break;
//        }
//        
//    }
//    
//    if (pOrder >= maxpOrder) {
//        std::cout << "Not converged for tolerance = " << tol << std::endl;
//        DebugStop();
//    }
//    else{
//        std::cout << "Finished successful!" << std::endl;
//    }
//}
//
//double Funcao2D (VecDouble &pt) {
//    return (3 + sin(4*pt[0])*cos(3*pt[1]));
//}
//
//double Funcao2D_NormadaDerivada (VecDouble &pt) {
//    return ( 4*cos(3*pt[1])*cos(4*pt[0])*4*cos(3*pt[1])*cos(4*pt[0]) + (-3*sin(3*pt[1])*sin(4*pt[0]))*(-3*sin(3*pt[1])*sin(4*pt[0])) );
//}
//
//double DetJacobian (VecDouble &pt) {
//    VecDouble x(2,0.);
//    
//    x[0] = 5*pt[0] + 0.5*sin(3*pt[1]);
//    x[1] = 4*pt[1] + 0.3*cos(3*pt[0]);
//    
//    Matrix jac(2, 2, 0.);
//    jac(0,0) = 5;
//    jac(0,1) = 1.5*cos(3*pt[1]);
//    jac(1,0) = -3.*sin(10*pt[0]);
//    jac(1,1) = 4;
//    
//    return jac(0,0)*jac(1,1)-jac(0,1)*jac(1,0);
//}
//
//// ----------------- TESTES DE INTEGRAÇÃO --------------------
//
//void IntegrationTests(){
//    
//    test1D();
//    testQuad();
//    testTriangle();
//}
//
//// Integração 1D
//
//int TestTIntRule1d()
//{
//    cout <<__PRETTY_FUNCTION__<<endl;
//    int result = 0;
//    
//    cout << "Teste1: Criacao de um IntRule1d com ordem negativa" <<endl;
//    try {
//        IntRule1d Integral(-1);
//    } catch (std::bad_exception) {
//        result = 1;
//    }
//    
//    if (result == 0)
//        cout << "Teste1: Errado." <<endl;
//    else
//        cout << "Teste1: Ok." <<endl;
//    
//    result = 0;
//    
//    cout << "Teste2: Criacao de um IntRulde1d com ordem > 19" <<endl;
//    try {
//        IntRule1d  Integral(20);
//    } catch (std::bad_exception) {
//        result = 1;
//    }
//    if (result == 1)
//        cout << "Teste2: Errado." <<endl;
//    else
//        cout << "Teste2: Ok." <<endl;
//    
//    result = 0;
//    
//    cout << "Teste3: Verificacao se a ordem dada atribui corretamente o numero de pesos" <<endl;
//    IntRule1d  Integral(7);
//    IntRule1d  Integral2(18);
//    
//    if (Integral.NPoints() == 4)
//        result++;
//    else
//        result = 0;
//    
//    if (Integral2.NPoints() == 10)
//        result++;
//    else
//        result = 0;
//    
//    if (result == 2)
//        cout << "Teste3: Ok." <<endl;
//    else
//        cout << "Teste3: Errado." <<endl;
//    
//    result = 0;
//    
//    cout << "Teste4: Verificacao dos pesos e coordenadas para uma ordem x." <<endl;
//    VecDouble testcoord(1);
//    double peso;
//    int pontoteste = 3;
//    
//    Integral.Point(pontoteste, testcoord, peso);
//    
//    if (testcoord[0] == 0.33998104358485626 && peso == 0.65214515486254609)
//        result++;
//    else
//        result = 0;
//    
//    pontoteste = 8;
//    Integral2.Point(pontoteste,testcoord,peso);
//    
//    if (testcoord[0] == -0.14887433898163122 && peso == 0.29552422471475287)
//        result++;
//    else
//        result = 0;
//    
//    if (result != 2) {
//        cout << "Teste1: Errado." <<endl;
//    }
//    if (result == 2) {
//        cout << "Teste4: Ok." <<endl;
//        result = 1;
//    }
//    
//    cout << "---------------//----------------" <<endl;
//    
//    return result;
//}
//
//int TestPointFunctionFull()
//{
//    cout <<__PRETTY_FUNCTION__<<endl;
//    int result = 0;
//    
//    cout << "Teste1: Teste completo de pesos e coordenadas para uma ordem x." <<endl;
//    IntRule1d  Teste(6);
//    VecDouble weights(4,0);
//    VecDouble valcoord(4.0);
//    
//    valcoord[0] = -0.86113631159405257;  weights[0] = 0.34785484513745385;
//    valcoord[1] = 0.86113631159405257;   weights[1] = 0.34785484513745385;
//    valcoord[2] = -0.33998104358485626;  weights[2] = 0.65214515486254609;
//    valcoord[3] = 0.33998104358485626;   weights[3] = 0.65214515486254609;
//    
//    VecDouble coord(1,0.);
//    double peso;
//    
//    for(int i=0; i<Teste.NPoints() ;i++)
//    {
//        Teste.Point(i, coord, peso);
//        
//        if (coord[0] == valcoord[i] && peso == weights [i]) {
//            result++;
//        }
//    }
//    
//    
//    
//    if (result != 4) {
//        cout << "Teste1: Errado." <<endl;
//    }
//    if (result == 4) {
//        cout << "Teste1: Ok." <<endl;
//        result = 1;
//    }
//    
//    cout << "---------------//----------------" <<endl;
//    return result;
//}
//
//double Poli1D(double coord, int polorder)
//{
//    double y;
//    
//    y = pow(coord, polorder);
//    
//    return y;
//}
//
//int TestPoliIntegration()
//{
//    cout <<__PRETTY_FUNCTION__<<endl;
//    int result = 1;
//    
//    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
//    
//    for (int order = 0; order <20; order++)
//    {
//        IntRule1d TestInt(order);
//        double val = 0;
//        double weight;
//        VecDouble Coord(1);
//        for (int polorder = 0; polorder <= order; polorder++)
//        {
//            val = 0.;
//            for (int i=0; i<TestInt.NPoints(); i++) {
//                TestInt.Point(i, Coord, weight);
//                val = val + weight*Poli1D(Coord[0], polorder);
//            }
//            double exact = 1./(polorder+1.)*(1-pow(-1, polorder+1));
//            double error = val-exact;
//            
//            if (fabs(error) >= 1.e-15) {
//                std::cout << "A regra de integracao de ordem " << order << " nao integrou um polinome de ordem " << polorder << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
//                result = 0;
//            }
//        }
//        
//    }
//    if (result == 1) {
//        cout << "Teste1: Ok." <<endl;
//        result = 1;
//    }
//    else
//        cout << "Teste1: Errado." <<endl;
//    
//    cout << "---------------//----------------" <<endl;
//    return result;
//}
//
//void test1D()
//{
//    int res = 0;
//    int count = 0;
//    
//    res = TestTIntRule1d();
//    if (res == 1) {
//        count++;
//    }
//    
//    res = TestPointFunctionFull();
//    if (res == 1) {
//        count++;
//    }
//    
//    res = TestPoliIntegration();
//    if (res == 1) {
//        count++;
//    }
//    
//    if (count == 0) {
//        cout <<"\nSeu codigo passou no teste de Integracao para funcoes de uma dimensao. (:D)" << "\n---------------//----------------";
//    }
//    
//}
//
//// Integração 2D
//
//double Poli2D(VecDouble &coord, int orderx, int ordery)
//{
//    return pow(coord[0],orderx) * pow(1.-coord[1],ordery);
//}
//
//double Poli2DQuad(VecDouble &coord, int orderx, int ordery)
//{
//    return pow(coord[0],orderx) * pow(coord[1],ordery);
//}
//
//int TestPoliIntegrationQuad()
//{
//    cout <<__PRETTY_FUNCTION__<<endl;
//    int result = 1;
//    
//    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
//    
//    for (int order = 0; order <20; order++)
//    {
//        IntRuleQuad TestInt(order);
//        
//        //TestInt.Print(std::cout);
//        
//        double val = 0;
//        double weight;
//        VecDouble Coord(2);
//        for (int polorderx = 0; polorderx <= order; polorderx++)
//        {
//            for (int polordery=0; polordery <= order; polordery++)
//            {
//                val = 0.;
//                int NPoint = TestInt.NPoints();
//                for (int i=0; i<NPoint; i++) {
//                    TestInt.Point(i, Coord, weight);
//                    val = val + weight*Poli2DQuad(Coord, polorderx,polordery);
//                }
//                double exact = 1./(polorderx+1.)*(1-pow(-1, polorderx+1))*1./(polordery+1.)*(1-pow(-1, polordery+1));
//                double error = val-exact;
//                
//                if (fabs(error) >= 1.e-14) {
//                    std::cout << "A regra de integracao quadrilatero de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
//                    result = 0;
//                }
//            }
//        }
//    }
//    if (result == 1) {
//        cout << "TesteQuad1: Ok." <<endl;
//        result = 1;
//    }
//    else
//        cout << "TesteQuad1: Errado." <<endl;
//    
//    cout << "---------------//----------------" <<endl;
//    return result;
//}
//
//int TestPoliIntegrationTriangle()
//{
//    cout <<__PRETTY_FUNCTION__<<endl;
//    int result = 1;
//    
//    cout << "Teste1: Integral de um polinomio e verifica resultado." <<endl;
//    
//    for (int order = 0; order <20; order++)
//    {
//        IntRuleTriangle TestInt(order);
//        double val = 0;
//        double weight;
//        VecDouble Coord(2);
//        for (int polorderx = 0; polorderx <= order; polorderx++)
//        {
//            for (int polordery=0; polordery <= order-polorderx; polordery++)
//            {
//                val = 0.;
//                int npoints = TestInt.NPoints();
//                for (int i=0; i<npoints; i++) {
//                    TestInt.Point(i, Coord, weight);
//                    val = val + weight*Poli2D(Coord, polorderx,polordery);
//                }
//                double exact = 1./(1.+polorderx)/(2.+polorderx+polordery);
//                double error = val-exact;
//                
//                if (fabs(error) >= 1.e-15) {
//                    std::cout << "A regra de integracao triangular de ordem " << order << " nao integrou um polinome de ordem " << polorderx << " , " << polordery << " exato " << exact << " calculado " << val << " erro " << error << std::endl;
//                    result = 0;
//                }
//            }
//        }
//    }
//    if (result == 1) {
//        cout << "TesteTriangle1: Ok." <<endl;
//        result = 1;
//    }
//    else
//        cout << "TesteTriangle1: Errado." <<endl;
//    
//    cout << "---------------//----------------" <<endl;
//    return result;
//}
//
//void testQuad()
//{
//    int result = 0;
//    result = TestPoliIntegrationQuad();
//    if (result) {
//        cout <<"\nSeu codigo passou no teste de Integracao sobre quadrilatero. (:D)" << "\n---------------//----------------";
//        
//    }
//}
//
//void testTriangle()
//{
//    int result = 0;
//    result = TestPoliIntegrationTriangle();
//    if (result) {
//        cout <<"\nSeu codigo passou no teste de Integracao sobre triangulo. (:D)" << "\n---------------//----------------";
//        
//    }
//}
//
//// ----------------- TESTES GMESH E CMESH --------------------
//
//void GmeshCmeshTests(){
//    GeoMesh * gmesh = new GeoMesh();
//    std::string filename("malhateste.msh");
//    ReadGmsh * teste = new ReadGmsh();
//    
//    // Teste do BuildConnectivity
//    teste->Read(*gmesh, filename);
//    gmesh->BuildConnectivity();
//    gmesh->Print(std::cout);
//    
//    
//    // Teste do Paraview
//    VTKGeoMesh * testevtk = new VTKGeoMesh();
//    testevtk->PrintGMeshVTK(gmesh, "teste.vtk");
//    
//    // Teste do AutoBuild do cmesh:
//    CompMesh * cmesh = new CompMesh(gmesh);
//    cmesh->AutoBuild();
//}
//
//// --------------- TESTE FUNÇÃO P CALCSTIFF ------------------
//
//void StiffnessMatrix1Test(){
//    
//    GeoMesh *gmesh = new GeoMesh();
//    VecDouble coord(3,0.);
//    double length = 1;
//    int64_t nelem = 1;
//    
//    gmesh->SetNumElements(nelem);
//    gmesh->SetDimension(2);
//    gmesh->SetNumNodes(4);
//    
//    VecInt nodeindices(4,0.);
//    
//    // Determinando coordenadas
//    for (int i=0; i<2; i++) {
//        for (int j=0; j<2; j++) {
//            int id = i + j*2;
//            coord[0] = 0 + i*length;
//            coord[1] = 0 + j*length;
//            gmesh->Node(id).SetCo(coord);
//        }
//    }
//    nodeindices[0] = 0;
//    nodeindices[1] = 1;
//    nodeindices[2] = 3;
//    nodeindices[3] = 2;
//    
//    int materialid = 1;
//    int index = 0;
//    GeoElement *gel = new GeoElementTemplate<GeomQuad>(nodeindices, materialid, gmesh, index);
//    gmesh->SetElement(0, gel);
//    
//    Matrix perm(2, 2, 1.);
//    perm(1,0) = 0.;
//    perm(0,1) = 0.;
//    Poisson * pos = new Poisson(1, perm);
//    pos->SetForceFunction(forcefunction);
//    
//    CompMesh * cmesh = new CompMesh(gmesh);
//    cmesh->SetDefaultOrder(1);
//    cmesh->SetNumberMath(1);
//    std::vector<MathStatement *> mathvec(1);
//    mathvec[0] = pos;
//    cmesh->SetMathVec(mathvec);
//    cmesh->AutoBuild();
//    
//    CompElement * cel = cmesh->GetElement(0);
//    Matrix ek;
//    Matrix ef;
//    cel->CalcStiff(ek, ef);
//    ek.Print();
//}
//
//// ----------------- TESTE DA ASSEMBLAGEM --------------------
//
//void AssembleTest(){
//    
//    GeoMesh *gmesh = new GeoMesh();
//    VecDouble coord(3,0.);
//    double length = 1;
//    int64_t nelem = 2;
//    
//    gmesh->SetNumElements(nelem);
//    gmesh->SetDimension(2);
//    gmesh->SetNumNodes(6);
//    
//    // Determinando coordenadas
//    for (int i=0; i<nelem+1; i++) {
//        for (int j=0; j<nelem+1; j++) {
//            int id = i*2 + j;
//            coord[0] = length + i*length/nelem;
//            coord[1] = 0 + j*length/nelem;
//            gmesh->Node(id).SetCo(coord);
//        }
//    }
//    
//    VecInt nodeindices(4,0.);
//    nodeindices[0] = 0;
//    nodeindices[1] = 2;
//    nodeindices[2] = 3;
//    nodeindices[3] = 1;
//    int materialid = 1;
//    int index = 0;
//    GeoElement *gel = new GeoElementTemplate<GeomQuad>(nodeindices, materialid, gmesh, index);
//    gmesh->SetElement(0, gel);
//    
//    nodeindices[0] = 2;
//    nodeindices[1] = 4;
//    nodeindices[2] = 5;
//    nodeindices[3] = 3;
//    index = 1;
//    gel = new GeoElementTemplate<GeomQuad>(nodeindices, materialid, gmesh, index);
//    gmesh->SetElement(1, gel);
//    
//    gmesh->BuildConnectivity();
//    
//    gmesh->Print(cout);
//    
//    CompMesh *cmesh = new CompMesh(gmesh);
//    cmesh->SetDefaultOrder(1);
//    
//    const std::vector<MathStatement *> mathvec(2);
//    Matrix perm(2,2,1.);
//    perm(0,1) = 0; perm(1,0) = 0;
//    Poisson * pos = new Poisson(1, perm);
//    pos->SetForceFunction(forcefunction);
//    cmesh->SetNumberMath(2);
//    cmesh->SetMathStatement(0, pos);
//    cmesh->SetMathStatement(1, pos);
//    
//    cmesh->AutoBuild();
//    
//    Assemble * assem = new Assemble(cmesh);
//    Matrix globmat;
//    Matrix rhs;
//    assem->Compute(globmat, rhs);
//    
//    globmat.Print();
//}
//
//// -------------- TESTE COND CONTORNO E SOLUTION -------------
//
//void SolutionTest(){
//    
//    GeoMesh *gmesh = new GeoMesh();
//    gmesh->SetDimension(2);
//    
//    double len_x = 1.;
//    double len_y = 1.;
//    int pOrder = 2;
//    
//    VecDouble error0(3,0.);
//    VecDouble error(3,0.);
//    
//    for (int64_t nelem_i = 1; nelem_i<8; nelem_i++) {
//        
//        double len_ix = len_x/(nelem_i);
//        double len_iy = len_y/(nelem_i);
//        int64_t nnodes = nelem_i+1;
//        
//        int64_t nelem = nelem_i*nelem_i;
//        gmesh->SetNumElements(nelem+nelem_i*4);
//        gmesh->SetNumNodes(nnodes*nnodes);
//        
//        // Enumeration: vertical order - from the below to the top, and from the left to the right
//        VecDouble coord(3,0.);
//        int64_t id;
//        
//        for(int i = 0; i < nnodes; i++){
//            for(int j = 0; j < nnodes; j++){
//                id = i*nnodes + j;
//                coord[0] = (j)*len_ix;
//                coord[1] = (i)*len_iy;
//                gmesh->Node(id).SetCo(coord);
//            }
//        }
//        
//        int materialid = 1;
//        for(int i = 0; i < nelem_i; i++){
//            for(int j = 0; j < nelem_i; j++){
//                int64_t id = i*nelem_i+j;
//                
//                VecInt connect(4,0); // sentido anti-horario
//                connect[0] = j*nnodes + i;
//                connect[1] = connect[0]+1;
//                connect[3] = connect[0]+nnodes;
//                connect[2] = connect[3]+1;
//                
//                GeoElement *gel = new GeoElementTemplate<GeomQuad>(connect, materialid, gmesh, id);
//                gmesh->SetElement(id, gel);
//            }
//        }
//        
//        int bcleft = -1;
//        int64_t elempos = nelem;
//        VecInt connect(2,0);
//        for (int64_t iel=0; iel<nelem_i; iel++) {
//            connect[0] = connect[1];
//            connect[1] = connect[0]+nnodes;
//            GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcleft, gmesh, elempos);
//            gmesh->SetElement(elempos, gelbc0);
//            elempos++;
//        }
//        
//        int bcright = -2;
//        connect[1] = nnodes-1;
//        for (int64_t iel=0; iel<nelem_i; iel++) {
//            connect[0] = connect[1];
//            connect[1] = connect[0]+nnodes;
//            GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcright, gmesh, elempos);
//            gmesh->SetElement(elempos, gelbc0);
//            elempos++;
//        }
//        
//        int bcbottom = -3;
//        connect[1] = 0;
//        for (int64_t iel=0; iel<nelem_i; iel++) {
//            connect[0] = connect[1];
//            connect[1] = connect[0]+1;
//            GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bcbottom, gmesh, elempos);
//            gmesh->SetElement(elempos, gelbc0);
//            elempos++;
//        }
//        
//        int bctop = -4;
//        connect[1] = nnodes*nnodes-nnodes;
//        for (int64_t iel=0; iel<nelem_i; iel++) {
//            connect[0] = connect[1];
//            connect[1] = connect[0]+1;
//            GeoElement *gelbc0 = new GeoElementTemplate<Geom1d>(connect, bctop, gmesh, elempos);
//            gmesh->SetElement(elempos, gelbc0);
//            elempos++;
//        }
//        
//        gmesh->BuildConnectivity();
//        
//        VTKGeoMesh * testevtk = new VTKGeoMesh();
//        testevtk->PrintGMeshVTK(gmesh, "teste.vtk");
//        
//        // Material:
//        Matrix perm(2, 2, 1.);
//        perm(1,0) = 0.;
//        perm(0,1) = 0.;
//        Poisson * pos = new Poisson(materialid, perm);
//        pos->SetForceFunction(&forcefunction);
//        
//        // Condição de contorno:
//        Matrix projection(2,2,0.);
//        int bctype = 0;
//        Matrix Val1(2,1,0.);
//        Matrix Val2(2,1,0.);
//        
//        L2Projection * bc0 = new L2Projection(bctype, bcleft, projection, Val1, Val2);
//        bc0->SetExactSolution(&solexact);
//        
//        L2Projection * bc1 = new L2Projection(bctype, bcright, projection, Val1, Val2);
//        bc1->SetExactSolution(&solexact);
//        
//        L2Projection * bc2 = new L2Projection(bctype, bcbottom, projection, Val1, Val2);
//        bc2->SetExactSolution(&solexact);
//        
//        L2Projection * bc3 = new L2Projection(bctype, bctop, projection, Val1, Val2);
//        bc3->SetExactSolution(&solexact);
//        
//        // Criando malha computacional:
//        CompMesh * cmesh = new CompMesh(gmesh);
//        cmesh->SetDefaultOrder(pOrder);
//        cmesh->SetNumberMath(5);
//        cmesh->SetMathStatement(0, pos);
//        cmesh->SetMathStatement(1, bc0);
//        cmesh->SetMathStatement(2, bc1);
//        cmesh->SetMathStatement(3, bc2);
//        cmesh->SetMathStatement(4, bc3);
//        cmesh->AutoBuild();
//        
//        // Solução:
//        Analysis * an = new Analysis(cmesh);
//        an->RunSimulation();
//        
//        VecDouble sol = cmesh->Solution();
//        
//        //        std::cout << "Elementos: " << nelem_i << std::endl;
//        //        for (int64_t i=0; i<sol.size(); i++) {
//        //            std::cout << sol[i] << std::endl;
//        //        }
//        //        std::cout << " * * * * * * * * * * * * * * * " << std::endl;
//        
//        PostProcessTemplate<Poisson> defPostProc = PostProcessTemplate<Poisson>(an);
//        defPostProc.SetExact(&solexact);
//        //        defPostProc.AppendVariable(Poisson::ESol);
//        //        defPostProc.AppendVariable(Poisson::EDSol);
//        VecDouble error = an->PostProcessError(std::cout, defPostProc);
//        
//        std::string filename("solution.vtk");
//        testevtk->PrintCMeshVTK(cmesh, 2, filename);
//        
//        std::cout << (log(error[0])-log(error0[0])) / (log(len_x/nelem_i)-log(len_x/(nelem_i-1))) << std::endl;
//        std::cout << (log(error[1])-log(error0[1])) / (log(len_x/nelem_i)-log(len_x/(nelem_i-1))) << std::endl;
//        std::cout << (log(error[2])-log(error0[2])) / (log(len_x/nelem_i)-log(len_x/(nelem_i-1))) << std::endl;
//        error0 = error;
//        
//    }
//    
//}
//
//void forcefunction(const VecDouble &co, VecDouble &result)
//{
//    result.resize(1);
//    
//    double xv = co[0];
//    double yv = co[1];
//    
//    double Pi = M_PI;
//    result[0] = 8.0*Pi*Pi*cos(2.0*Pi*yv)*sin(2.0*Pi*xv);
//    
//    //    result[0] = -2*sin(xv)*cos(yv);
//}
//
//void solexact(const VecDouble &co, VecDouble &result, Matrix &deriv)
//{
//    deriv.Resize(2,1);
//    result.resize(1);
//    
//    double xv = co[0];
//    double yv = co[1];
//    
//    double Pi = M_PI;
//    
//    double v_x = cos(2*Pi*yv)*sin(2*Pi*xv);
//    
//    result[0] = v_x;
//    
//    deriv(0,0)=  2*Pi*cos(2*Pi*yv)*cos(2*Pi*xv);
//    deriv(1,0)= -2*Pi*sin(2*Pi*xv)*sin(2*Pi*yv);
//    
//    //    double v_x = sin(xv)*cos(yv);
//    //
//    //    result[0] = v_x;
//    //
//    //    deriv(0,0)=  cos(xv)*cos(yv);
//    //    deriv(1,0)= -sin(xv)*sin(yv);
//}
