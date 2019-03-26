//
//  ShapeTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "ShapeTriangle.h"
#include "DataTypes.h"
#include "tpanic.h"
#include "Shape1d.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
    
    if (orders[6]<1) {
        std::cout << "ShapeTriangle::Shape, bad parameter order" << std::endl;
        DebugStop();
    }
    if (orders[6]>2) {
        std::cout << "ShapeTriangle::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    
//    int nf  = NShapeFunctions(orders);
//    phi.resize(nf);
//    dphi.Resize(2, nf);
//    
//    phi[0]=1-xi[0]-xi[1];
//    phi[1]=xi[0];
//    phi[2]=xi[1];
//    
//    dphi(0,0)=-1;
//    dphi(0,1)=1;
//    dphi(0,2)=0;
//    dphi(1,0)=-1;
//    dphi(1,1)=0;
//    dphi(1,2)=1;
//    
//    if (orders[6]==2) {
//        
//        phi[0]=2*phi[0]*(phi[0]-0.5);
//        phi[1]=2*phi[1]*(phi[1]-0.5);
//        phi[2]=2*phi[2]*(phi[2]-0.5);
//        phi[3]=4*phi[0]*phi[1];
//        phi[4]=4*phi[1]*phi[2];
//        phi[5]=4*phi[2]*phi[0];
//        
//        dphi(0,0) = -2*(0.5-xi[1]-xi[0]) - 2*(1-xi[1]-xi[0]);
//        dphi(0,1) = 2*(-0.5+xi[0])+2*xi[0];
//        dphi(0,3) = 4*(1-xi[1]-xi[0]) - 4*xi[0];
//        dphi(0,4) = 4*xi[1];
//        dphi(0,5) = -4*xi[1];
//        
//        dphi(1,0) = -2*(0.5-xi[1]-xi[0]) - 2*(1-xi[1]-xi[0]);
//        dphi(1,2) = 2*(-0.5+ xi[1]) + 2*xi[1];
//        dphi(1,3) = -4*xi[0];
//        dphi(1,4) = 4*xi[0];
//        dphi(1,5) = -4*xi[1] + 4*(1-xi[1]-xi[0]);
//        
//    }
    
    phi[0] =  1.-xi[0]-xi[1];
    phi[1] =  xi[0];
    phi[2] =  xi[1];
    dphi(0,0) = -1.;
    dphi(1,0) = -1.;
    dphi(0,1) =  1.;
    dphi(1,1) =  0.;
    dphi(0,2) =  0.;
    dphi(1,2) =  1.;
    
    if (orders[6] == 2) {
        int is;
        for(is=3; is<6; is++)
        {
            int is1 = is%3;
            int is2 = (is+1)%3;
            phi[is] = phi[is1]*phi[is2];
            dphi(0,is) = dphi(0,is1)*phi[is2]+phi[is1]*dphi(0,is2);
            dphi(1,is) = dphi(1,is1)*phi[is2]+phi[is1]*dphi(1,is2);
        }
        int is1 = 0;
        int is2 = 1;
        int is3 = 2;
        phi[is] = phi[is1]*phi[is2]*phi[is3];
        dphi(0,is) = dphi(0,is1)*phi[is2]*phi[is3]+phi[is1]*dphi(0,is2)*phi[is3]+phi[is1]*phi[is2]*dphi(0,is3);
        dphi(1,is) = dphi(1,is1)*phi[is2]*phi[is3]+phi[is1]*dphi(1,is2)*phi[is3]+phi[is1]*phi[is2]*dphi(1,is3);
        
        // Make the generating shape functions linear and unitary
        double mult[] = {1.,1.,1.,4.,4.,4.,27.};
        for(is=3;is<6; is++)
        {
            phi[is] *= mult[is];
            dphi(0,is) *= mult[is];
            dphi(1,is) *= mult[is];
        }
    }
    
    
}

/// returns the number of shape functions associated with a side
int ShapeTriangle::NShapeFunctions(int side, int order){
    switch(side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
        case 6:
            return order-1;
    }
    
    DebugStop();
    std::cout << "ShapeTriangle::NShapeFunctions, bad parameter side " << std::endl;
    return 0;
}

/// returns the total number of shape functions
int ShapeTriangle::NShapeFunctions(VecInt &orders){
    
    int res=3;
    for(int in=3; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
    
}
