//
//  ShapeQuad.cpp
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include "DataTypes.h"
#include "tpanic.h"
#include "Shape1d.h"
#include "ShapeQuad.h"

/// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
    
    if (orders[0]>2) {
        std::cout << "ShapeQuad::Shape, only implemented until order = 2" << std::endl;
        DebugStop();
    }
    
    int nf = phi.size();
    
    VecDouble phi_xi(nf), phi_eta(nf);
    Matrix dphi_xi(2,nf), dphi_eta(2,nf);

//    VecInt orders_aux(3, 1);
//    for (int isides=0; isides<3; isides++) {
//        orders_aux[isides] = orders[isides];
//    }
//    int indexes[3][3] = {{0, 7, 3},{4, 8, 6},{1, 5, 2}};
//    int id = (nf%2 == 0) ? 1 : 0;
//    
//    VecDouble etav(1, xi[1]);
//    
//    for (int kxi=0; kxi < 2; kxi++) {
//        Shape1d::Shape(xi, orders_aux, phi_xi, dphi_xi);
//        
//        for (int eta=0; eta < 2; eta++) {
//            Shape1d::Shape(etav, orders_aux, phi_eta, dphi_eta);
//            
//            int index = indexes[(1+id*kxi)*kxi][(1+id*eta)*eta];
//            
//            phi[index] = phi_xi[kxi]*phi_eta[eta];
//            dphi(0,index) = dphi_xi(0,kxi)*phi_eta[eta];
//            dphi(1,index) = dphi_eta(0,eta)*phi_xi[kxi];
//        }
//    }
    double x[2],dx[2],y[2],dy[2];
    x[0]  =  (1.-xi[0])/2.;
    x[1]  =  (1.+xi[0])/2.;
    dx[0] = -0.5;
    dx[1] =  0.5;
    y[0]  =  (1.-xi[1])/2.;
    y[1]  =  (1.+xi[1])/2.;
    dy[0] = -0.5;
    dy[1] =  0.5;
    phi[0]  = x[0]*y[0];
    phi[1]  = x[1]*y[0];
    phi[2]  = x[1]*y[1];
    phi[3]  = x[0]*y[1];
    dphi(0,0) = dx[0]*y[0];
    dphi(1,0) = x[0]*dy[0];
    dphi(0,1) = dx[1]*y[0];
    dphi(1,1) = x[1]*dy[0];
    dphi(0,2) = dx[1]*y[1];
    dphi(1,2) = x[1]*dy[1];
    dphi(0,3) = dx[0]*y[1];
    dphi(1,3) = x[0]*dy[1];

    if (orders[0]==2) {
        int is;
        for(is=4; is<8; is++)
        {
            phi[is] = phi[is%4]*phi[(is+1)%4];
            dphi(0,is) = dphi(0,is%4)*phi[(is+1)%4]+phi[is%4]*dphi(0,(is+1)%4);
            dphi(1,is) = dphi(1,is%4)*phi[(is+1)%4]+phi[is%4]*dphi(1,(is+1)%4);
        }
        phi[8] = phi[0]*phi[2];
        dphi(0,8) = dphi(0,0)*phi[2]+phi[0]*dphi(0,2);
        dphi(1,8) = dphi(1,0)*phi[2]+phi[0]*dphi(1,2);
        
        // Make the generating shape functions linear and unitary
        for(is=4; is<8; is++)
        {
            phi[is] += phi[8];
            dphi(0,is) += dphi(0,8);
            dphi(1,is) += dphi(1,8);
            phi[is] *= 4.;
            dphi(0,is) *= 4.;
            dphi(1,is) *= 4.;
        }
        phi[8] *= 16.;
        dphi(0,8) *= 16.;
        dphi(1,8) *= 16.;
    }
    
}

/// returns the number of shape functions associated with a side
int ShapeQuad::NShapeFunctions(int side, int order){
    if(side<4)
        return 1;//0 a 4
    else if(side<8)
        return (order-1);//6 a 14
    else if(side==8)
        return ((order-1)*(order-1));
    
    std::cout << "ShapeQuad::NShapeFunctions, bad parameter side " << side << std::endl;
    DebugStop();
    
    return 0;
}

/// returns the total number of shape functions
int ShapeQuad::NShapeFunctions(VecInt &orders){
    
    int res=4;
    for(int in=4; in<orders.size(); in++) {
        res += NShapeFunctions(in, orders[in]);
    }
    
    return res;
}
