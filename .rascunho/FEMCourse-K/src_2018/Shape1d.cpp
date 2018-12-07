//
//  Shape1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include <cmath>
#include <math.h>
#include "tpanic.h"
#include "DataTypes.h"
#include "Shape1d.h"

void Shape1d::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi){
    
    if (orders[0] < 0 || orders[1] < 0) {
        std::cout << "Shape1d::Shape: Invalid dimension for arguments: order" << std::endl;
        DebugStop();
    }
    phi.resize(orders[2]+1);
    dphi.Resize(1,orders[2]+1);
    
    if ( orders[0] == 0)
    {
        phi[0] = 1.;
        dphi(0,0) = 0.;
    }
    
//    int order = orders[0]-1;
//    for (int i=0; i<order+1; i++) {
//        double epsi=-1.+i*2./order;
//
//        for (int j=0; j<order+1; j++) {
//
//            Matrix axdphi(1,order+1);
//            if (i!=j) {
//                double epsj=-1.+j*2./order;
//                phi[i]*=(xi[0]-epsj)/(epsi-epsj);
//                axdphi(0,i)=1/(epsi-epsj);
//
//                for (int k=0; k<order+1; k++) {
//                    if (k!=i&&k!=j) {
//                        epsj=-1.+k*2./order;
//                        axdphi(0,i)*=(xi[0]-epsj)/(epsi-epsj);
//                    }
//                }
//                dphi(0,i)+=axdphi(0,i);
//            }
//        }
//    }
    
    phi[0] = (1-xi[0])/2.;
    phi[1] = (1+xi[0])/2.;
    dphi(0,0) = -0.5;
    dphi(0,1)= 0.5;
    if (orders[0] == 2) {
        phi[2]= phi[0]*phi[1];
        dphi(0,2) = dphi(0,0)*phi[1]+phi[0]*dphi(0,1);
        
        phi[2] *= 4.;
        dphi(0,2) *= 4.;
    }
    
    
}

/// returns the number of shape functions associated with a side
int Shape1d::NShapeFunctions(int side, int order){
    int nsf = 0;
    
    if(side<2){
        nsf += 1;
    }
    else if(side<3){
        nsf += order-1;
    }
    
    return nsf;
}

/// returns the total number of shape functions
int Shape1d::NShapeFunctions(VecInt &orders) {
    
    int nsf_tot = 0;
    for (int is=0; is<3; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}
