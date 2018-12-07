//
//  TIntRule1d.cpp
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 4/3/18.
//
//

#include <cmath>
#include "DataTypes.h"
#include "IntRule.h"
#include "tpanic.h"

IntRule::IntRule() : fOrder(0), fPoints(), fWeights(0){
    
}

IntRule::IntRule(int order){
    fPoints.Resize(0, 0);
    fWeights.resize(0);
    SetOrder(order);
};

IntRule::~IntRule(){
    
};

IntRule &IntRule::operator=(const IntRule &copy){
    fOrder = copy.fOrder;
    fPoints = copy.fPoints;
    fWeights = copy.fWeights;
    return *this;
};

IntRule::IntRule(const IntRule &copy){
    fOrder = copy.fOrder;
    fPoints = copy.fPoints;
    fWeights = copy.fWeights;
};

int IntRule::NPoints() const {
    return fWeights.size();
};

void IntRule::Point(int p, VecDouble &co, double &weight) const {
    
    if (p > NPoints()) {
        std::cout << "IntRule:: Point id higher than number of integration points" << std::endl;
        DebugStop();
    }
    if(co.size() > 3 )
    {
        std::cout << "IntRule:: Point coordinates must have 3 values" << std::endl;
        DebugStop();
    }
    if (fOrder < 0) {
        std::cout << "IntRule:: NULL integration order\n";
        DebugStop();
    }
    
    for (int idim=0; idim<fPoints.Cols(); idim++)
        co[idim] = fPoints.GetVal(p, idim);
    
    weight=fWeights[p];
};

void IntRule::Print(std::ostream &out) const {
    
    int npts = fWeights.size();
    out << "Cubature rule Gauss-Legendre " << std::endl;
    out << "Order: " << fOrder << " Number of points: " << npts << std::endl;
    
    double w;
    VecDouble pts(3,0.);
    for(int ipts=0; ipts<npts; ipts++) {
        Point(ipts, pts, w);
        out << "ip " << ipts << " pos " << pts[0] << ", " << pts[1] << " w " << w << std::endl;
    }
}
