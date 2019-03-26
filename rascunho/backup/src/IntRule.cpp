/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 * 
 * Pedro Lima 2018
 */

#include "IntRule.h"
#include "tpanic.h"
#include "DataTypes.h"



IntRule::IntRule() : fPoints(1, 1), fWeights(0)
{
    fOrder = 0;
}

IntRule::IntRule(int order)
{
    fOrder = order;
}

IntRule::~IntRule(){}

IntRule &IntRule::operator=(const IntRule &copy)
{
    fOrder = copy.fOrder;
    fPoints = copy.fPoints;
    fWeights = copy.fWeights;
    return *this;
}

IntRule::IntRule(const IntRule& copy)
{
    fOrder = copy.fOrder;
    fPoints = copy.fPoints;
    fWeights = copy.fWeights;
}

int IntRule::NPoints() const
{
    return fWeights.size();
}

void IntRule::Print(std::ostream &out)
{
    int npts = fWeights.size();
    out << "Gauss-Legendre \n";
    out << "Order: " << fOrder << " Number of points: " << npts << std::endl;

    double w;
    VecDouble pts(3, 0.);
    for (int i = 0; i < npts; i++)
    {
        Point(i, pts, w);
        out << "ip " << i << "; pos (" << pts[0] << ", " << pts[1] << "); w = " << w << std::endl;
    }
}

void IntRule::Point(int p, VecDouble &co, double &w) const
{
    for (int j = 0; j < fPoints.Cols(); j++)
    {
        co[j] = fPoints.GetVal(j, 0);
    }
    //w = fWeights(j);
}
