//
//  IntRule.h
//  FEMCoursePZ
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FEMCoursePZ__IntRule__
#define __FEMCoursePZ__IntRule__

#include <cmath>
#include <stdio.h>
#include "DataTypes.h"

class IntRule
{
  
protected:
    
    // Polynomial order of the integration rule
    int fOrder;
    
    // Number of integration points for this object
    Matrix fPoints;
    
    // Weight of the integration point
    VecDouble fWeights;

public:
  
    // Default Constructor of integration rule
    IntRule(){};
    
    // Constructor of integration rule
    IntRule(int order){};
    
    // Destructor of integration rule
    ~IntRule();
    
    // Operator of copy
    virtual IntRule &operator=(const IntRule &copy);
    
    // Copy constructor of integration rule
    IntRule(const IntRule &copy);
    
    // Method to set polynomial order of the integration rule
    virtual void SetOrder(int order)
    {
        fOrder=order;
    }
    
    // Method to get polynomial order of the integration rule
    virtual int GetOrder()
    {
        return fOrder;
    }
    
    // Method to return the number of integration points
    virtual int NPoints() const;
    
    // Function returning coordinates and weights of integration points
    virtual void Point(int p, VecDouble &co, double &weight) const;
    
    // Function for printing results
    virtual void Print(std::ostream &out) const;
    
};


#endif /* defined(__FEMCoursePZ__TIntRule1d__) */