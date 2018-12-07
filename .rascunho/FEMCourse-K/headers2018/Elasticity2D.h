//
//  Header.h
//  FemSC
//
//  Created by Karolinne Oliveira Coelho on 7/5/18.
//
//

#ifndef Elasticity2D_h
#define Elasticity2D_h

#include "MathStatement.h"
#include "DataTypes.h"
#include "IntPointData.h"
#include <functional>

class Elasticity2D : public MathStatement
{
    // Force function
    VecDouble  ff;
    
    // Elasticity
    double fE;
    
    // Poisson coefficient
    double fnu;
    
    // Lame coefficients
    double flambda;
    double fmu;
    
    /** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
    int fPlaneStress;
    
    // Force funtion
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
    std::function<void(const VecDouble &loc, VecDouble &result, Matrix &deriv)> SolutionExact;
    
public:
    
    enum PostProcVar {ENone, ESol, ESigmaX, ESigmaY, EEpsilonX, EEpsilonY, ESolExact, EDSolExact};
    
    // Default constructor of Poisson
    Elasticity2D();
    
    // Constructor of Poisson
    Elasticity2D(int materialid, double fE, double fnu);
    
    // Copy constructor of Poisson
    Elasticity2D(const Elasticity2D &copy);
    
    // Operator of copy
    Elasticity2D &operator=(const Elasticity2D &copy);
    
    // Method for creating a copy of the element
    virtual Elasticity2D *Clone() const;
    
    // Destructor of Poisson
    virtual ~Elasticity2D();
    
    // Return the permeability matrix
    double GetYoungModulus() const;
    
    // Set the permeability matrix
    void SetYoungModulus(const double &elast);
    
    // Return the permeability matrix
    double GetPoisson() const;
    
    // Set the permeability matrix
    void SetPoisson(const double &nu);
    
    // Return the force function related to Poisson math statement
    std::function<void(const VecDouble &co, VecDouble &result)> GetForceFunction() const
    {
        return forceFunction;
    }
    
    // Set the force function related to Poisson math statement
    void SetForceFunction(const std::function<void(const VecDouble &co, VecDouble &result)> &f)
    {
        forceFunction = f;
    }
    
    // Set the exact solution related to Poisson math statement
    void SetExactSolution(const std::function<void(const VecDouble &loc, VecDouble &result, Matrix &deriv)> &Exact)
    {
        SolutionExact = Exact;
    }
    
    // returns the integrable dimension of the material
    int Dimension() const {return 2;}
    
    virtual int NEvalErrors() const;
    
    // Return the number of state variables
    virtual int NState() const;
    
    virtual int VariableIndex(const PostProcVar var) const;
    
    // Return the variable index associated with the name
    virtual PostProcVar VariableIndex(const std::string &name);
    
    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    virtual int NSolutionVariables(const PostProcVar var);
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const;
    
    // Method to implement error over element's volume
    virtual void ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const;
    
    // Prepare and print post processing data
    virtual void PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const;
    
    virtual double Inner(Matrix &S,Matrix &T) const;
    
    void SetPlaneStressState(int state);
    
};
#endif /* Elasticity2D_h */
