#ifndef PARAMETERS_H
    #define PARAMETERS_H
    #include "parameters.h"
#endif

#ifndef REALDATA_H
    #define REALDATA_H
    #include "realdata.h"
#endif

#ifndef SIMULATION_H
    #define SIMULATION_H
    #include "simulation.h"
#endif

#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif

#include <vector>
#include <string>

class Minimizer
{

public: 
      
    // constructor . Loads all parameter and data and initializes simulation and histograms
    Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int iterations, double stepSize, double delta);

    /****************************************************************************
     * Calculates S = sum( (data_i - simulation_i )^2 ) using their histograms
     * Arguments:
     *      bool : calculates S for aux simulation if true
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateS(bool aux = false);

    /****************************************************************************
     * Calculates the derivative  of S with respect to a variable parameters
     * Arguments:
     *      int : parameter index
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateDerivative(int numParam);

    /****************************************************************************
     * Calculates the gradient of S with respect to the variable parameters
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateGradient();
    
    // returns data member
    inline RealData GetData() { return data; }

    // returns current simulation member
    inline Simulation GetSim() { return sim; }

    inline TH1D* GetDataHist() { return &dataHist; }
    inline TH1D* GetSimHist() { return &simHist; }


private:
    Simulation sim;
    RealData data;
    TH1D simHist;
    TH1D dataHist;
    FixedParameters fparams;
    VariableParameters vparams;
    std::vector<double> vparams_vector;

    // aux minimization variables
    double S;
    std::vector<double> gradient;
    double aux_S;
    VariableParameters aux_vparams;
    std::vector<double> aux_vparams_vector;
    TH1D aux_simHist;

    // minimization parameters
    int iterations;
    double stepSize;
    int currentStep;
    double delta;

};