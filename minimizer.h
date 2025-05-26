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
    Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int iterations, double stepPerc, double minmimUpdatePerc);

    /****************************************************************************
     * Calculates S = sum( (data_i - simulation_i )^2 ) using their histograms
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateS();

    /****************************************************************************
     * Calculates S = sum( (data_i - simulation_i )^2 ) using an external histogram
     * Arguments:
     *      TH1D : external histogram
     *      std::string : option to use: plus or minus
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateSExternal(TH1D hist, std::string option);

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

    /****************************************************************************
     * Updates the variable parameters and the simulation with the new parameters
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void Update();

    /****************************************************************************
     * Runs minimizer for a number of iterations
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void RunMinimizer();
    
    // returns data member
    inline RealData GetData() { return data; }

    // returns current simulation member
    inline Simulation GetSim() { return sim; }

    // returns data histogra
    inline TH1D* GetDataHist() { return &dataHist; }

    // returns simulation histogram
    inline TH1D* GetSimHist() { return &simHist; }

    // returns S value for each iteration
    inline std::vector<double> GetSVector() { return svector; }

    // returns S value change for each iteration
    inline std::vector<double> GetSChangeVector() { return schangevector; }


private:
    Simulation sim;
    RealData data;
    TH1D simHist;
    TH1D dataHist;
    FixedParameters fparams;
    VariableParameters vparams;
    std::vector<double> vparams_vector;
    std::vector<double> svector;
    std::vector<double> schangevector;

    // aux minimization variables
    double S;
    std::vector<double> gradient;
    double aux_S_plus;
    double aux_S_minus;
    VariableParameters aux_vparams;
    std::vector<double> aux_vparams_vector;
    TH1D aux_simHist;

    // minimization parameters
    int iterations;
    double stepPerc;
    int currentStep;

    // --- Adam optimizer state ---
    std::vector<double> m;  // First moment vector
    std::vector<double> v;  // Second moment vector
    int t = 0;              // Time step

    // Adam hyperparameters
    double beta1 = 0.7;         // standard is 0.9, but lower is better for noisier gradients
    double beta2 = 0.999;
    double epsilon = 1e-8;
    double minimUpdatePerc;
};