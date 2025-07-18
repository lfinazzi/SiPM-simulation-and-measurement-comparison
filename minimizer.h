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
    Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int iterations, double derivStepPerc, double startingLearningRate);

    /****************************************************************************
     * Calculates the negative log likelihood of measuring the sim data
     * given the measurement data (which is treated as the real distribution).
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateNLL();

    /****************************************************************************
     * Calculates the negative log likelihood of measuring the sim data
     * given the measurement data (which is treated as the real distribution)
     * using an external histogram.
     * Arguments:
     *      TH1D : external histogram
     *      std::string : option to use: plus or minus
     * Returns:
     *      None
     ****************************************************************************/
    void CalculateNLLExternal(TH1D hist, std::string option);

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
     * uses gradientIters iterations to average the gradient
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

    /****************************************************************************
     * Normalizes parameters so all of them are in the range [0, 1] for 
     * minimization purposes.
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void NormalizeParameters();

    /****************************************************************************
     * Denormalizes parameters after minimzation to rerun simulations.
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void DenormalizeParameters();

    /****************************************************************************
     * Prints variable parameters to console
     * Arguments:
     *      None
     * Returns:
     *      None
     ****************************************************************************/
    void PrintParameters();
    
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

    std::vector<bool> vparams_update = {true, true, true, true, true, true, true, true, true, true, true, true, true};

    std::vector<double> svector;
    std::vector<double> schangevector;

    // aux minimization variables
    double NLL;
    std::vector<double> gradient;
    double aux_NLL_plus;
    double aux_NLL_minus;
    VariableParameters aux_vparams;
    std::vector<double> aux_vparams_vector;
    TH1D aux_simHist;

    // minimization parameters
    int iterations;
    double derivStepPerc;
    int currentStep;
    int gradientIters = 1;

    // --- Adam optimizer state ---
    std::vector<double> m;  // First moment vector
    std::vector<double> v;  // Second moment vector
    int t = 0;              // Time step

    // Adam hyperparameters
    double beta1 = 0.7;         // standard is 0.9, but lower is better for noisier gradients
    double beta2 = 0.999;
    double epsilon = 1e-8;
    double startingLearningRate;

    // hist parameters
    double maxHistVal = 2.75E-12;
    int histBins = 400;
};