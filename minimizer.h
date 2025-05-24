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
      
    Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data)
        : sim(_fparameters, _vparameters), data(_data), fparams(_fparameters), vparams(_vparameters) 
    {
        sim.Simulate();
        vparams_vector.resize(16);
        ParameterToVector(vparams, vparams_vector);
    }  
    
    // returns data member
    inline RealData GetData() { return data; }

    // returns current simulation member
    inline Simulation GetSim() { return sim; }


private:
    Simulation sim;
    RealData data;
    FixedParameters fparams;
    VariableParameters vparams;
    std::vector<double> vparams_vector;

};