#include "TH1D.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#define EVENT_NUMBER_IN_BYTES 16

#ifndef PARAMETERS_H
    #define PARAMETERS_H
    #include "parameters.h"
#endif

class RealData{

public:
    // constructor --> initialized with path of .bin file
    RealData(std::string pathToFile, int eventsToLoad);


    /****************************************************************************
     * Loads SiPM data from binary files in ./data/ folder
     * Arguments:
     *      std::string : path to .bin file
     *      int : number of events to load
     * Returns:
     *      None
     ****************************************************************************/
    void LoadData(std::string path, int eventsToLoad);


    /****************************************************************************
     * Returns chargeValues member
     * Arguments:
     *      None
     * Returns:
     *      std::vector<double> : charge values
     ****************************************************************************/
    inline std::vector<double> GetChargeData() { return chargeValues; }


    /****************************************************************************
     * Adjusts real measurements to compensate for AFE gain
     * Arguments:
     *      FixedParameters : fixed parameters, which contain AFE gain adjustment
     * Returns:
     *      None
     ****************************************************************************/
    void Adjust(FixedParameters fparams);

private:
    std::vector<double> chargeValues;

};

