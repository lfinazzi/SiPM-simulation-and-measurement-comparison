#include "TH1D.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#define EVENT_NUMBER_IN_BYTES 16

class RealData{

public:
    // constructor --> initialized with path of .bin file
    RealData(std::string pathToFile);


    /****************************************************************************
     * Loads SiPM data from binary files in ./data/ folder
     * Arguments:
     *      std::string : path to .bin file
     * Returns:
     *      None
     ****************************************************************************/
    void LoadData(std::string path);


    /****************************************************************************
     * Returns chargeValues member
     * Arguments:
     *      None
     * Returns:
     *      std::vector<double> : charge values
     ****************************************************************************/
    inline std::vector<double> GetChargeData() { return chargeValues; }

private:
    std::vector<double> chargeValues;

};

