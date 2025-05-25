#define SPEEDOFLIGHT 299792458
#define PLANKCONSTANT 6.626E-34
#define ELECTRONCHARGE 1.6E-19

struct FixedParameters{
    /***********************************************
    * General Parameters
    ************************************************/

    // Incoming photon rate [1/s]
    double photonRate                                       = 0 * 1E3;

    // Integration time/measurement window [s]
    double T                                                = 30;

    // Time bin or time discretization [s]
    double timeBin                                          = 500E-12;

    // Analog front end amplification gain [C]
    // this is the amplification gain of measurements. AFE = TIA + Voltage amplifier
    // Number between (...) was simulated with pyhton script -> TODO: Measure it
    // (...) is the integral of voltage pulse
    // 11 is the V/V gain of AFE
    // 10E3 is the TIA gain of AFE
    // 1.8E-10 is the 1 p.e. charge measurement
    // 40E-15 is the charge measurement sensitivity
    double AFEadjustment                                    = (4 * 22.5E-9) / 11 / 10E3 / 1.8E-10 * 40E-15; 
    
};

struct VariableParameters{
    /***********************************************
    * SiPM internal parameters
    ************************************************/

    // Dark count rate [1/s]
    double DCR                                              = 1.11E5;

    // Photon detection efficiency
    double pde                                              = 0.272;

    // Probability of short afterpulsing
    double pAPSHORT                                         = 0.02;

    // Time constant of short afterpulsing [s]
    double tAPSHORT                                         = 10E-9;

    // Probability of long afterpulsing
    double pAPLONG                                          = 0.07;

    // Time constant of long afterpulsing [s]
    double tAPLONG                                          = 100E-9;

    // Probability of crosstalk event
    double pCT                                              = 0.15;

    // Time constant of crosstalk event [s]
    double tCT                                              = 1E-9;

    // Timing jitter standard deviation [s]
    double jitter                                           = 565E-12 / 2.355;

    // Timing jitter mean [s]
    double jitterLoc                                        = 0;

    // Pulse width of microcell detection [s]
    double pulseWidth                                       = 11E-9;

    // Microcell recharge time [s]
    double rechargeTime                                     = 82E-9;

    // SiPM gain
    double gain                                             = 5E6;

    // Microcell standard deviation
    double gainStd                                          = 3.5E5;
        
    // charge integration gate
    double gate                                             = 75E-9;

};