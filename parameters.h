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
    double gain                                             = 6E6;

    // Microcell standard deviation
    double gainStd                                          = 3.5E5;
        
    // charge integration gate
    double gate                                             = 75E-9;

    // Analog front end amplification gain
    // this is the amplification gain of measurements. AFE = TIA + Voltage amplifier
    double AFEAmpGainFactor                                 = 200;

};