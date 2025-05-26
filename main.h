#ifndef REALDATA_H
    #define REALDATA_H
    #include "realdata.h"
#endif

#ifndef SIMULATION_H
    #define SIMULATION_H
    #include "simulation.h"
#endif

#ifndef MINIMIZER_H
    #define MINIMIZER_H
    #include "minimizer.h"
#endif

#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif

#include "TCanvas.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TGraph.h"


/****************************************************************************
 * Plots data from measurements and simulation before minimization
 * Arguments:
 *      Minimizer : object containing RealData and Simulation objects
 * Returns:
 *      None
 ****************************************************************************/
void PlotData(Minimizer minimizer);


/****************************************************************************
 * Plots data from measurements and simulation after minimization
 * Arguments:
 *      Minimizer : object containing RealData and Simulation objects
 * Returns:
 *      None
 ****************************************************************************/
void PlotData2(Minimizer minimizer);


/****************************************************************************
 * Plots data from measurements and simulation before minimization, and
 * S value during minimization in a 2x2 canvas.
 * Arguments:
 *      Minimizer : object containing RealData and Simulation objects before
 *               minimization
 *      Minimizer : object containing RealData and Simulation objects after
 *               minimization
 * Returns:
 *      None
 ****************************************************************************/
void PlotAll(Minimizer minimizerBefore, Minimizer minimizerAfter);