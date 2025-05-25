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

#include "TCanvas.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TStyle.h"

void PlotData(RealData data)
{
    TH1D* chargeHistogram = new TH1D("charge histogram", "charge histogram", 100, 0, 1E-11);
    std::vector<double> charges = data.GetChargeData();

    for(size_t i = 0; i < charges.size(); i++) {
        chargeHistogram->Fill(charges[i]);
    }

    TCanvas *c = new TCanvas("c", "Charge Histogram", 800, 600);
    chargeHistogram->Draw();
    c->SetLogy();
    c->Update();

    return;
}

void PlotData(Simulation sim)
{
    TH1D* chargeHistogram = new TH1D("charge histogram", "charge histogram", 100, 0, 1E-11);

    std::vector<double> charges = sim.SimDataCharge(0);

    for(size_t i = 0; i < charges.size(); i++) {
        chargeHistogram->Fill(charges[i]);
    }

    TCanvas *c = new TCanvas("c", "Charge Histogram", 800, 600);
    chargeHistogram->Draw();
    c->SetLogy();
    c->Update();

    return;
}

void PlotData(Minimizer minimizer)
{
    TH1D* chargeHistogramData = minimizer.GetDataHist();
    TH1D* chargeHistogramSim = minimizer.GetSimHist();
    
    TCanvas *c = new TCanvas("c", "Charge Histogram", 800, 600);
    
    chargeHistogramData->Draw();
    chargeHistogramSim->Draw("SAME");
    c->SetLogy();
    c->Update();

    return;
}


int main(int argc, char *argv[])
{
	TApplication *myApp = new TApplication("myApp", &argc, argv, 0, -1);
    gRandom->SetSeed(0);
	gStyle->SetOptFit(1111);

    // Load measurements here
    RealData data("./data/measurement.bin");

    // initialize parameters with parameters.h file
    FixedParameters fparams;
    VariableParameters vparams;

    // Adjust data to compensate for AFE gain at acquisition (this allows simulation comparison)
    data.Adjust(fparams);

    // Initializes simulation
    int iterations = 10;
    double stepSize = 1;
    double delta = 0.001;
    Minimizer minimizer(fparams, vparams, data, iterations, stepSize, delta);

    // plots  finger spectrum for loaded data and simulation
    //PlotData(minimizer.GetSim());
    PlotData(minimizer);

    minimizer.CalculateDerivative(12);
    minimizer.CalculateDerivative(12);


    myApp->Run();
    
    return 0;
}
