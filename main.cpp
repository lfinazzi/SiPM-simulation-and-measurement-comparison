#ifndef REALDATA_H
    #define REALDATA_H
    #include "realdata.h"
#endif

#ifndef SIMULATION_H
    #define SIMULATION_H
    #include "simulation.h"
#endif

#include "TCanvas.h"
#include "TApplication.h"
#include "TRandom.h"
#include "TStyle.h"

void PlotData(RealData data)
{
    TH1D* chargeHistogram = new TH1D("charge histogram", "charge histogram", 100, 0, 1E-9);
    std::vector<double> charges = data.GetChargeData();

    for(size_t i = 0; i < charges.size(); i++) {
        chargeHistogram->Fill(charges[i] * 40E-15);     // acquired with 40 fC/LSB
    }

    TCanvas *c = new TCanvas("c", "Charge Histogram", 800, 600);
    chargeHistogram->Draw();
    c->SetLogy();
    c->Update();

    return;
}

void PlotData(Simulation sim)
{
    TH1D* chargeHistogram = new TH1D("charge histogram", "charge histogram", 100, 0, 1E-9);

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

void PlotData(RealData data, Simulation sim)
{
    TH1D* chargeHistogramData = new TH1D("charge histogram data", "charge histogram sim", 100, 0, 1E-9);
    TH1D* chargeHistogramSim = new TH1D("charge histogram data", "charge histogram sim", 100, 0, 1E-9);
    
    std::vector<double> chargesData = data.GetChargeData();
    std::vector<double> chargesSim = sim.SimDataCharge(0);

    for(size_t i = 0; i < chargesData.size(); i++) {
        chargeHistogramData->Fill(chargesData[i] * 40E-15);     // acquired with 40 fC/LSB
    }


    for(size_t i = 0; i < chargesSim.size(); i++) {
        chargeHistogramSim->Fill(chargesSim[i]);
    }

    TCanvas *c = new TCanvas("c", "Charge Histogram", 800, 600);
    
    chargeHistogramSim->SetLineColor(kRed);
    
    chargeHistogramData->Draw("E");
    chargeHistogramSim->Draw("E SAME");
    c->SetLogy();
    c->Update();

    return;
}

int main(int argc, char *argv[])
{
	TApplication *myApp = new TApplication("myApp", &argc, argv, 0, -1);
    gRandom->SetSeed(0);
	gStyle->SetOptFit(1111);

    // Load real measurements here

    // initialize parameters with parameters.h file
    FixedParameters fparams;
    VariableParameters vparams;

    // Initializes simulation
    Simulation sim(fparams, vparams);

    // Runs simulation
    sim.Simulate();

    // Loads real data
    RealData data("./data/measurement.bin");

    // plots  finger spectrum for loaded data and simulation
    PlotData(data, sim);


    myApp->Run();
    
    return 0;
}
