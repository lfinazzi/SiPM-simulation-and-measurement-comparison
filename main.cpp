#ifndef MAIN_H
    #define MAIN_H
    #include "main.h"
#endif

void PlotData(Minimizer minimizer)
{
    TH1D* chargeHistogramData = minimizer.GetDataHist();
    TH1D* chargeHistogramSim = minimizer.GetSimHist();
    
    TCanvas *c = new TCanvas("c", "Charge Histogram before minimization", 800, 600);
    
    chargeHistogramData->Draw();
    chargeHistogramSim->Draw("SAME");
    c->SetLogy();
    c->Update();

    return;
}

void PlotData2(Minimizer minimizer)
{
    TH1D* chargeHistogramData = minimizer.GetDataHist();
    TH1D* chargeHistogramSim = minimizer.GetSimHist();
    
    TCanvas *c2 = new TCanvas("c2", "Charge Histogram after minimization", 800, 600);
    
    chargeHistogramData->Draw();
    chargeHistogramSim->Draw("SAME");
    c2->SetLogy();
    c2->Update();

    return;
}

void PlotAll(Minimizer minimizerBefore, Minimizer minimizerAfter)
{
    // Histograms before and after
    TH1D* histDataBefore = minimizerBefore.GetDataHist();
    TH1D* histSimBefore = minimizerBefore.GetSimHist();

    TH1D* histDataAfter = minimizerAfter.GetDataHist();
    TH1D* histSimAfter = minimizerAfter.GetSimHist();

    // S values during minimization
    std::vector<double> Svec = minimizerAfter.GetSVector();
    int nPoints = Svec.size();
    std::vector<double> x(nPoints), y(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        x[i] = i;
        y[i] = Svec[i];
    }

    TGraph* gS = new TGraph(nPoints, x.data(), y.data());

    // S change values during minimization
    std::vector<double> SChangevec = minimizerAfter.GetSChangeVector();
    nPoints = SChangevec.size();
    std::vector<double> xc(nPoints), yc(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        xc[i] = i;
        yc[i] = SChangevec[i];
    }

    TGraph* gChangeS = new TGraph(nPoints, xc.data(), yc.data());

    // --- Canvas and pads ---
    TCanvas* c = new TCanvas("cAll", "Minimization Overview", 1000, 800);
    c->Divide(2, 2);

    histDataBefore->SetTitle("Charge Histogram before minimization;Charge [C];Counts");
    // --- Pad 1: Histogram before ---
    c->cd(1);
    histDataBefore->SetLineColor(kBlack);
    histSimBefore->SetLineColor(kRed);
    histDataBefore->Draw();
    histSimBefore->Draw("SAME");
    gPad->SetLogy();

    histDataAfter->SetTitle("Charge Histogram after minimization;Charge [C];Counts");
    // --- Pad 2: Histogram after ---
    c->cd(2);
    histDataAfter->SetLineColor(kBlack);
    histSimAfter->SetLineColor(kRed);
    histDataAfter->Draw();
    histSimAfter->Draw("SAME");
    gPad->SetLogy();

    // --- Pad 3: S vs iteration ---
    c->cd(3);
    gS->SetTitle("S value during minimization;n;S_{n}");
    gS->SetLineColor(kBlue);
    gS->SetLineWidth(2);
    gS->Draw("AL");
    gPad->SetLogy();

    // --- Pad 4: change in S vs iteration ---
    c->cd(4);
    gChangeS->SetTitle("Change in S value during minimization;n;S_{n} - S_{n-1}");
    gChangeS->SetLineColor(kBlue);
    gChangeS->SetLineWidth(2);
    gChangeS->Draw("AL");

    c->Update();
}


int main(int argc, char *argv[])
{
	TApplication *myApp = new TApplication("myApp", &argc, argv, 0, -1);
    gRandom->SetSeed(0);
	gStyle->SetOptFit(1111);

    // Load measurements here
    int eventsToLoad = 20E3;
    RealData data("./data/measurement.bin", eventsToLoad);

    // initialize parameters with parameters.h file
    FixedParameters fparams;

    VariableParameters vparams;
    double randomScale = 0.15; // scale for randomization of parameters (to test Minimizer with slightly incorrect parameters)
    RandomizeParameters(vparams, randomScale);

    // Adjust data to compensate for AFE gain at acquisition (this allows simulation comparison)
    data.Adjust(fparams);

    // Initializes simulation
    int iterations = 350;
    double stepPerc = 5E-3;
    double minimUpdatePerc = 0.05;
    
    Minimizer minimizer(fparams, vparams, data, iterations, stepPerc, minimUpdatePerc);

    // plot before minimization
    Minimizer min_before = minimizer;
    minimizer.RunMinimizer();

    // plots: charge histogram before and after minimization, and S as a function of iteration step
    PlotAll(min_before, minimizer);

    myApp->Run();
    
    return 0;
}
