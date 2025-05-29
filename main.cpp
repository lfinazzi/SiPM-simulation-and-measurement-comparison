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

    histSimBefore->Scale(1.0/histSimBefore->GetEntries());      // normalizes the simHist to compare with normalized dataHist
    histSimAfter->Scale(1.0/histSimAfter->GetEntries());      // normalizes the simHist to compare with normalized dataHist

    // NLL values during minimization
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
    gS->SetTitle("NLL value during minimization;n;NLL");
    gS->SetLineColor(kBlue);
    gS->SetLineWidth(2);
    gS->Draw("AL");
    //gPad->SetLogy();

    // --- Pad 4: change in S vs iteration ---
    c->cd(4);
    gChangeS->SetTitle("Change in NLL value during minimization;n;#Delta NLL");
    gChangeS->SetLineColor(kBlue);
    gChangeS->SetLineWidth(2);
    gChangeS->Draw("AL");

    c->Update();
}

void PlotAllTiming(Minimizer minimizerBefore, Minimizer minimizerAfter)
{
    // Histograms before and after
    TH1D* histDataBeforeTiming = minimizerBefore.GetDataHistTiming();
    TH1D* histSimBeforeTiming = minimizerBefore.GetSimHistTiming();

    TH1D* histDataAfterTiming = minimizerAfter.GetDataHistTiming();
    TH1D* histSimAfterTiming = minimizerAfter.GetSimHistTiming();

    histSimBeforeTiming->Scale(1.0/histSimBeforeTiming->GetEntries());      // normalizes the simHist to compare with normalized dataHist
    histSimAfterTiming->Scale(1.0/histSimAfterTiming->GetEntries());      // normalizes the simHist to compare with normalized dataHist

    // NLL values during minimization
    std::vector<double> SvecTiming = minimizerAfter.GetSVectorTiming();
    int nPoints = SvecTiming.size();
    std::vector<double> xTiming(nPoints), yTiming(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        xTiming[i] = i;
        yTiming[i] = SvecTiming[i];
    }

    TGraph* gSTiming = new TGraph(nPoints, xTiming.data(), yTiming.data());

    // S change values during minimization
    std::vector<double> SChangevecTiming = minimizerAfter.GetSChangeVectorTiming();
    nPoints = SChangevecTiming.size();
    std::vector<double> xcTiming(nPoints), ycTiming(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        xcTiming[i] = i;
        ycTiming[i] = SChangevecTiming[i];
    }

    TGraph* gChangeSTiming = new TGraph(nPoints, xcTiming.data(), ycTiming.data());

    // --- Canvas and pads ---
    TCanvas* c2 = new TCanvas("cAll Timing", "Minimization Overview (Timing)", 1000, 800);
    c2->Divide(2, 2);

    histDataBeforeTiming->SetTitle("Timing Histogram before minimization;Time [s];Counts");
    // --- Pad 1: Histogram before ---
    c2->cd(1);
    histDataBeforeTiming->SetLineColor(kBlack);
    histSimBeforeTiming->SetLineColor(kRed);
    histDataBeforeTiming->Draw();
    histSimBeforeTiming->Draw("SAME");
    gPad->SetLogx();
    gPad->SetLogy();

    histDataAfterTiming->SetTitle("Timing Histogram after minimization;Time [s];Counts");
    // --- Pad 2: Histogram after ---
    c2->cd(2);
    histDataAfterTiming->SetLineColor(kBlack);
    histSimAfterTiming->SetLineColor(kRed);
    histDataAfterTiming->Draw();
    histSimAfterTiming->Draw("SAME");
    gPad->SetLogx();
    gPad->SetLogy();

    // --- Pad 3: S vs iteration ---
    c2->cd(3);
    gSTiming->SetTitle("NLL value during minimization (Timing);n;NLL");
    gSTiming->SetLineColor(kBlue);
    gSTiming->SetLineWidth(2);
    gSTiming->Draw("AL");
    //gPad->SetLogy();

    // --- Pad 4: change in S vs iteration ---
    c2->cd(4);
    gChangeSTiming->SetTitle("Change in NLL value during minimization (Timing);n;#Delta NLL");
    gChangeSTiming->SetLineColor(kBlue);
    gChangeSTiming->SetLineWidth(2);
    gChangeSTiming->Draw("AL");

    c2->Update();
}


int main(int argc, char *argv[])
{
	TApplication *myApp = new TApplication("myApp", &argc, argv, 0, -1);
    gRandom->SetSeed(0);
	gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    // Load measurements here
    int eventsToLoad = 3E6;
    RealData data("./data/measurement.bin", eventsToLoad);

    // initialize parameters with parameters.h file
    FixedParameters fparams;

    VariableParameters vparams;
    double randomScale = 0; // scale for randomization of parameters (to test Minimizer with slightly incorrect parameters)
    RandomizeParameters(vparams, randomScale);

    // Adjust data to compensate for AFE gain at acquisition (this allows simulation comparison)
    data.Adjust(fparams);

    // Initializes simulation
    int iterations = 50;
    double derivStepPerc = 1E-3;
    double startingLearningRate = 5E-3;
    
    Minimizer minimizer(fparams, vparams, data, iterations, derivStepPerc, startingLearningRate);

    // plot before minimization
    Minimizer min_before = minimizer;
    minimizer.RunMinimizer();

    // plots: charge histogram before and after minimization, and S as a function of iteration step
    PlotAll(min_before, minimizer);
    PlotAllTiming(min_before, minimizer);

    myApp->Run();
    
    return 0;
}
