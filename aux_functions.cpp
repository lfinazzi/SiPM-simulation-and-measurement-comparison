#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif

#include <iostream>
#include <fstream>

#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"


double Minimum(double value1, double value11, double value2, double value22)
{
    if(value1 < value2)
        return value11;
    else
        return value22;
}


std::vector<std::vector<int>> GenerateCTs(double mu)
{
    int n = gRandom->Poisson(mu);
    std::vector<std::vector<int>> output;
    std::vector<int> firstCT = {1};
    for(int i = 0; i < n; i++){
        output.push_back(firstCT);
    }

    for(int i = 0; i < n; i++){
        output[i].push_back(gRandom->Poisson(mu));
    }
    return output;
}


double Gain(VariableParameters params, double t)
{
    return params.gain*(1 - TMath::Exp(-t/params.rechargeTime));
}


std::vector<double> GenerateTriggers(FixedParameters fparams, VariableParameters vparams)
{
    std::vector<double> photons;
    std::vector<double> dcr;

    double t = 0;
    if(fparams.photonRate == 0)
        goto dcr;

    while(t < fparams.T){
            double time = gRandom->Exp(1/(fparams.photonRate*vparams.pde));
            photons.push_back(t + time);
            t += time;
    }   

    
    t = 0;
    dcr:
    while(t < fparams.T){
        double time = gRandom->Exp(1/vparams.DCR);
        dcr.push_back(t + time);
        t += time;
    }
    photons.insert(photons.end(), dcr.begin(), dcr.end());
    std::sort(photons.begin(), photons.end());
    return photons;

}


Measurement CreateMeasurement(double _time, int _number, double _charge, std::string _label)
{
    Measurement output;
    output.time = _time;
    output.number = _number;
    output.charge = _charge;
    output.label = _label;
    return output;
}


std::vector<double> Diffs(std::vector<Measurement> measurements)
{
    std::vector<double> timeDiffs;
    if(measurements.size() < 2)
        return {0};

    for(uint i = 0; i < measurements.size() - 1; i++)
    {
        timeDiffs.push_back(measurements[i + 1].time - measurements[i].time);
    }
    return timeDiffs;
}


std::vector<double> Diffs(std::vector<double> measurements)
{
    std::vector<double> timeDiffs;
    if(measurements.size() < 2)
        return {0};

    for(uint i = 0; i < measurements.size() - 1; i++)
    {
        timeDiffs.push_back(measurements[i + 1] - measurements[i]);
    }
    return timeDiffs;    
}

std::vector<double> LogSpace(double start, double stop, int N)
{
    double realStart = TMath::Power(10, start);
    double realBase = TMath::Power(10, (stop - start)/N);

    std::vector<double> retval;
    retval.reserve(N);
    std::generate_n(std::back_inserter(retval), N, Logspace<>(realStart, realBase));
    return retval;
}


void FillHist(TH1F* hist, std::vector<double> counts)
{
    for(uint i = 0; i < counts.size(); i++)
    {
        hist->Fill(counts[i]);
    }
    return;
}


void FillHist2D(TH2F* hist, std::vector<double> counts1, std::vector<double> counts2)
{
    for(uint i = 0; i < counts1.size(); i++)
    {
        hist->Fill(counts1[i], counts2[i]);
    }
    return;
}


SimulationOutput SortOutput(SimulationOutput simOut)
{
    std::vector<int> indexVec;
    for (std::size_t i = 0; i != simOut.measurements.size(); ++i) { 
        indexVec.push_back(i); 
    }

    std::sort(indexVec.begin(), indexVec.end(), [&](std::size_t i, std::size_t j) { return simOut.measurements[i].time < simOut.measurements[j].time; });

    SimulationOutput sortedOut;
    for (std::size_t i = 0; i != indexVec.size(); ++i)
    {
        sortedOut.measurements.push_back(simOut.measurements[indexVec[i]]);
    }
    
    return sortedOut;
}

/*
double SaturationCount(Parameters params)
{
    return params.NMicrocells*( 1 - TMath::Exp(-(params.photonRate*params.pde + params.DCR)*(1 + params.pCT) / params.NMicrocells) );
}
*/

double Var(std::vector<double> v)
{
    double sum = 0;
    uint N = v.size();

    double mean = 0;
    for(uint i = 0; i < v.size(); i++){
        mean += v[i] / N;
    }
    for(uint i = 0; i < v.size(); i++){
        sum += (v[i] - mean)*(v[i] - mean);
    }

    return sum / (N - 1);
}


double Mean(std::vector<double> v)
{
    uint N = v.size();

    double mean = 0;
    for(uint i = 0; i < v.size(); i++){
        mean += v[i] / N;
    }

    return mean;
}

/*
std::vector<double> RemoveSaturatedCounts(Parameters params, std::vector<double> v)
{
    int photons = params.photonRate*params.pde*params.rechargeTime;
    std::vector<double> output;

    for(uint i = 0; i < v.size(); i++)
    {
        if(gRandom->Uniform() > photons/double(params.NMicrocells)){
            output.push_back(v[i]);
        }
    }

    return output;
}
*/

void PlotHistogram(std::string name, int numBins, double min, double max)
{
    std::fstream ifs;
    double d;
    std::vector<double> data;
    ifs.open(name);
    if(ifs.is_open())
    {
        while(ifs >> d)
        {
            data.push_back(d);
        }
    }
    else
        std::cout << "Couldn't open specified file.\n";

    TCanvas* customCanvas = new TCanvas("Custom canvas", "Custom canvas", 800, 600);
    TH1F* customHist = new TH1F("customHist", "customHist", numBins, min, max);

    FillHist(customHist, data);
    customHist->Draw("E1");
    customHist->GetXaxis()->SetTitle("Time [s]");
    customHist->GetYaxis()->SetTitle("Entries");

    // Custom fit function for histogram - Currently to test hypothesis that timing arrival distribution is not uniform
    TF1* fitFunc = new TF1("fitFunc", "[0]", min, max);
    fitFunc->SetParameters(0);
    //customHist->Fit(fitFunc, "R");

    customCanvas->Update();
    
    return;
}

std::vector<double> ChargeOutput(std::vector<Measurement> input, VariableParameters vparams)
{
    std::vector<double> charge;
    std::vector<double> trigger;
    for(uint i = 0; i < input.size(); i++)
    {
        charge.push_back(input[i].charge);
        trigger.push_back(input[i].time);
    }

    // Now the charge needs to be integrated in a time bin

    double q = 0;
    std::vector<double> output;

    uint pointer = 0;
    for(uint i = 0; i < trigger.size() - 1; i++)
    {
        if(pointer >= trigger.size() - 1)
            break;

        q = charge[pointer];
        int j = 1;
        while(trigger[pointer + j] - trigger[pointer] < vparams.gate)
        {
            q += charge[pointer + j];
            j += 1;
        }
        pointer += j;
        output.push_back(q * vparams.AFEAmpGainFactor);
    }


    std::cout << "Done!\n";
    return output;
}