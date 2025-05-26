#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif


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

std::vector<double> ChargeOutput(std::vector<Measurement> input, VariableParameters vparams)
{
    std::vector<double> charge;
    std::vector<double> trigger;
    for(uint i = 0; i < input.size(); i++)
    {
        charge.push_back(input[i].charge);
        trigger.push_back(input[i].time);
    }

    // The charge needs to be integrated in a time bin

    double q = 0;
    std::vector<double> output;

    uint pointer = 0;
    // TODO: fix potential out of bounds error - For now, restricted for loop to not include the last 10 elements
    for(uint i = 0; i < trigger.size() - 10; i++)
    {
        if(pointer >= trigger.size() - 10)
            break;

        q = charge[pointer];
        int j = 1;
        while(trigger[pointer + j] - trigger[pointer] < vparams.gate)
        {
            q += charge[pointer + j];
            j += 1;
        }
        pointer += j;
        output.push_back(q);
    }

    return output;
}

std::vector<double> ParameterToVector(VariableParameters &vparams)
{
    std::vector<double> output;
    output.push_back(vparams.DCR);
    output.push_back(vparams.pde);
    output.push_back(vparams.pAPSHORT);
    output.push_back(vparams.tAPSHORT);
    output.push_back(vparams.pAPLONG);
    output.push_back(vparams.tAPLONG);
    output.push_back(vparams.pCT);
    output.push_back(vparams.tCT);
    output.push_back(vparams.jitter);
    output.push_back(vparams.pulseWidth);
    output.push_back(vparams.rechargeTime);
    output.push_back(vparams.gain);
    output.push_back(vparams.gainStd);
    output.push_back(vparams.gate);
    return output;
}

void UpdateParameters(std::vector<double> &vparams_vector, VariableParameters &vparams)
{
    vparams.DCR = vparams_vector[0];
    vparams.pde = vparams_vector[1];
    vparams.pAPSHORT = vparams_vector[2];
    vparams.tAPSHORT = vparams_vector[3];
    vparams.pAPLONG = vparams_vector[4];
    vparams.tAPLONG = vparams_vector[5];
    vparams.pCT = vparams_vector[6];
    vparams.tCT = vparams_vector[7];
    vparams.jitter = vparams_vector[8];
    vparams.pulseWidth = vparams_vector[9];
    vparams.rechargeTime = vparams_vector[10];
    vparams.gain = vparams_vector[11];
    vparams.gainStd = vparams_vector[12];
    vparams.gate = vparams_vector[13];

    return;
}

void PrintVector(std::vector<double> &vparams_vector, std::string legend)
{
    std::cout << legend;
    for(uint i = 0; i < vparams_vector.size(); i++)
    {
        std::cout << vparams_vector[i] << ", ";
    }
    std::cout << std::endl;
    return;
}

void RandomizeParameters(VariableParameters &vparams, double scale)
{
    vparams.DCR += gRandom->Gaus(0, vparams.DCR * scale);
    vparams.pde += gRandom->Gaus(0, vparams.pde * scale);
    vparams.pAPSHORT += gRandom->Gaus(0, vparams.pAPSHORT * scale);
    vparams.tAPSHORT += gRandom->Gaus(0, vparams.tAPSHORT * scale);
    vparams.pAPLONG += gRandom->Gaus(0, vparams.pAPLONG * scale);
    vparams.tAPLONG += gRandom->Gaus(0, vparams.tAPLONG * scale);
    vparams.pCT += gRandom->Gaus(0, vparams.pCT * scale);
    vparams.tCT += gRandom->Gaus(0, vparams.tCT * scale);
    vparams.jitter += gRandom->Gaus(0, vparams.jitter * scale);
    vparams.pulseWidth += gRandom->Gaus(0, vparams.pulseWidth * scale);
    vparams.rechargeTime += gRandom->Gaus(0, vparams.rechargeTime * scale);
    vparams.gain += gRandom->Gaus(0, vparams.gain * scale);
    vparams.gainStd += gRandom->Gaus(0, vparams.gainStd * scale);
    vparams.gate += gRandom->Gaus(0, vparams.gate * scale);

    return;
}


