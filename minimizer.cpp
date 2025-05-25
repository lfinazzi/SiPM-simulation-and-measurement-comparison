#ifndef MINIMIZER_H
    #define MINIMIZER_H
    #include "minimizer.h"
#endif


Minimizer::Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int _iterations, double _stepSize, double _delta)
    : sim(_fparameters, _vparameters), data(_data), fparams(_fparameters), vparams(_vparameters), iterations(_iterations), stepSize(_stepSize), delta(_delta)
{
    sim.Simulate();
    vparams_vector = ParameterToVector(vparams);

    dataHist.SetTitle("charge histogram data");
    simHist.SetTitle("charge histogram sim");

    dataHist.SetBins(300, 0, 1E-11);
    simHist.SetBins(300, 0, 1E-11);

    currentStep = 0;

    std::vector<double> chargesData = data.GetChargeData();
    std::vector<double> chargesSim = sim.SimDataCharge(currentStep);

    for(size_t i = 0; i < chargesData.size(); i++) {
        dataHist.Fill(chargesData[i]);
    }


    for(size_t i = 0; i < chargesSim.size(); i++) {
        simHist.Fill(chargesSim[i]);
    }

    simHist.SetLineColor(kRed);
    
    currentStep++;
    CalculateS(false);

    aux_vparams = vparams;
    aux_vparams_vector = ParameterToVector(aux_vparams);

    aux_simHist.SetTitle("aux charge histogram data");
    aux_simHist.SetTitle("aux charge histogram sim");

    aux_simHist.SetBins(300, 0, 1E-11);

}  

void Minimizer::CalculateS(bool aux)
{
    double s = 0;
    double simBinContent = 0;
    double dataBinContent = 0;
    for(int i = 1; i <= dataHist.GetNbinsX(); i++)
    {
        dataBinContent = dataHist.GetBinContent(i);
        if(aux == true)
            simBinContent = aux_simHist.GetBinContent(i);
        else
            simBinContent = simHist.GetBinContent(i);
        s += (dataBinContent - simBinContent)*(dataBinContent - simBinContent);
    }

    if(aux == true){
        aux_S = s;  // stores the current S value to aux_S
        std::cout << "aux_S (param + delta): " << aux_S << std::endl;
    }
    else{
        S = s;      // stores the current S value to class member  
        std::cout << "S (param): " << s << std::endl;
    }
       
}

void Minimizer::CalculateDerivative(int numParam)
{
    aux_vparams_vector[numParam] += delta;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    Simulation aux_sim(fparams, aux_vparams);
    aux_sim.Simulate();

    aux_simHist.Reset();
    std::vector<double> chargesSim = aux_sim.SimDataCharge(0);

    for(size_t i = 0; i < chargesSim.size(); i++) {
        aux_simHist.Fill(chargesSim[i]);
    }    

    CalculateS(false);
    CalculateS(true);

    double derivative = (aux_S - S) / (2 * delta);
    
    std::cout << "Derivative: " << derivative << std::endl;
}