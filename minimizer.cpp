#ifndef MINIMIZER_H
    #define MINIMIZER_H
    #include "minimizer.h"
#endif


Minimizer::Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int _iterations, double _stepPerc, double _mup)
    : sim(_fparameters, _vparameters), data(_data), fparams(_fparameters), vparams(_vparameters), iterations(_iterations), stepPerc(_stepPerc), minimUpdatePerc(_mup)
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
    CalculateS();

    aux_vparams = vparams;
    aux_vparams_vector = ParameterToVector(aux_vparams);

    aux_simHist.SetTitle("aux charge histogram data");
    aux_simHist.SetTitle("aux charge histogram sim");

    aux_simHist.SetBins(300, 0, 1E-11);

    gradient.resize(vparams_vector.size(), 0);

}  

void Minimizer::CalculateS()
{
    double s = 0;
    double simBinContent = 0;
    double dataBinContent = 0;
    for(int i = 1; i <= dataHist.GetNbinsX(); i++)
    {
        dataBinContent = dataHist.GetBinContent(i);
        simBinContent = simHist.GetBinContent(i);
        s += (dataBinContent - simBinContent)*(dataBinContent - simBinContent);
    }

    S = s;      // stores the current S value to class member  
       
}

void Minimizer::CalculateSExternal(TH1D hist, std::string option)
{
    double s = 0;
    double simBinContent = 0;
    double dataBinContent = 0;
    for(int i = 1; i <= hist.GetNbinsX(); i++)
    {
        dataBinContent = dataHist.GetBinContent(i);
        simBinContent = hist.GetBinContent(i);
        s += (dataBinContent - simBinContent)*(dataBinContent - simBinContent);
    }

    if(option == "plus")
        aux_S_plus = s;  // stores the current S value to aux_S_plus
    else if(option == "minus")
        aux_S_minus = s; // stores the current S value to aux_S_minus
    else
        std::cout << "Invalid option for CalculateSExternal" << std::endl;
}

void Minimizer::CalculateDerivative(int numParam)
{
    double deltaStep = vparams_vector[numParam] * stepPerc;
    aux_vparams_vector[numParam] += deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    Simulation aux_sim_plus(fparams, aux_vparams);
    aux_sim_plus.Simulate();

    TH1D aux_simHist("aux_simHist", "aux_simHist", 300, 0, 1E-11); 
    std::vector<double> chargesSim = aux_sim_plus.SimDataCharge(0);

    for(size_t i = 0; i < chargesSim.size(); i++) {
        aux_simHist.Fill(chargesSim[i]);
    }    

    CalculateSExternal(aux_simHist, "plus");

    aux_vparams_vector[numParam] -= 2 * deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    Simulation aux_sim_minus(fparams, aux_vparams);
    aux_sim_minus.Simulate();

    TH1D aux_simHist_m("aux_simHist_m", "aux_simHist_m", 300, 0, 1E-11); 
    std::vector<double> chargesSim_m = aux_sim_minus.SimDataCharge(0);

    for(size_t i = 0; i < chargesSim_m.size(); i++) {
        aux_simHist_m.Fill(chargesSim_m[i]);
    }    

    CalculateSExternal(aux_simHist_m, "minus");


    // reset aux_vparams to the original value
    aux_vparams_vector[numParam] += deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    return;
}

void Minimizer::CalculateGradient()
{
    for(size_t i = 0; i < vparams_vector.size(); i++)
    {
        CalculateDerivative(i);
        gradient[i] = (aux_S_plus - aux_S_minus) / (2 * vparams_vector[i] * stepPerc);
        //std::cout << "Gradient component: " << gradient[i] << std::endl;
    }
}

void Minimizer::Update()
{
    std::cout << "Updating parameters with ADAM..." << std::endl;
    // Initialize Adam state if needed
    if (m.empty()) {
        m.resize(vparams_vector.size(), 0.0);
        v.resize(vparams_vector.size(), 0.0);
    }

    // learning rate decays with each iteration
    double currentUpdatePerc = minimUpdatePerc / (1.0 + 1E-4 * t);
    t += 1;

    // --- Adam parameter update ---
    for (size_t i = 0; i < vparams_vector.size(); ++i)
    {
        // Compute biased first moment estimate
        m[i] = beta1 * m[i] + (1.0 - beta1) * gradient[i];

        // Compute biased second raw moment estimate
        v[i] = beta2 * v[i] + (1.0 - beta2) * gradient[i] * gradient[i];

        // Compute bias-corrected estimates
        double m_hat = m[i] / (1.0 - std::pow(beta1, t));
        double v_hat = v[i] / (1.0 - std::pow(beta2, t));


        // Parameter update
        vparams_vector[i] -= currentUpdatePerc * abs(vparams_vector[i]) * m_hat / (std::sqrt(v_hat) + epsilon);
    }
    
    UpdateParameters(vparams_vector, vparams);
    
    // update simulation
    sim.SetVParameters(vparams);
    sim.Simulate();
    
    // update histograms
    simHist.Reset();
    std::vector<double> chargesSim = sim.SimDataCharge(currentStep);
    
    for(size_t i = 0; i < chargesSim.size(); i++) {
        simHist.Fill(chargesSim[i]);
    }

    CalculateS();

    aux_vparams = vparams;
    aux_vparams_vector = ParameterToVector(aux_vparams);
    
    currentStep++;
}

void Minimizer::RunMinimizer()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    double prevS = 0;
    std::vector<double> prevParams;
    for(int i = 0; i < iterations; i++)
    {
        std::cout << "Iteration: " << i + 1 << std::endl;
        CalculateGradient();
        prevS = S;
        svector.push_back(S);   // stores S for later visualization
        prevParams = vparams_vector; // store previous parameters
        std::cout << "Previous S: " << prevS << std::endl;
        PrintVector(prevParams, "Previous parameters: ");
        std::cout << "----------------------------------------" << std::endl;
        PrintVector(gradient, "Gradient: ");
        std::cout << "----------------------------------------" << std::endl;
        Update();
        schangevector.push_back(S - prevS);   // stores S for later visualization
        std::cout << "current S: " << S << std::endl;
        PrintVector(vparams_vector, "Current parameters: ");
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "----------------------------------------\n" << std::endl;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << iterations << " iterations completed.\nMinimization completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()/1000 << " seconds.\n";

}