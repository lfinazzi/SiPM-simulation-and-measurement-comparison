#ifndef MINIMIZER_H
    #define MINIMIZER_H
    #include "minimizer.h"
#endif


Minimizer::Minimizer(FixedParameters _fparameters, VariableParameters _vparameters, RealData _data, int _iterations, double _derivStepPerc, double _slr)
    : sim(_fparameters, _vparameters), data(_data), fparams(_fparameters), vparams(_vparameters), iterations(_iterations), derivStepPerc(_derivStepPerc), startingLearningRate(_slr)
{
    sim.Simulate();
    vparams_vector = ParameterToVector(vparams);

    dataHist.SetTitle("charge histogram data");
    simHist.SetTitle("charge histogram sim");

    dataHist.SetBins(histBins, 0, maxHistVal);
    simHist.SetBins(histBins, 0, maxHistVal);

    currentStep = 0;

    std::vector<double> chargesData = data.GetChargeData();
    std::vector<double> chargesSim = sim.SimDataCharge(currentStep);

    for(size_t i = 0; i < chargesData.size(); i++) {
        dataHist.Fill(chargesData[i]);
    }
    dataHist.Scale(1.0/dataHist.GetEntries());      // treats data histogram as the real probability distribution


    for(size_t i = 0; i < chargesSim.size(); i++) {
        simHist.Fill(chargesSim[i]);
    }

    simHist.SetLineColor(kRed);
    
    currentStep++;
    CalculateNLL();

    aux_vparams = vparams;
    aux_vparams_vector = ParameterToVector(aux_vparams);

    aux_simHist.SetTitle("aux charge histogram data");
    aux_simHist.SetTitle("aux charge histogram sim");

    aux_simHist.SetBins(histBins, 0, maxHistVal);

    gradient.resize(vparams_vector.size(), 0);

    // this leaves DCR fixed in minimization (TODO: Add variable to FixedParameters)
    //vparams_update[0] = false; // DCR, currently minimized

}  

void Minimizer::CalculateNLL()
{
    double nll = 0;
    double simBinContent = 0;
    double dataBinContent = 0;
    for(int i = 1; i <= dataHist.GetNbinsX(); i++)
    {
        dataBinContent = dataHist.GetBinContent(i);
        simBinContent = simHist.GetBinContent(i);
        //nll -= dataBinContent * log(simBinContent + 1E-8) - simBinContent * ( log(simBinContent + 1E-8) - 1 ); 
        nll -= dataBinContent * log(simBinContent + 1E-8);
    }

    NLL = nll;      // stores the current NLL value to class member  
       
}

void Minimizer::CalculateNLLExternal(TH1D hist, std::string option)
{
    double nll = 0;
    double simBinContent = 0;
    double dataBinContent = 0;
    for(int i = 1; i <= hist.GetNbinsX(); i++)
    {
        dataBinContent = dataHist.GetBinContent(i);
        simBinContent = hist.GetBinContent(i);        
        //nll -= dataBinContent * log(simBinContent + 1E-8) - simBinContent * ( log(simBinContent + 1E-8) - 1 ); 
        nll -= dataBinContent * log(simBinContent + 1E-8);
    }

    if(option == "plus")
        aux_NLL_plus = nll;  // stores the current NLL value to aux_NLL_plus
    else if(option == "minus")
        aux_NLL_minus = nll; // stores the current NLL value to aux_NLL_minus
    else
        std::cout << "Invalid option for CalculateNLLExternal" << std::endl;
}

void Minimizer::CalculateDerivative(int numParam)
{
    double deltaStep = vparams_vector[numParam] * derivStepPerc;
    aux_vparams_vector[numParam] += deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    Simulation aux_sim_plus(fparams, aux_vparams);
    aux_sim_plus.Simulate();

    TH1D aux_simHist("aux_simHist", "aux_simHist", histBins, 0, maxHistVal); 
    std::vector<double> chargesSim = aux_sim_plus.SimDataCharge(0);

    for(size_t i = 0; i < chargesSim.size(); i++) {
        aux_simHist.Fill(chargesSim[i]);
    }    

    CalculateNLLExternal(aux_simHist, "plus");

    aux_vparams_vector[numParam] -= 2 * deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    Simulation aux_sim_minus(fparams, aux_vparams);
    aux_sim_minus.Simulate();

    TH1D aux_simHist_m("aux_simHist_m", "aux_simHist_m", histBins, 0, maxHistVal); 
    std::vector<double> chargesSim_m = aux_sim_minus.SimDataCharge(0);

    for(size_t i = 0; i < chargesSim_m.size(); i++) {
        aux_simHist_m.Fill(chargesSim_m[i]);
    }    

    CalculateNLLExternal(aux_simHist_m, "minus");

    // reset aux_vparams to the original value
    aux_vparams_vector[numParam] += deltaStep;
    UpdateParameters(aux_vparams_vector, aux_vparams);

    return;
}

void Minimizer::CalculateGradient()
{
    for(size_t i = 0; i < vparams_vector.size(); i++)
    {
        std::vector<double> gradMeans;
        for(int j = 0; j < gradientIters; j++)
        {
            CalculateDerivative(i);
            gradMeans.push_back( (aux_NLL_plus - aux_NLL_minus) / (2 * vparams_vector[i] * derivStepPerc) );
        }
        gradient[i] = Mean(gradMeans); // calculate the mean of the gradient values
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
    double lr = startingLearningRate;// / (1.0 + 1E-4 * t);
    t += 1;

    NormalizeParameters(); // normalize parameters for minimization

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
        if(vparams_update[i] == false)
            lr = 0;
        else
            lr = startingLearningRate;  // (1.0 + 1E-4 * t); // learning rate decays with each iteration
        vparams_vector[i] -= lr * m_hat / (std::sqrt(v_hat) + epsilon);
    }

    DenormalizeParameters(); // denormalize parameters for simulation
    
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

    CalculateNLL();

    aux_vparams = vparams;
    aux_vparams_vector = ParameterToVector(aux_vparams);
    
    currentStep++;
}

void Minimizer::RunMinimizer()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    double prevNLL = 0;
    double gradNorm = 0;
    std::vector<double> prevParams;
    std::vector<double> nllVals;
    int i = 0;
    for(i = 0; i < iterations; i++)
    {
        std::cout << "Iteration: " << i + 1 << std::endl;
        CalculateGradient();
        prevNLL = NLL;
        svector.push_back(NLL);   // stores NLL for later visualization
        prevParams = vparams_vector; // store previous parameters
        std::cout << "Previous NLL: " << prevNLL << std::endl;
        PrintVector(prevParams, "Previous parameters: ");
        std::cout << "----------------------------------------" << std::endl;
        //PrintVector(gradient, "Gradient: ");
        
        gradNorm = 0;
        for (double g : gradient) gradNorm += g * g;
        gradNorm = std::sqrt(gradNorm);
        std::cout << "||Gradient|| : " << gradNorm << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        
        Update();
        schangevector.push_back(NLL - prevNLL);   // stores NLL for later visualization
        std::cout << "current NLL: " << NLL << std::endl;
        PrintVector(vparams_vector, "Current parameters: ");
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "----------------------------------------\n" << std::endl;

        bool cancelRun = true;
        if(i >= 5){
            for(int j = 0; j < 5; j++){
                cancelRun &= (svector[i - j] > svector[i - j - 1]);
                
            }
            if (cancelRun){
                std::cout << "NLL is increasing, stopping minimization." << std::endl;
                goto end;
            }
        }
    }

    end:
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << i << " iterations completed.\nMinimization completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()/1000 << " seconds.\n";

    // prints final parameters to console
    std::cout << std::endl;
    PrintParameters();

}

void Minimizer::NormalizeParameters()
{
    vparams_vector[0] /= 1E6;       // DCR    
    vparams_vector[3] /= 100E-9;    // tAPSHORT
    vparams_vector[5] /= 1E-6;      // tAPLONG
    vparams_vector[7] /= 10E-9;     // tCT
    vparams_vector[8] /= 1E-9;      // jitter
    vparams_vector[9] /= 100E-9;    // pulseWidth
    vparams_vector[10] /= 100E-9;   // rechargeTime
    vparams_vector[11] /= 1E7;      // gain
    vparams_vector[12] /= 1E6;      // gainStd
}

void Minimizer::DenormalizeParameters()
{
    vparams_vector[0] *= 1E6;       // DCR    
    vparams_vector[3] *= 100E-9;    // tAPSHORT
    vparams_vector[5] *= 1E-6;      // tAPLONG
    vparams_vector[7] *= 10E-9;     // tCT
    vparams_vector[8] *= 1E-9;      // jitter
    vparams_vector[9] *= 100E-9;    // pulseWidth
    vparams_vector[10] *= 100E-9;   // rechargeTime
    vparams_vector[11] *= 1E7;      // gain
    vparams_vector[12] *= 1E6;      // gainStd
}

void Minimizer::PrintParameters()
{
    std::cout << "Parameters: " << std::endl;
    std::cout << "DCR: " << vparams_vector[0] / 1E6 << " [MHz]" << std::endl;
    std::cout << "pde: " << vparams_vector[1] << std::endl;
    std::cout << "pAPSHORT: " << vparams_vector[2] << std::endl;
    std::cout << "tAPSHORT: " << vparams_vector[3] / 1E-9 << " [ns]" << std::endl;
    std::cout << "pAPLONG: " << vparams_vector[4] << std::endl;
    std::cout << "tAPLONG: " << vparams_vector[5] / 1E-9 << " [ns]" << std::endl;
    std::cout << "pCT: " << vparams_vector[6] << std::endl;
    std::cout << "tCT: " << vparams_vector[7] / 1E-9 << " [ns]" << std::endl;
    std::cout << "jitter: " << vparams_vector[8] / 1E-9 << " [ns]" << std::endl;
    std::cout << "pulseWidth: " << vparams_vector[9] / 1E-9 << " [ns]" << std::endl;
    std::cout << "rechargeTime: " << vparams_vector[10] / 1E-9 << " [ns]" << std::endl;
    std::cout << "gain: " << vparams_vector[11] << std::endl;
    std::cout << "gainStd: " << vparams_vector[12] << std::endl;
}