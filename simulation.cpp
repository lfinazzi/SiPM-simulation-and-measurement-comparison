#ifndef SIMULATION_H
    #define SIMULATION_H
    #include "simulation.h"
#endif

#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif


void Simulation::Simulate()
{
    //auto start = std::chrono::high_resolution_clock::now();
    //std::cout << "\nStarting simulation...\n";

    double currentTime = 0;
    SimulationOutput simOut;

    // Generates photon events
    std::vector<double> triggers = GenerateTriggers(fparams, vparams);
    triggers.pop_back();

    for(uint l = 0; l < triggers.size(); l++)
    {
        currentTime = triggers[l];

        double detTime = currentTime + gRandom->Gaus(fparams.jitterLoc, vparams.jitter);

        double charge = gRandom->Gaus(vparams.gain, vparams.gainStd)*ELECTRONCHARGE;
        Measurement meas = CreateMeasurement(detTime, 1, charge, "P");
        
        //if(gRandom->Uniform(0, params.NMicrocells) > SaturationCount(params))
            simOut.measurements.push_back(meas);


        // Generates CTs
        double lamb = -TMath::Log(1 - vparams.pCT);
        std::vector<std::vector<int>> CTs = GenerateCTs(lamb);
        int numCrosstalks = CTs.size();
        
        for(uint i = 0; i < CTs.size(); i++){
            numCrosstalks += CTs[i][1];
        }

        for(uint i = 0; i < CTs.size(); i++){

            // Primary
            double CTtime = currentTime + gRandom->Exp(vparams.tCT) + gRandom->Gaus(fparams.jitterLoc, vparams.jitter);
            double charge = gRandom->Gaus(vparams.gain, vparams.gainStd)*ELECTRONCHARGE;
            Measurement meas = CreateMeasurement(CTtime, 1, charge, "CT1");

            //if(gRandom->Uniform(0, params.NMicrocells) > SaturationCount(params))
                simOut.measurements.push_back(meas);

            // Secondary
            for(int k = 0; k < CTs[i][1]; k++){
                double CTtime2 = currentTime + gRandom->Exp(vparams.tCT) + gRandom->Gaus(fparams.jitterLoc, vparams.jitter);
                double charge2 = gRandom->Gaus(vparams.gain, vparams.gainStd)*ELECTRONCHARGE;
                Measurement meas = CreateMeasurement(CTtime2, 1, charge2, "CT2");

                //if(gRandom->Uniform(0, params.NMicrocells) > SaturationCount(params))
                    simOut.measurements.push_back(meas);
            }
        }
        
        // Generates AP events
        for(int i = 0; i < 1 + numCrosstalks; i++)
        {
            // Short APs
            int nAPSHORT = gRandom->Poisson(-TMath::Log(1 - vparams.pAPSHORT));
            for(int j = 0; j < nAPSHORT; j++)
            {
                double apTime = gRandom->Exp(vparams.tAPSHORT);
                double charge = gRandom->Gaus(Gain(vparams, apTime), vparams.gainStd)*ELECTRONCHARGE;
                Measurement meas = CreateMeasurement(currentTime + apTime + gRandom->Gaus(fparams.jitterLoc, vparams.jitter), 1, charge, "APSHORT");
                simOut.measurements.push_back(meas);
            }

            // Long APs
            int nAPLONG = gRandom->Poisson(-TMath::Log(1 - vparams.pAPLONG));
            for(int j = 0; j < nAPLONG; j++)
            {
                double apTime = gRandom->Exp(vparams.tAPLONG);
                double charge = gRandom->Gaus(Gain(vparams, apTime), vparams.gainStd)*ELECTRONCHARGE;
                Measurement meas = CreateMeasurement(currentTime + apTime + gRandom->Gaus(fparams.jitterLoc, vparams.jitter), 1, charge, "APLONG");
                simOut.measurements.push_back(meas);
            }            
        }

    }

    SimulationOutput sorted = SortOutput(simOut);

    simData.push_back(sorted);
    numSimulations += 1;

    simData[simData.size()-1].charges = ChargeOutput(sorted.measurements, vparams);

    // Time after simulation execution
    //auto stop = std::chrono::high_resolution_clock::now();
    //std::cout << "Simulation completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " milliseconds.\n";

    return;
}

void Simulation::Simulate(int number)
{
    for(int i = 0; i < number; i++){   
        Simulate();
    }
    return;
}
