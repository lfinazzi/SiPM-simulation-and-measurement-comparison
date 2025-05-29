#ifndef REALDATA_H
    #define REALDATA_H
    #include "realdata.h"
#endif

#ifndef AUX_FUNCTIONS_H
    #define AUX_FUNCTIONS_H
    #include "aux_functions.h"
#endif

RealData::RealData(std::string pathToFile, int eventsToLoad)
{
    LoadData(pathToFile, eventsToLoad);
}

void RealData::LoadData(std::string pathToFile, int eventsToLoad)
{
	/* Opens binary files */
	FILE *f0 = fopen(pathToFile.c_str(), "rb");

    if (!f0) {
        std::cerr << "Error opening file/s!" << std::endl;
        return;
    }

	fseek(f0, 0, SEEK_END);		// seek to end of file
	int fileSize0 = ftell(f0);	// get current file pointer
	fseek(f0, 0, SEEK_SET);		// seek back to beginning of file

    int totalEvents = fileSize0 / EVENT_NUMBER_IN_BYTES;

    // prevents out of bounds error
    int events0 = (totalEvents > eventsToLoad) ? eventsToLoad : totalEvents;

	std::cout << "File opened (" << fileSize0 / 1E9 << " GB, " << events0 / 1E6 << "M events will be loaded)" << std::endl;

    int iter0 = 1;
    uint16_t int16_temp;
    uint16_t int16_temp2;
    uint32_t int32_temp;
    bool ret = false;

    std::vector<double> timestamps;

    while(iter0 < events0){
        // ETT
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // TT
        ret = fread(&int32_temp, sizeof(uint32_t), 1, f0);
        // fine TT
        ret = fread(&int16_temp2, sizeof(uint16_t), 1, f0);

        // TODO: check if this is the correct formula
        double timetag = (double)(int16_temp << 16 | int32_temp) + (double)int16_temp2 / 1024.0;
        timestamps.push_back( timetag * 1E-9 ); // convert to seconds

        // chage (short gate)
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // charge (long gate)
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);

        chargeValues.push_back((double)int16_temp);
        
        // baseline
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // baseline std    
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);

        iter0++;
    }

    timingValues = Diffs(timestamps);

    std::cout << "data loading successful (ret code: " << ret << ")\n";
	return;
}

void RealData::Adjust(FixedParameters fparams)
{
    for(size_t i = 0; i < chargeValues.size(); i++) {
        chargeValues[i] *= fparams.AFEadjustment;
    }
    return;
}