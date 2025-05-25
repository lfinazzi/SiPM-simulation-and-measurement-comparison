#ifndef REALDATA_H
    #define REALDATA_H
    #include "realdata.h"
#endif

RealData::RealData(std::string pathToFile)
{
    LoadData(pathToFile);
}

void RealData::LoadData(std::string pathToFile)
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

    int events0 = fileSize0 / EVENT_NUMBER_IN_BYTES;

	std::cout << "File opened (" << fileSize0 / 1E9 << " GB, " << events0 / 1E6 << "M events)" << std::endl;

    int iter0 = 1;
    uint16_t int16_temp;
    uint32_t int32_temp;
    bool ret = false;

    while(iter0 < events0){
        // ETT
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // TT
        ret = fread(&int32_temp, sizeof(uint32_t), 1, f0);
        // fine TT
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // chage (short gate)
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        chargeValues.push_back((double)int16_temp);
        // charge (long gate)
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // baseline
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);
        // baseline std    
        ret = fread(&int16_temp, sizeof(uint16_t), 1, f0);

        iter0++;
    }

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