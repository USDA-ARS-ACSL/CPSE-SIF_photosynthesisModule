/*! @file
*  Original version: No sif and no sunlit/shaded version
* -----------------------------------------------------------------  
* SIF-ENABLED-PHOTOSYNTHESIS
* MLR-SIF model integration into C3 and C4 photosynthesis models: VERSION 0.1.0-Dev
* -----------------------------------------------------------------
* Revised (current): Modified by Sahila Beegum, Christine Chang
* Basic Refernces:  https://www.mdpi.com/2223-7747/9/10/1358
* https://academic.oup.com/jxb/article/72/17/6003/6309859?login=true#301357310
* https://lieth.ucdavis.edu/pub/Pub058_AnnBot91-771_KimLieth.pdf
* -----------------------------------------------------------------

*  
*  Can simulate sif based C3, C4 photosynthesis
*  Have sunlit and shaded components
*  Can run multiple files at once: climate * sif* crop 
* 
* 
*  Defines the entry point for the console application.
* 
*  This program reads simulation run definitions from "run.dat" (each line formatted as:
*  ClimateIn.dat, maize, SIF1.dat).
*  ClimateIn.dat, Maize, SIF3.dat
*  ClimateIn.dat, Potato, SIF2.dat
*  ClimateIn.dat, Rice, SIF.dat
* 
*  It reads the entire "parameters.csv" file once into memory,
*  then for each run it looks up the species-specific parameters from that stored data.
* 
*  Next, it reads the climate and SIF data, runs the photosynthesis model for sunlit and shaded, and writes result
*  to "Results.dat".  
*/


#include "stdafx.h"       // Precompiled header (if used)
#include "pch.h"          // Precompiled header (if used)
#include "gas_exchange.h" // for CGasExchange and its tParms structure
#include "radtrans.h"     // for CRadTrans (radiation transfer)
#include "solar.h"        // for CSolar (solar related)
#include <iostream>       // For standard I/O (cout, cerr)
#include <fstream>        // For file I/O
#include <sstream>        // For stringstream
#include <vector>         // For std::vector
#include <algorithm>      // For std::remove and std::transform
#include <cctype>         // For std::isdigit
#include <string>         // For std::string
#include <chrono>
///#include "../FvCB_NotManaged/gas_exchange.h"


using namespace std;
using namespace photomod; // Bring the photosynthesis model namespace into scope

//units
//temperature (C)
//PAR (umol photons m-2 s-1)
//CO2 content (umol mol-1)
//humidity (%)
//wind (m s-1)
//a flag (0,1) to tell the program if constant temperature is used (or let temperature of leaf
 
//PAR                   (umol photons m-2 s-1)
//Net Photosynthesis    (umol CO2 m-2 s-1)
//Gross Photosynthesis
//VPD                   (kPa)
//LeafTemperature       (C)
//BoundaryLayerConductance (mol m-2 s-1)
//Internal CO2          (umol mol-1)
//Respiration           (umolCO2 m-2 s-1)
//Transpiration         (umol H2O m-2 s-1) \n


//---------------------------------------------------------------------
// Global constant: delimiter used for token CSV lines
const char* pDelim = ","; // pDelim is declared here

//---------------------------------------------------------------------
// Helper function -convert a string to lower case
// This function uses std::transform and ::tolower to convert all characters in the string
static inline string toLower(const string& s)
{
    string result = s;
    transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result; // return lowecase string
}


//---------------------------------------------------------------------
// Helper function: trim string spaces
// removes all spaces from anywhere in the string
// we can modify it if you only want to remove leading/trailing spaces
static inline string trim(const string& s)
{
    string result = s;
    result.erase(remove(result.begin(), result.end(), ' '), result.end());
    return result; 
}

//---------------------------------------------------------------------
// Function: parseParameterLine
// Parse one line from "parameters.csv" and store its tokens in a tParms structure.
// 24 Items:
//   1) ID, 2) Species, 3) type, 4) Vcm25, 5) Jm25, 6) Vpm25, 7) TPU25, 8) Rd25,
//   9) theta, 10) EaVc, 11) Eaj, 12) Hj, 13) Sj, 14) HV, 15) EaVp, 16) Sv,
//   17) EAP, 18) EaR, 19) g0, 20) g1, 21) stomRatio, 22) lfwidth, 23) lfAngleFactor, 24) remark.
bool parseParameterLine(const string& line, CGasExchange::tParms& param)
{
    //create a vector to hold individual variable from Parameters.csv
    vector<string> tokens;
    tokens.reserve(24); // Reserve space for 24 tokens


    // Use a stringstream to split the line by commas
    stringstream ss(line);
    string token;
    while (getline(ss, token, ','))  
    {
        // Trim spaces from each token and add to the vector
        tokens.push_back(trim(token));
    }


    // Verify that the line has 24 col
    if (tokens.size() != 24)
    {
        cerr << "Skipping line (does not have 24 columns): " << line << endl;
        return false;
    }


    // Assign tokens to the corresponding fields in tParms
    // The first field, ID, is converted to lowercase for case-insensitive matching.
    param.ID = toLower(tokens[0]); // ( "maize")
    param.species = tokens[1];          //!<Species name ( Zea mays L.)
    param.Type = tokens[2];          // Plant type ("C3" or "C4")
    param.Vcm25 = atof(tokens[3].c_str());
    param.Jm25 = atof(tokens[4].c_str());
    param.Vpm25 = atof(tokens[5].c_str());
    param.TPU25 = atof(tokens[6].c_str());
    param.Rd25 = atof(tokens[7].c_str());
    param.Theta = atof(tokens[8].c_str());
    param.EaVc = atof(tokens[9].c_str());
    param.Eaj = atof(tokens[10].c_str());
    param.Hj = atof(tokens[11].c_str());
    param.Sj = atof(tokens[12].c_str());
    param.Hv = atof(tokens[13].c_str());
    param.EaVp = atof(tokens[14].c_str());
    param.Sv = atof(tokens[15].c_str());
    param.Eap = atof(tokens[16].c_str());
    param.Ear = atof(tokens[17].c_str());
    param.g0 = atof(tokens[18].c_str());
    param.g1 = atof(tokens[19].c_str());
    param.stomaRatio = atof(tokens[20].c_str());
    param.LfWidth = atof(tokens[21].c_str());
    param.LfAngFact = atof(tokens[22].c_str());
    param.Remark = tokens[23]; // Remark field



    //                parameter        description
    //Vcm25	        4		Photosynthetic Rubisco Capacity at 25C(umol m - 2 s - 1)
    //Jm25	        5		Potential Rate of electron transport at 25C(umol m - 2 s - 1)
    //Vpm25	        6		C4 Carboxylation rate at 25C(C4, umol m - 2 s - 1)
    //TPU25	        7		Rate if Triose Phosphate Utilization at 25C(C3, umol m - 2 s - 1)
    //RD25	        8		Mitochondrial respiration in the light at 25C(umol m - 2 s - 1)
    //Theta         9		Initial slope of CO2 response(umol m2 s - 1) - de Pury(1997)
    //EaVc          10		Activation energy for Arrhenius function used to calculate temperature dependence for Vcmax(kJ mol - 1)
    //Eaj           11		Activation energy for Arrhenius function used to calculate temperature dependence for J(kJ mol - 1)
    //Hj            12		Curvature parameter of the temperature dpendence of Jmax(kJ mol - 1)
    //Sj	        13		Electron transport temperature response parameter for Jmax(J mole - 1 K - 1)
    //Hv            14       Curvature parameter of the temperature dependence of Vcmax(J mole - 1)
    //EaVp          15		Activation energy for Arrhenius function used to calculate temperature dependence for Vpmax(kJ mol - 1)
    //Sv	        16		Electron transport temperature response parameter for Vcmax(J mole - 1 K - 1)
    //EAP	        17		Activation energy for Arrhenius function used to calculate temperature dependence for TPU(kJ mol - 1)
    //EAR	        18	 	Activation energy for Arrhenius function used to calculate temperature dependence for respiration(kJ mol - 1)
    //g0	        19		Minimum stomatal conductance to water vapor at the light compensation point in the BWB model(mol m - 2 s - 1)
    //g1	        20   	Empirical coefficient for the sensitivity of StomatalConductance to A, Cs and hs in BWB model(no units)
    //StomRatio	    21	    Stomatal Ratio(fraction)
    //LfWidth       22	    Leaf Width(m)
    //LfAngFact	    23	    Leaf Angle Factor
    //Remark		24	    Text

    return true;
}


//---------------------------------------------------------------------
// Main function
int _tmain(int argc, _TCHAR* argv[])
{
   
    
    //-----------------------------------------------------------------
    // Step 1: Declare variables for environmental and other simulation data.
    //-----------------------------------------------------------------
    //PhotoFluxDensity	Photosynthetic Flux Density(umol Quanta m - 2 s - 1)
    //Tair	Air Temperature(C)
    //RH	Relative Humidity(%)
    //wind	Windspeed at 2.5 m, m s - 1
    //Press	Atmospheric pressure(kpa m - 2)
    cerr << "*******************************************************************\n"
        << "  SIF-ENABLED-PHOTOSYNTHESIS, VERSION 0.1.0 - Dev\n"
        << "  MLR - SIF model integration into C3 and C4 photosynthesis models\n"
        << "******************************************************************** \n\n";
    string DataLine; // DataLine is declared here for reading lines from files
    double PFD; //!<jj
    double Temperature, RelativeHumidity, Wind, CO2, LAI, ql_sunlit_shaded = 0.0;
    double kdf_sunlit_shaded = 0.0;
    double PhiPS2_sunlit_shaded=0.0;
    double Pressure = 100; // Pressure  
    double DayOfYear, Time, Lat, Lon, Alt, SolRad;
    double SifUnitConv=1, SIFtotalFactor=0.0;
    const double pi = 3.141592653589793;
 
    // SIF variables declared:
        //To simulate the photosythesis based on the SIF : steps for A using SIF is given in
        //https ://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18045

        //Paper for the SIF J equation : https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.15796
        //GasExchange simulate the shaded and sunlit part of the plant
        //There are some specific inputs for shaded and sunlit
        //Reading the inputs for those here : Jday, Date, Hour, ql_sunlit, ql_shaded,
        //Kdf_sunlit, Kdf_shaded, PhiPS2max sunlit, PhiPS2max shaded, e(will modify these)
 

        //light response curves for many crops : https://leafweb.org/information/data-publications/
        //qL can be estimated based on par.the equation can be a quadratic function(PAR)
        //qlshaded = a_shaded.x2 + b_shaded.x + c_shaded, x = par
        //qlsunlit = a_sunlit.x2 + b_sublit.x + c_sunlit, x = par
        //Inputs needed : are   a_shaded / sunlit, b_shaded / sunlit, c_shaded / sunlit

        //PhiPS2maxShaded and PhiPS2maxSunlit  can be constant for now - read as input
        //PhiPS2max default is 0.83 https ://link.springer.com/article/10.1007/bf00402983
        //https ://www.sciencedirect.com/science/article/pii/S0168192324003058#sec0026

        //kdf = KD / KF, Kd is constant, Kf is function of PAR. (beta * APAR / Fm - 1)
        //beta usually approximated as ~0.5 = the split between PSII and PSI
        //https ://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/nph.15796    
        //sifvalue total full band SIF emitted by PSII[μmol photons m - 2s - 1]

    double sifSunlit = 0.0, sifShaded = 0.0, ql_sunlit = 0.0, ql_shaded = 0.0, kdf_sunlit = 0.0, kdf_shaded = 0.0;
    double PhiPS2_sunlit = 0.0, PhiPS2_shaded = 0.0, fesc = 0.0;
    int SIF_GasEx = 0; // SIF_GasEx flag 
    double sunlitLAI = 0.0, shadedLAI = 0.0;
    double sunlitPFD = 0.0, shadedPFD = 0.0;
    double fesc_sun = 0.0, fesc_sha = 0.0;
    // double Anet = 0.0, Agross = 0.0, VPD = 0.0, LeafTemperature = 0.0;
    //double BoundaryLayerConductance = 0.0, InternalCO2 = 0.0, Respiration = 0.0;
    //double StomatalConductance = 0.0, Transpiration = 0.0;
   // double AjLeaf = 0.0, AcLeaf = 0.0, AvLeaf = 0.0, Jval =0.0;
   // double GammaVal = 0.0, GammaC4Val = 0.0, Cs_stom_conductanceVal = 0.0;
    bool ConstantTemperature; // ConstantTemperature flag
    FILE* pFile = nullptr;
    //ConstantTemperature boolian if true, leaf temperature=air temperature when calculating gas exchange (check)


    //-----------------------------------------------------------------
    // Step 2: Declare default file names and default species. (not sure if we need this)
    //-----------------------------------------------------------------
    string DataFileName = "ClimateIn.dat";
    string Species = "Maize"; // Default species  
    string sunlitShadedFlag = "SunlitShadedNo";
    string sifFileName = "";
    string outputFile = "Results.dat";  // default


    //-----------------------------------------------------------------
    // Step 3: Read the entire parameters.csv file in vector
    //-----------------------------------------------------------------
    vector<CGasExchange::tParms> allParams; // Vector to hold all param struct
    ifstream paramFile("parameters.csv"); //open file
    if (!paramFile)
    {
        cerr << "Parameter file not found.\n";
        return 1;
    }
    // Read and skip the header line
    getline(paramFile, DataLine);
    // Read each subsequent lines
    while (getline(paramFile, DataLine))
    {
        if (DataLine.empty()) // if empty skip it
            continue;
        // not sure why there was a line with 1,2,3 in the Parameters.csv
        // Skip lines that look like numeric index lines ("1,2,3,...,24").
        if (isdigit(DataLine[0]) && DataLine.find(',') != string::npos)
        {
            cerr << "Skipping numeric index line: " << DataLine << endl;
            continue;
        }
        CGasExchange::tParms param;// temporary tparms strucutre

        if (parseParameterLine(DataLine, param)) // if success,
        {
            allParams.push_back(param);
        }
    }
    paramFile.close();


    //-----------------------------------------------------------------
    // Step 4: Open the run file ("run.dat") and process each run.
    //-----------------------------------------------------------------
    ifstream runFile("run.dat");
    if (!runFile)
    {
        cerr << "Run file 'run.dat' not found. Exiting.\n";
        return 1;
    }


    //-----------------------------------------------------------------
    // Step 5: Process each run defined in run.dat.
    // Each run line expected: ClimateIn.dat, species, SIF.dat
    //-----------------------------------------------------------------
    string runLine;
    while (getline(runFile, runLine))
    {

        auto start = std::chrono::high_resolution_clock::now();
        if (runLine.empty())
            continue;


        // Use strtok_s to tokenize the run line
        char* runCtx = nullptr;
        string climateFile, runSpecies, sunlitShadedFlag = "SunlitShadedNo";
        string sifFileName = "";


        // Token 1: Climate file
        char* pnt = strtok_s((char*)runLine.c_str(), pDelim, &runCtx);
        if (pnt != nullptr)
            climateFile = trim(pnt);


        // Token 2: Species
        pnt = strtok_s(NULL, pDelim, &runCtx);
        if (pnt != nullptr)
            runSpecies = toLower(trim(pnt));
        else
            runSpecies = "maize";


        // Token 3: SunlitShadedYes or SunlitShadedNo
        pnt = strtok_s(NULL, pDelim, &runCtx);
        if (pnt != nullptr) {
            string token = toLower(trim(pnt));
            if (token == "sunlitshadedyes" || token == "sunlitshadedno")
                sunlitShadedFlag = (token == "sunlitshadedyes") ? "SunlitShadedYes" : "SunlitShadedNo";
        }


        // Token 4: SIF file name (can be empty)
        pnt = strtok_s(NULL, pDelim, &runCtx);
        if (pnt != nullptr)
            sifFileName = trim(pnt);


        // Token 5: Output file name
        pnt = strtok_s(NULL, pDelim, &runCtx);
        if (pnt != nullptr)
            outputFile = trim(pnt);
        else
            outputFile = "Results.dat"; // fallback


        string sFlag = sifFileName.empty() ? "SIFno" : "SIFyes";


        // Debug print
        cout << "Processing run: Species: " << runSpecies
            << ", Climate file: " << climateFile
            << ", SunlitShaded: " << sunlitShadedFlag
            << ", SIF file: " << (sifFileName.empty() ? "None" : sifFileName)
            << ", Output file: " << outputFile << endl;

        //-----------------------------------------------------------------
        // Step 6: Open the output file for writing summary report.
        //-----------------------------------------------------------------
        FILE* pFile = nullptr;
        if (fopen_s(&pFile, outputFile.c_str(), "w") != 0 || !pFile)
        {
            cerr << "Could not open output file.\n";
            return 1;
        }
        // Write the header line to the output file.
        fprintf(pFile,
            "Species, DayOfYear, Time, Lat, Long, Alt, SolRad, PFD, Temp, CO2, RH, Wind, LAI, ConstTemp, "
            "sifSunlit, sifShaded,ql_sunlit, ql_shaded, kdf_sunlit, kdf_shaded, PhiPS2_sunlit, PhiPS2_shaded, fesc, "
            "ANet, AGross, VPD, LeafTemp, BoundaryL_Conductance, InternalCO2, Respiration, StomatalConduct, "
            "Transpiration,J,sunlitPFD,shadedPFD,sunlitLAI,shadedLAI,"
            "GammaValue(gamma_C3 or GammaStar_C4), Gamma_C4,Cs_stomatalConductance, ci, AjLeaf,AcLeaf,AvLeaf,ApLeaf,"
            "Azimuth, DayLength, Declination, FracPARDiffuse, FracPARDirect, NIR,"
            "NIRFraction, NIRTotal, PAR ,SolarNoon, Sunrise, Sunset, SinElevation, SolarElevation,"
            "PotSolarDiffuse, PotSolarTotal, PotSolarDirect, MeasuredSolarRad, PARFraction,"
            "PFD, PFDDiffuse, PFDDirect, PFDTotal, Zenith, PotParDirect,PotParDiffuse, PotNIRDirect,PotNIRDiffuse\n"
        );



 


        //-----------------------------------------------------------------
        // Step 7: Look up the parameters for the species from the vector allParams.
        //-----------------------------------------------------------------
        CGasExchange::tParms chosenParams;
        bool foundParams = false;
        for (size_t i = 0; i < allParams.size(); i++)
        {
            if (toLower(allParams[i].ID) == runSpecies)
            {
                chosenParams = allParams[i]; // Copy parameters.
                foundParams = true;
                break;
            }
        }
        if (!foundParams)
        {
            cerr << "No parameters found for species '" << runSpecies << "' - skipping this run.\n";
            continue;
        }


        //-----------------------------------------------------------------
        // Step 8: Open the climate input datafile
        //-----------------------------------------------------------------
        ifstream dataFile(climateFile.c_str());
        if (!dataFile)
        {
            cerr << "Could not open climate file " << climateFile << " - skipping this run.\n";
            continue;
        }


        //-----------------------------------------------------------------
        // Step 9: Open the SIF file if SIF data are used.
        //-----------------------------------------------------------------
        ifstream sifHandle;
        if (sFlag == "SIFyes")
        {
            sifHandle.open(sifFileName.c_str());
            if (!sifHandle)
            {
                cerr << "Could not open SIF file " << sifFileName << " - skipping this run.\n";
                continue;
            }
        }


        //-----------------------------------------------------------------
        // Step 10: Create gas exchange objects for sunlit and shaded parts
        //-----------------------------------------------------------------
        CGasExchange* sunlit = new CGasExchange();
        CGasExchange* shaded = new CGasExchange();

        CGasExchange* sunlit_shaded = new CGasExchange();

        sunlit->SetParams(&chosenParams);
        shaded->SetParams(&chosenParams);
        sunlit_shaded->SetParams(&chosenParams);


        //-----------------------------------------------------------------
        // Step 11: Process each line of the climate data file
        // Climate: 
        // Day of year
        // Time  
        // Latitude 
        // Longitude 
        // altitude (m)
        // solar radiation: Radiation incident at soil surface(Wm - 2)
        // PFD, (umol s-1 m-2) (note: par *4.6; // conversion from PAR in Wm-2 to umol s-1 m-2)
        // temperature (C) 
        // co2, CO2 concentration (umol mol-1)
        // Relative humidity (%)
        // wind (m/sec)
        // LAI
        //-----------------------------------------------------------------
        while (getline(dataFile, DataLine))
        {
            if (DataLine.empty())
                continue;


            // Tokenize the climate data line using strtok_s
            char* ctx = nullptr;
            pnt = strtok_s((char*)DataLine.c_str(), pDelim, &ctx);
            DayOfYear = atof(pnt);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Time = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Lat = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Lon = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Alt = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx); 
            double SolRad = (pnt ? atof(pnt) : 0.0); //Total Radiation incident at soil surface(Wm - 2)
            pnt = strtok_s(NULL, pDelim, &ctx);
            double PFD = (pnt ? atof(pnt) : 0.0); //umol s-1 m-2 (note: par *4.6; // conversion from PAR in Wm-2 to umol s-1 m-2)
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Temp = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double CO2 = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double RH = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double Wind = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            double LAI = (pnt ? atof(pnt) : 0.0);
            pnt = strtok_s(NULL, pDelim, &ctx);
            ConstantTemperature = (pnt ? atoi(pnt) : 0);


            //-----------------------------------------------------------------
            // Step 12: Create solar geometry and radiation trans obj.
            //-----------------------------------------------------------------
            CSolar* sun = new CSolar();
            CRadTrans* light = new CRadTrans();
            sun->SetVal(DayOfYear, Time, Lat, Lon, Alt, SolRad, PFD);
            light->SetVal(*sun, LAI, chosenParams.LfAngFact);
            sunlitPFD = light->Qsl();
            shadedPFD = light->Qsh();
            sunlitLAI = light->LAIsl();
            shadedLAI = light->LAIsh();
         

            //-----------------------------------------------------------------
            // Step 13: Read SIF data if Sif yes
            //-----------------------------------------------------------------
            if (sFlag == "SIFyes")
            {
                string sifLine;
                if (!getline(sifHandle, sifLine))
                {
                    cerr << "SIF file ended before climate data finished\n";
                    break;
                }
                stringstream sifStream(sifLine);
                string tk;
                getline(sifStream, tk, ','); sifSunlit = atof(tk.c_str());
                getline(sifStream, tk, ','); sifShaded = atof(tk.c_str());
                getline(sifStream, tk, ','); ql_sunlit = atof(tk.c_str());
                getline(sifStream, tk, ','); ql_shaded = atof(tk.c_str());
                getline(sifStream, tk, ','); kdf_sunlit = atof(tk.c_str());
                getline(sifStream, tk, ','); kdf_shaded = atof(tk.c_str());
                getline(sifStream, tk, ','); PhiPS2_sunlit = atof(tk.c_str());
                getline(sifStream, tk, ','); PhiPS2_shaded = atof(tk.c_str());
                getline(sifStream, tk, ','); fesc_sun = atof(tk.c_str());
                getline(sifStream, tk, ','); fesc_sha = atof(tk.c_str());
                getline(sifStream, tk, ','); SIFtotalFactor = atof(tk.c_str());

                //SIFdata(m) is observed top - of - canopy SIF at 760 nm[mW m - 2 nm - 1 sr - 1]
                //SIF from most sensors are in units mW m - 2 nm - 1 sr - 1 
                //SIFdata from SCOPE that will be used here is  F761  [W m - 2 um - 1 sr - 1]
                // both units are same
                //https://www.sciencedirect.com/science/article/pii/S0168192324003058#sec0002
                //convert it to  μmol photons m - 2s - 1nm - 1[see eqn 15 in paper above]
                    //h = 6.626E-34 !Planck’s constant(W·s²)
                    //c = 3.0E8 !Speed of light(m / s)
                    //Na = 6.022E23 !Avogadro constant(mol⁻¹)
                    //lambda = 7.6E-7 !Wavelength(m)
                    //PI = 3.14

                SifUnitConv = (0.001 * pi * 1e6) / ((6.626e-34 * 3.0e8 * 6.022e23) / 7.6e-7);

                //fPSII(%) represents the contribution of  SIF_TOC_760 emitted by PSII is assumed to be 0.6676 for SIF Toc760
                //[see eqn 13 in paper above]
                SifUnitConv = SifUnitConv * 0.6676;
                //For fPSII, we could also use eqn in this: https://www.mdpi.com/2072-4292/16/16/2874      

            //Finally we need is total full band SIF emitted by PSII[μmol photons m - 2s - 1]
                //for that we need to multiply by f(lamda), is the conversion factor used to calculate SIF at lamda
                //nm by SIF at 760 nm, which is constant for the specific wavelength
                //we used https ://data.caltech.edu/records/xc9rx-8qs95 the FluorescenceSpectra_Greenhouse_Lightresponse.csv 
            //User can enter tehier estimate of the integrated sif instead of 136
                SifUnitConv = SifUnitConv * 136.0;
                //SifUnitConv = SifUnitConv  * SIFtotalFactor;

                sifSunlit = sifSunlit * SifUnitConv;
                sifShaded = sifShaded * SifUnitConv;
                SIF_GasEx = 1;


                //// Debug output: Print the climate and SIF data
                //cout << "Climate Data: Day=" << DayOfYear << ", Time=" << Time
                //    << ", Lat=" << Lat << ", Lon=" << Lon
                //    << ", Alt=" << Alt << ", SolRad=" << SolRad
                //    << ", PFD=" << PFD << ", Temp=" << Temp
                //    << ", CO2=" << CO2 << ", RH=" << RH
                //    << ", Wind=" << Wind << ", LAI=" << LAI << endl;
                //cout << "Species: " << runSpecies << endl;
                //cout << "SIF Data: sifSunlit=" << sifSunlit
                //    << ", sifShaded=" << sifShaded
                //    << ", ql_sunlit=" << ql_sunlit
                //    << ", ql_shaded=" << ql_shaded
                //    << ", kdf_sunlit=" << kdf_sunlit
                //    << ", kdf_shaded=" << kdf_shaded
                //    << ", PhiPS2_sunlit=" << PhiPS2_sunlit
                //    << ", PhiPS2_shaded=" << PhiPS2_shaded
                //    << ", fesc_sun=" << fesc_sun
                //    << ", fesc_sha=" << fesc_sha << endl;

                if (sunlitShadedFlag == "SunlitShadedYes")
                {
                    // Pass the SIF and climate parameters to the gas exchange objects.
                    sunlit->SetVal(sunlitPFD, Temp, CO2, RH, Wind, sunlitLAI, Pressure, ConstantTemperature,
                        SIF_GasEx, sifSunlit, ql_sunlit, kdf_sunlit, PhiPS2_sunlit, fesc_sun);
                    shaded->SetVal(shadedPFD, Temp, CO2, RH, Wind, shadedLAI, Pressure, ConstantTemperature,
                        SIF_GasEx, sifShaded, ql_shaded, kdf_shaded, PhiPS2_shaded, fesc_sha);
                }
                else
                {
                   // No sunlit/shaded specific- all read as sunlit (to have minimal changes in the input

				   //SCOPE can have each component as sunlit shaded (e.g., fesc, ql, kdf, etc.)
                   // ql_sunlit = (ql_shaded * shadedLAI + ql_sunlit * sunlitLAI);
                   //sunlit_shaded->SetVal(sunlitPFD + shadedPFD, Temp, CO2, RH, Wind, LAI, Pressure, ConstantTemperature,
                   // SIF_GasEx, sifSunlit * sunlitLAI + sifShaded * shadedLAI, ql_sunlit, (kdf_sunlit + kdf_shaded) / 2, PhiPS2_sunlit, fesc_sun * sunlitLAI + fesc_sha * shadedLAI);
                   // Note: The below is when canopy level measurements (e.g, field level, EC and SIF tower)
                   // No sunlit/shaded specific- SIF canopy readed as sunlit (to have minimal changes in the input)
    
                    if (runSpecies == "maize")
                    {
                       
                        ql_sunlit = 0.83 * exp(-0.00063 * sunlitPFD);  // Table S3  in Han et al. (2022))
                        ql_shaded = 0.83 * exp(-0.00063 * shadedPFD);
                        ql_sunlit_shaded = (ql_sunlit * sunlitLAI + ql_shaded * shadedLAI) / LAI;
                        kdf_sunlit_shaded = 13.5;
                        PhiPS2_sunlit_shaded = 0.83;
                    }
                    
                    if (runSpecies == "soybean")
                    {
                        // For other species, we can use the provided values directly.
                        ql_sunlit = 0.80 * exp(-0.00095 * sunlitPFD); //Table S3 in Han et al. (2022))
                        ql_shaded = 0.80 * exp(-0.00095 * shadedPFD);
                        ql_sunlit_shaded = (ql_sunlit * sunlitLAI + ql_shaded * shadedLAI) / LAI; 
                        kdf_sunlit_shaded = 13;
                        PhiPS2_sunlit_shaded = 0.81;


                    }

 
                  
                    sunlit_shaded->SetVal(sunlitPFD+shadedPFD, Temp, CO2, RH, Wind, LAI, Pressure, ConstantTemperature,
                    SIF_GasEx, sifSunlit, ql_sunlit_shaded, kdf_sunlit_shaded, PhiPS2_sunlit_shaded, fesc_sun);
                }
               
            }
            else
            {
                // If no SIF data are provided, set all SIF parameters to zero.
                sifSunlit = 0.0;
                sifShaded = 0.0;
                ql_sunlit = 0.0;
                ql_shaded = 0.0;
                kdf_sunlit = 0.0;
                kdf_shaded = 0.0;
                PhiPS2_sunlit = 0.0;
                PhiPS2_shaded = 0.0;
                fesc = 0.0;
                SIF_GasEx = 0;
                if (sunlitShadedFlag == "SunlitShadedYes")
                {
                    sunlit->SetVal(sunlitPFD, Temp, CO2, RH, Wind, sunlitLAI, Pressure, ConstantTemperature,
                        SIF_GasEx, 0.0, 0.0, 0.0, 0.0, 0.0);
                    shaded->SetVal(shadedPFD, Temp, CO2, RH, Wind, shadedLAI, Pressure, ConstantTemperature,
                        SIF_GasEx, 0.0, 0.0, 0.0, 0.0, 0.0);
                }
                else
                { 
                sunlit_shaded->SetVal(PFD, Temp, CO2, RH, Wind, LAI, Pressure, ConstantTemperature,
                    SIF_GasEx, 0.0, 0.0, 0.0, 0.0, 0.0);
                }

                // Debug output: Print the climate and SIF data
                //cout << "Climate Data: Day=" << DayOfYear << ", Time=" << Time
                //    << ", Lat=" << Lat << ", Lon=" << Lon
                //    << ", Alt=" << Alt << ", SolRad=" << SolRad
                //    << ", PFD=" << PFD << ", Temp=" << Temp
                //    << ", CO2=" << CO2 << ", RH=" << RH
                //    << ", Wind=" << Wind << ", LAI=" << LAI << endl;
                //cout << "Species: " << runSpecies << endl;
                //cout << "SIF Data: sifSunlit=" << sifSunlit
                //    << ", sifShaded=" << sifShaded
                //    << ", ql_sunlit=" << ql_sunlit
                //    << ", ql_shaded=" << ql_shaded
                //    << ", kdf_sunlit=" << kdf_sunlit
                //    << ", kdf_shaded=" << kdf_shaded
                //    << ", PhiPS2_sunlit=" << PhiPS2_sunlit
                //    << ", PhiPS2_shaded=" << PhiPS2_shaded
                //    << ", fesc_sun=" << fesc_sun
                //    << ", fesc_sha=" << fesc_sha << endl;

            }



            //-----------------------------------------------------------------
            // Step 14: model outputs from the sunlit/shaded gas exchange object (need to check the details)
            //-----------------------------------------------------------------
            // Model output variables
            double Anet = 0.0, Agross = 0.0, VPD = 0.0, LeafTemperature = 0.0;
            double BoundaryLayerConductance = 0.0, InternalCO2 = 0.0;
            double Respiration = 0.0, StomatalConductance = 0.0, Transpiration = 0.0;
            double AjLeaf = 0.0, AcLeaf = 0.0, AvLeaf = 0.0, ApLeaf = 0.0,Jval = 0.0;
            double GammaVal = 0.0, GammaC4Val = 0.0, Cs_stom_conductanceVal = 0.0;


            // this is for sunlit/shaded condition
            if (sunlitShadedFlag == "SunlitShadedYes")
            {
                Anet = (sunlit->get_ANet() * sunlitLAI + shaded->get_ANet() * shadedLAI); //(umol CO2 m - 2 s - 1)
                Agross = (sunlit->get_AGross() * sunlitLAI + shaded->get_AGross() * shadedLAI); //(umol CO2 m - 2 s - 1)
                VPD = sunlit->get_VPD();//vapor pressure deficit (kpa)
                LeafTemperature = sunlit->get_LeafTemperature();//leaf temperature (C)
                BoundaryLayerConductance = sunlit->get_BoundaryLayerConductance(); //boundary layer conductance(mol m - 2 s - 1)
                InternalCO2 = (sunlit->get_Ci() * sunlitLAI +shaded->get_Ci() * shadedLAI)/(sunlitLAI+ shadedLAI);//internal CO2 concentration (umol mol-1)
                Respiration = (sunlit->get_Respiration() * sunlitLAI + shaded->get_Respiration() * shadedLAI); //(umol CO2 m-2 s-1)
                StomatalConductance = __max(0, ((sunlit->get_StomatalConductance() * sunlitLAI + shaded->get_StomatalConductance() * shadedLAI) ));//average stomatal conductance to water vapor (mol m-2 s-1)
                Transpiration = (sunlit->get_Transpiration() * sunlitLAI + shaded->get_Transpiration() * shadedLAI);//(umol H2O m-2 s-1)

                AjLeaf = (sunlit->get_Aj() * sunlitLAI + shaded->get_Aj() * shadedLAI); //(umol CO2 m - 2 s - 1)
                AcLeaf = (sunlit->get_Ac() * sunlitLAI + shaded->get_Ac() * shadedLAI); //(umol CO2 m - 2 s - 1)
                AvLeaf = (sunlit->get_Av() * sunlitLAI + shaded->get_Av() * shadedLAI); //(umol CO2 m - 2 s - 1)
                ApLeaf = (sunlit->get_Ap() * sunlitLAI + shaded->get_Ap() * shadedLAI); //(umol CO2 m - 2 s - 1)
                Jval = (sunlit->get_J() * sunlitLAI + shaded->get_J() * shadedLAI);
                sunlitPFD = sunlitPFD * sunlitLAI;
                shadedPFD = shadedPFD * shadedLAI;
                GammaVal = (sunlit->get_GammaValueC3_C4() * sunlitLAI + shaded->get_GammaValueC3_C4() * shadedLAI);
                GammaC4Val = sunlit->get_GammaC4();               // sunlit and shaded can have same gammaC4
                Cs_stom_conductanceVal = sunlit->get_Cs_stom_conductance();

            }
            else
            {

				Anet = sunlit_shaded->get_ANet(); //(umol CO2 m - 2 s - 1) ! SIF is already at canopy level, so no need to multiply by LAI
                Agross = sunlit_shaded->get_AGross(); //(umol CO2 m - 2 s - 1)
                VPD = sunlit_shaded->get_VPD();//vapor pressure deficit (kpa)
                LeafTemperature = sunlit_shaded->get_LeafTemperature();//leaf temperature (C)
                BoundaryLayerConductance = sunlit_shaded->get_BoundaryLayerConductance(); //boundary layer conductance(mol m - 2 s - 1)
                InternalCO2 = (sunlit_shaded->get_Ci());//internal CO2 concentration (umol mol-1)
                Respiration = (sunlit_shaded->get_Respiration()); //(umol CO2 m-2 s-1)
                StomatalConductance = __max(0, (sunlit_shaded->get_StomatalConductance()))/LAI;//average stomatal conductance to water vapor (mol m-2 s-1)
                Transpiration = (sunlit_shaded->get_Transpiration());//(umol H2O m-2 s-1)

                AjLeaf = sunlit_shaded->get_Aj(); //(umol CO2 m - 2 s - 1)
                AcLeaf = sunlit_shaded->get_Ac(); //(umol CO2 m - 2 s - 1)
                AvLeaf = sunlit_shaded->get_Av(); //(umol CO2 m - 2 s - 1)
                ApLeaf = sunlit_shaded->get_Ap(); //(umol CO2 m - 2 s - 1)
                Jval = sunlit_shaded->get_J()/LAI;
                GammaVal = sunlit_shaded->get_GammaValueC3_C4();
                GammaC4Val = sunlit_shaded->get_GammaC4();                // sunlit and shaded can have same gammaC4
                Cs_stom_conductanceVal = sunlit_shaded->get_Cs_stom_conductance();

            }
            // if no sunlit shaded condition

           



            //-----------------------------------------------------------------
            // Step 15: output file.
            //-----------------------------------------------------------------
            fprintf(pFile, "%s, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, "
                "%8.2f, %8.2f, %8.2f, %d, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,%8.2f, %8.2f, %8.2f, %8.2f,"
                "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,%8.2f,%8.2f,%8.2f,%8.2f, "
                "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,%8.2f,%8.2f,%8.2f,%8.2f, "
                "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,%8.2f,%8.2f,%8.2f,%8.2f, "
                "%8.2f, % 8.2f,  % 8.2f,% 8.2f,%8.2f\n",


                runSpecies.c_str(),   // Species used for this run.
                DayOfYear,            // Climate parameter: Day of Yea
                Time,              // Climate parameter: Time
                Lat,                  // Climate parameter: Latitude
                Lon,                  // Climate parameter: Longitude
                Alt,                  // Climate parameter: Altitude
                SolRad,                 // Climate parameter: Solar Radiation
                PFD,                  // Climate parameter: PFD
                Temp,                 // Climate parameter: Temperature
                CO2,                 // Climate parameter: CO2 concentration.(umol mol-1)
                RH,                // Climate parameter: Relative Humidity
                Wind,              // Climate parameter: Wind speed.(m/s)
                LAI,               // Climate parameter: LAI.
                (int)ConstantTemperature, // Climate parameter: Constant Temperature flag.
                (sifSunlit <= 0.0 || SifUnitConv <= 0.0 || std::isnan(sifSunlit) || std::isnan(SifUnitConv) ? 0.0 : sifSunlit / 1.81),             // SIF parameter: sifvalue /1.81 prints in mW/m2/nm/s.
                (sifShaded <= 0.0 || SifUnitConv <= 0.0 || std::isnan(sifShaded) || std::isnan(SifUnitConv) ? 0.0 : sifShaded / 1.81),             // SIF parameter: sifvalue.
                ql_sunlit,            // SIF parameter: ql_sunlit.
                ql_shaded,            // SIF parameter: ql_shaded
                kdf_sunlit,           // SIF parameter: kdf_sunlit
                kdf_shaded,           // SIF parameter: kdf_shaded.
                PhiPS2_sunlit,        // SIF parameter: PhiPS2 _sunlit
                PhiPS2_shaded,          // SIF parameter: PhiPS2_shded.
                fesc_sun+fesc_sha,                 // SIF parameter: fesc
                Anet,                 // Model output: Net photosynthesis
                Agross,               // Model output: Gross photosynthesis
                VPD,                  // Model output: Vapor Pressure Deficit
                LeafTemperature,      // Model output: Leaf Temperature.
                BoundaryLayerConductance,    // Model output: Boundary Layer Conductance
                InternalCO2,                 // Model output: Internal CO2
                Respiration,                 // Model output: Respiration
                StomatalConductance,         // Model output: Stomatal Conductance
                Transpiration,                 // Model output: Transpiration
                Jval,
                sunlitPFD,
                shadedPFD,
                sunlitLAI,
                shadedLAI,
                //for debug purpose 
                GammaVal,
                GammaC4Val,                // sunlit and shaded can have same gammaC4
                Cs_stom_conductanceVal,
                InternalCO2,
                AjLeaf,
                AcLeaf,
                AvLeaf,
                ApLeaf,
                sun->GetAzimuth(),
                sun->GetDayLength(),
                sun->GetDeclination(),
                //sun->GetFracNIRDirect(),
                sun->GetFracPARDiffuse(),
                sun->GetFracPARDirect(),
                sun->GetNIR(),
                //sun->GetNIRDiffuse(),
                //sun->GetNIRDirect(),
                sun->GetNIRFraction(),
                sun->GetNIRTotal(),
                sun->GetPAR(),
                sun->GetSolarNoon(),
                sun->GetSunrise(),
                sun->GetSunset(),
                sun->GetSinElevation(),
                sun->GetSolarElevation(),
                sun->GetPotentialSolarDiffuse(),
                sun->GetPotentialSolarTotal(),
                sun->GetPotentialSolarDirect(),
                sun->GetSolarRadiation(),
                sun->GetPARFraction(),
                sun->GetPFD(),
                sun->GetPFDDiffuse(),
                sun->GetPFDDirect(),
                sun->GetPFDTotal(),
                light->GetZenith(), 
                sun->GetPotentialPARDirect(),
                sun->GetPotentialPARDiffuse(),
                sun->GetPotentialNIRDirect(),
                sun->GetPotentialNIRDiffuse()
              

            );


            // Clean temporary solar and radiation objs.
            delete sun;
            delete light;
        } // End of climate data loop.


        //-----------------------------------------------------------------
        // Step 16: Clean gas exchange obj's and close files for the current run.
        //-----------------------------------------------------------------
        delete sunlit;
        delete shaded;

        dataFile.close();
        if (sFlag == "SIFyes")
            sifHandle.close(); 
        fclose(pFile);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Total execution time: " << elapsed.count() << " seconds." << std::endl;
    } // End loop of run file


    runFile.close();



    cout << "Done." << endl;
   

    return 0;
}





