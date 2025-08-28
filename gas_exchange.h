﻿/*! @file 
      function headers for calculations 
* @mainpage Photosynthesis Module Overview 
*  @section intro_sec Introduction
* This is a coupled model of leaf gas-exchange and photosynthesis-stomatal conductance-energy balance for a maize leaf. 
* The C++ GasExchange class encapsulates all the calculations to estimate 
* CO2 assimilation rate,  stomatal conductance, leaf temperature, and transpiration for a square meter of leaf.
* The calculations are based on von Caemmerer(2000) C4 model, Kim and Leith (2003) C3 rose model, BWB stomatal
* conductance(Ball et al., 1987) and energy balance model as described in Campbell and Norman (1998) 


@authors Soo-Hyung Kim, Univ. of Washington
@authors Dennis Timlin, USDA-ARS
@authors David Fleisher, USDA-ARS
@version 1.0
@date June 2016
@license This project is released under the GNU Public License

 <b>Bibliography </b>
\li Ball J.T., Woodrow I.E., Berry J.A. 1987. A model predicting stomatal conductance and its contribution to the control of photosynthesis  under different environmental conditions. In: Biggens J, ed. Progress in photosynthesis research. The Netherlands: Martinus Nijhoff Publishers. 
\li Campbell, G.S., and J.M. Norman. 1998. The light environment of plant canopies. An Introduction to Environmental Biophysics. Springer, New York. pp: 247-281. 
\li Kim, S.-H., and J.H. Lieth. 2003. A coupled model of photosynthesis, stomatal conductance and transpiration for a rose leaf (Rosa hybrida L.). Ann. Bot. 91:771-781. 
\li Kim, S.-H., D.C. Gitz, R.C. Sicher, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2007. Temperature dependence of growth, development, and photosynthesis in maize under elevated CO2. Env. Exp. Bot. 61:224-236. 
\li Kim, S.-H., R.C. Sicher, H. Bae, D.C. Gitz, J.T. Baker, D.J. Timlin, and V.R. Reddy. 2006. Canopy photosynthesis, evapotranspiration, leaf nitrogen, and transcription  
\li von Caemmerer, S., 2000. Biochemical models of leaf photosynthesis. CSIRO Publishing, Collingwood, Australia. \n
*
*
\n
*
*
*
*
*
<b> Source Code </b>
<BR>
// GitHub https://github.com/ARS-CSGCL-DT
*
*
*
*
\n
\n
\n
\n
\n

*/


/*!

<P>
* \page page1 A flow chart of the model
* \image rtf Flowchart.tif "Flow Chart for model" width=8cm
* \image html flowchart.tif "Flow Chart for model" width=8cm
</P>
\n
\n
\n
\n
\n
*/


#include "stdafx.h"

using namespace std;
namespace photomod   //Photosynthesis Model NameSpace
{
	/*!<p><!-- pagebreak --></p>
	  \page page4  Model Documentation 
	   \section model_sec Description of Modules
	   These pages describe the model in detail */


/*! \class CGasExchange 
* \brief Class for gas exchange calculations\n
* This class simulates gas exchange in plant leaves for C3 and C4 plants. \n
* \par Usage
	- Use <b>SetParams</b> to initialize the GasExchange object with parameters for a specific variety of plant
	- Use <b>SetVal</b> to pass environmental variables for a simulation and return a structure with output. \n
* See \ref Interface for details on input and output	
    

	*/

 
	class CGasExchange
	{


	public:
		CGasExchange();

		~CGasExchange(void);



		/*!    \struct tparms gas_exchange.h
			   \brief Structure to hold plant parameters for the photosynthesis model. \n
			   Some parameters are specific for C3 or C4 type Plants \n
		*/

		/*!
		//                parameter        description
		@param		ID		1		Name of plant (string)
		@param 		species	2 		Species Name  (string)
		@param 		type	3			'C3' or 'C4' (string)
		@param 		Vcm25	4		Photosynthetic Rubisco Capacity at 25C (umol m-2 s-1)
		@param 		Jm25	5		Potential Rate of electron transport at 25C  (umol m-2 s-1)
		@param 		Vpm25	6		C4 Carboxylation rate at 25C (C4, umol m-2 s-1)
		@param 		TPU25	7		Rate if Triose Phosphate Utilization at 25C (C3, umol m-2 s-1)
		@param 		RD25	8		Mitochondrial respiration in the light at 25C (umol m-2 s-1)
		@param 		Theta   9		Initial slope of CO2 response (umol m2 s-1) - de Pury (1997)
		@param 		EaVc   10		Activation energy for Arrhenius function used to calculate temperature dependence for Vcmax (kJ mol-1)
		@param 		Eaj    11		Activation energy for Arrhenius function used to calculate temperature dependence for J (kJ mol-1)
		@param 		Hj     12		Curvature parameter of the temperature dpendence of Jmax (kJ mol-1)
		@param		Sj	   13		Electron transport temperature response parameter for Jmax (J mole-1 K-1)
		@param      Hv     14       Curvature parameter of the temperature dependence of Vcmax (J mole-1)
		@param 	    EaVp   15		Activation energy for Arrhenius function used to calculate temperature dependence for Vpmax (kJ mol-1)
		@param 	    Sv	   16		Electron transport temperature response parameter for Vcmax (J mole-1 K-1)
		@param 		EAP	   17		Activation energy for Arrhenius function used to calculate temperature dependence for TPU (kJ mol-1)
		@param 		EAR	   18	 	Activation energy for Arrhenius function used to calculate temperature dependence for respiration (kJ mol-1)
		@param 		g0	   19		Minimum stomatal conductance to water vapor at the light compensation point in the BWB model (mol m-2 s-1)
		@param 	    g1	   20   	Empirical coefficient for the sensitivity of StomatalConductance to A, Cs and hs in BWB model (no units)
		@param 		StomRatio	21	Stomatal Ratio (fraction)
		@param  	LfWidth     22	Leaf Width (m)
		@param 		LfAngFact	23	Leaf Angle Factor
		@param 		Remark		24	Text
		*/

		struct tParms    //*!< holds parameters for the model
		{
			string   ID;
			std::string species;
			std::string Type;

			double
				Vcm25,
				Jm25,
				Vpm25,
				TPU25,
				Rd25,
				Theta,
				EaVc,
				Eaj,
				Hj,
				Sj,
				Hv,
				EaVp,
				Sv,
				Eap,
				Ear,
				g0,
				g1,
				stomaRatio,
				LfWidth,
				LfAngFact;
			std::string Remark;
		}sParms; //!< Variable of type tParms 


		//These functions return results of calculations 

		double get_ANet() { return AssimilationNet; }      //!< return net photosynthesis (umol CO2 m-2 s-1) 
		double get_AGross() { return AssimilationGross; }  //!< return gross photosynthesis  (umol CO2 m-2 s-1)
		double get_Transpiration() { return Transpiration; }   //!< return transpiration rate (umol H2O m-2 s-1)
		double get_LeafTemperature() { return Tleaf; } //!< return leaf temperature (C)
		double get_Ci() { return Ci; }        //!< return internal CO2 concentration (umol mol-1)
		double get_StomatalConductance() { return StomatalConductance; }    //!< return stomatal conductance to water vapor (mol m-2 s-1)
		double get_BoundaryLayerConductance() { return BoundaryLayerConductance; }    //!< return boundary layer conductance (mol m-2 s-1)
		double get_Respiration() { return DarkRespiration; }   //!< return respiration rate (umol CO2 m-2 s-1)
		double get_VPD() { return VPD; }         //!< return vapor pressure deficit (kpa)
		void  SetParams(tParms* sParms);              //!<used to pass structure with parameters from the interface to the calculator
		void SetVal(double PhotoFluxDensity, double Tair, double CO2, double RH,
			double wind, double LAI, double Press, bool ConstantTemperature,
			int SIF_GasEx, double sifvalue, double qL, double  kdf, double phiPS2max, double fesc);
		//!< sets input values for calculations for a particular set of environmental variables
		//! Debug purpose
		double get_GammaC4()
		{
			return (PlantType == "C3") ? 0.0 : Gamma_C4;
		}

		double get_Cs_stom_conductance() { return Cs_stom_conductance; }
		double get_GammaValueC3_C4() { return GammaValue; }
		double get_J() { return Jvalue; }

		double get_Aj() { return AjValue; }
		double get_Ac() { return AcValue; }
		double get_Av() { return AvValue; }
		double get_Ap() { return ApValue; }
	protected:
		void GasEx();  //!<Main module to calculate gas exchange rates
		void PhotosynthesisC3(double Ci);  //!<Function to calculate C3 photosynthesis
		void PhotosynthesisC4(double Ci);  //!<Function to calculate C4 photosynthesis
		void EnergyBalance();     //!<calculates leaf temperature and transpiration
		double SearchCi(double CO2i);          //!< called iterively to find optimal internal CO2 concentration returns optimal internal CO2 concentration (CO2i) 
		double EvalCi(double Ci);   //!<Calls photosynthesis modules to evaluate Ci dependent calculations returns difference between old Ci and new Ci
		double CalcStomatalConductance();             //!< Stomatal conductance (mol m-2 s-1)
		double CalcTurbulentVaporConductance();             //!<  Conductance for turbulant vapor transfer in air - forced convection (mol m-2 s-1)
		double Es(double Temperature);      //!< Saturated vapor pressure at given temperature. kPa
		double Slope(double Temperature);  //!<  Slope of the vapor pressure curve
		double QuadSolnUpper (double a, double b, double c ); //!< Upper part of quadratic equation solution
		double QuadSolnLower (double a, double b, double c ); //!< Lower part of quadratic equation solution
		double minh(double fn1,double fn2,double theta2); //!< Hyperbolic minimum
		double get_CiCaRatio()  {return Ci_Ca;} //!< Ratio between internal and atmospheric CO2

		double  PhotoFluxDensity;  //!< Photosynthetic Flux Density (umol photons m-2 s-1 
		double 	R_abs; //!< Absorbed incident radiation (watts m-2)        
		double 	Tair;  //!< Air temperature at 2m, (C) 
		double 	CO2;   //!< CO2 concentration (umol mol-1 air) 
		double 	RH;   //!<  Relative Humidity (%, i.e., 80) 
		double 	wind; //!<  Windspeed at 2 meters (km s-1) 
		double 	width; //!< Leaf width (m) 
		double 	Press;  //!<  Air pressure (kPa) 
		
		double Anet = 0.0, Agross = 0.0,  LeafTemperature = 0.0;
		double InternalCO2 = 0.0, Respiration = 0.0;
	    double AjLeaf = 0.0, AcLeaf = 0.0, AvLeaf = 0.0, ApLeaf = 0.0;
	    double GammaVal = 0.0, GammaC4Val = 0.0, Cs_stom_conductanceVal = 0.0;



		int SIF_GasEx;//!< switch for SIF yes or No
		double sifvalue;//!< total full band SIF emitted by PSII [μmol photons m-2s-1]
		double qL;//!< fraction of open PSII reaction centers
		double kdf; //!<  kDF=kD/kF. kD and kF are the rate constants of constitutive thermal dissipation and fluorescence
		double phiPS2max;//!< maximum photochemical quantum yield of PSII
		double fesc;//!< canopy escape probability of SIF photons
	
		// These variables hold the output of calculations 
		double AssimilationNet,    //!< Net photosynthesis (umol CO2 m-2 s-1)
			AssimilationGross, //!< Gross photosynthesis (umol CO2 m-2 s-1) (Adjusted for respiration)
			Transpiration,     //!< Transpiration mol H2O m-2 s-1 
			Tleaf,  //!< Leaf temperature C 
			Ci,     //!< Internal CO2 concentration umol mol-1 
			StomatalConductance,     //!< Stomatal conductance umol m-2 s-1 
			BoundaryLayerConductance,    //!< Boundary layer conductance umol m-2 s-1 
			DarkRespiration,    //!< Plant respiration    umol m-2 s-1 
			VPD,    //!< Vapor Pressure Density, kPa */
			Ci_Ca,//!< Ratio of internal to external CO2, unitless
			GammaValue, Jval;

		//for debug gamma (from c3), GammaStar (from c4), Gamma (from c4) and Cs (from StomatalConductance)
		double Gamma_C4;
		double Cs_stom_conductance;
		double Jvalue;
		double AcValue;
		double AjValue;
		double AvValue;
		double ApValue;
		double errTolerance; /*!< Error tolerance for iterations */
		double eqlTolerance; /*!< Equality tolerance */

		int iter_total;      //!< Total number of iterations */
		///  @param PlantType string that holds the type of plant, C3 or C4

		std::string PlantType;       //!< For C3 or C4 designation */
		bool ConstantLeafTemperature; //!< if true, uses constant temperature - if true, does not solve for leaf temperature */
		int iter1, iter2;         //!< Iteration counters
		int  iter_Ci;   /*!< Iteration value for Ci umol mol-1, internal CO2 concentration */
		bool isCiConverged; /*!< True if Ci (internal CO2 concentration) iterations have converged */
		inline double Square(double a) { return a * a; } /*!< Squares number */ 
		inline double Min(double a, double b, double c) {return (__min(__min(a,b),c));} /*!< Finds minimum of three numbers */


	};


} //end namespace photomod
 
