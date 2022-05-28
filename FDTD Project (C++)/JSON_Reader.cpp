/**
 * @file JSON_Reader.cpp
 * @author FDTD-Lab (https://github.com/FDTD-Lab)
 * @version 0.1
 * @date 2022-01-19
 *
 * @brief JSON Reader file, housing function to read in JSON simulation
 * control files, including custom DFT controls.
 *
 * @details This file allows for the simulation controls to be given
 * via JSON format. Has external dependencies. JSON inputs are assume to
 * be given wrapped in a list ( mean '[ { } ]' ) structure. Can also set
 * custom DFT monitors.
 *
 * @copyright Copyright (c) 2022
 */
#include <stdio.h>
#include <json-c/json.h>
#include <cstring>
#include "externs.h"


/**
 * @brief Reads in the JSON simulation control file and unpack into 'pids'
 *
 * @param strPath [in] path to auxilliary files to be read (e.g. chi tensors)
 * @param pids global simulation control structure
 * @param FailureType [out] location to store error, should one occur
 * @return true when error has occured
 * @return false when no error was raised
 */
bool ReadJSONDataFile(const string &strPath, INI_DATA_STRUCT *&pids, FAILURE_TYPES &FailureType)
{
	fstream pfin;
	string strS;
	FILE *fp;
	char buffer[102422];

	int NonLoc, G_D, i;
	string NL_ANTENNA_MAT, NL_FIRST_MAT, NL_MIDDLE_MAT, NL_SECOND_MAT; // Izzat
	REAL E_FIELD_MAX;

	// JSON object variables for the reading in and unpacking
	struct json_object *jsonControls;	/** JSON Controls object */
	struct json_object *dummyValue; 	/** placeholder variable for reading in single values */
	struct json_object *dummyArray; 	/** placeholder variable for reading in arrays */

	// Read in JSON file to buffer
	string controlPathJSON = strPath + "/pphinfoini.json";
	fp = fopen(controlPathJSON.c_str(), "r");
	size_t successReadCount = fread(buffer, 102422, 1, fp);
	fclose(fp);

	// Create JSON object from buffer
	jsonControls = json_tokener_parse(buffer);

	// -----------------------------------------------------------------------
	// 							SIMULATION CONTROLS
	// -----------------------------------------------------------------------

	// Unpack "Domain Size" --> nCells in each direction
	json_object_object_get_ex(jsonControls, "Domain Size", &dummyArray);
	pids->NRCELLSTOTX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
	pids->NRCELLSTOTY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
	pids->NRCELLSTOTZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));

	// Unpack "Domain Decomposition" --> nProcesses in each direction
	json_object_object_get_ex(jsonControls, "Domain Decomposition", &dummyArray);
	pids->XPROCSSIZE = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
	pids->YPROCSSIZE = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
	pids->ZPROCSSIZE = json_object_get_int(json_object_array_get_idx(dummyArray, 2));

	// Number of time steps to take for simulation
	json_object_object_get_ex(jsonControls, "Number of Time Steps", &dummyValue);
	pids->NRMAXTIMESTEPS = (dummyValue == NULL) ? 100000 : json_object_get_int(dummyValue);

	// Unpack "PML Box" --> nCells_PML at each wall, each direction
	json_object_object_get_ex(jsonControls, "PML Box", &dummyArray);
	if (dummyArray == NULL)
	{
		// Default PML size
		pids->NRCELLSCPMLX = 20;
		pids->NRCELLSCPMLY = 20;
		pids->NRCELLSCPMLZ = 20;
	}
	else
	{
		pids->NRCELLSCPMLX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		pids->NRCELLSCPMLY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		pids->NRCELLSCPMLZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	// Unpack "TFSF Box" --> nCells for TFSF Zone, each wall, each direction
	json_object_object_get_ex(jsonControls, "TFSF Box", &dummyArray);
	if (dummyArray == NULL)
	{
		// Default TFSF size
		pids->NRINFLATEX = 11;
		pids->NRINFLATEY = 11;
		pids->NRINFLATEZ = 11;
	}
	else
	{
		pids->NRINFLATEX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		pids->NRINFLATEY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		pids->NRINFLATEZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	// Unpack "Scattering Box" --> nCells from boundary for scatter calc
	json_object_object_get_ex(jsonControls, "Scattering Box", &dummyArray);
	if (dummyArray == NULL)
	{
		// Default Scattered Field calculation
		pids->NRQSCAX = pids->NRINFLATEX - 5;
		pids->NRQSCAY = pids->NRINFLATEY - 5;
		pids->NRQSCAZ = pids->NRINFLATEZ - 5;
	}
	else
	{
		pids->NRQSCAX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		pids->NRQSCAY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		pids->NRQSCAZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	// Unpack "Total Box" --> nCells from boundary for total calc
	json_object_object_get_ex(jsonControls, "Total Box", &dummyArray);
	if (dummyArray == NULL)
	{
		// Default Scattered Field calculation
		pids->NRQTOTX = pids->NRINFLATEX + 5;
		pids->NRQTOTY = pids->NRINFLATEY + 5;
		pids->NRQTOTZ = pids->NRINFLATEZ + 5;
	}
	else
	{
		pids->NRQTOTX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		pids->NRQTOTY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		pids->NRQTOTZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	// Unpack "Step Size" --> cell dimensions in each direction
	json_object_object_get_ex(jsonControls, "Step Size", &dummyArray);
	pids->DX = json_object_get_double(json_object_array_get_idx(dummyArray, 0));
	pids->DY = json_object_get_double(json_object_array_get_idx(dummyArray, 1));
	pids->DZ = json_object_get_double(json_object_array_get_idx(dummyArray, 2));

	// Origin of simulation domain
	json_object_object_get_ex(jsonControls, "Center Position", &dummyArray);
	if (dummyArray == NULL)
	{
		// Default center positions
		pids->CENTERPOSITIONX = int(0.5 * (pids->NRCELLSTOTX - 2 * pids->NRCELLSCPMLX));
		pids->CENTERPOSITIONY = int(0.5 * (pids->NRCELLSTOTY - 2 * pids->NRCELLSCPMLY));
		pids->CENTERPOSITIONZ = int(0.5 * (pids->NRCELLSTOTZ - 2 * pids->NRCELLSCPMLZ));
	}
	else
	{
		pids->CENTERPOSITIONX = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		pids->CENTERPOSITIONY = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		pids->CENTERPOSITIONZ = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	// Periodic boundary flag in X direction
	json_object_object_get_ex(jsonControls, "Periodic X", &dummyValue);
	pids->PBC_X = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Periodic boundary flag in Z direction
	json_object_object_get_ex(jsonControls, "Periodic Z", &dummyValue);
	pids->PBC_Z = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// -----------------------------------------------------------------------
	// 							SOURCE CONTROLS
	// -----------------------------------------------------------------------

	// Determine Signal type via integer or string
	json_object_object_get_ex(jsonControls, "Signal Type", &dummyValue);
	if (dummyValue == NULL)
	{
		// Default Signal type
		pids->SIGNAL_TYPE = 1;
	}
	else
	{
		if (strcmp(json_object_get_string(dummyValue), "Raised Cosine") == 0)
		{
			pids->SIGNAL_TYPE = 1;
		}
		if (strcmp(json_object_get_string(dummyValue), "Gaussian") == 0)
		{
			pids->SIGNAL_TYPE = 6;
			CENTERFREQ = 1000e-9; // default
			PULSEWIDTH = 5e-15;	  // default
			json_object_object_get_ex(jsonControls, "Center Wavelength", &dummyValue);
			CENTERFREQ = json_object_get_double(dummyValue);
			json_object_object_get_ex(jsonControls, "Pulse Width", &dummyValue);
			PULSEWIDTH = json_object_get_double(dummyValue);
		}
		if (strcmp(json_object_get_string(dummyValue), "Read in Dipole") == 0)
		{
			pids->SIGNAL_TYPE = 100;
		}
	}

	// wMax of Raised Cosine
	json_object_object_get_ex(jsonControls, "Maximum Frequency", &dummyValue);
	pids->MAXFREQUENCY = (dummyValue == NULL) ? 2.0e15 : json_object_get_double(dummyValue);

	// Source type? Like, TFSF box, Gaussian, etc
	json_object_object_get_ex(jsonControls, "TFSF", &dummyValue);
	pids->TFSF = (dummyValue == NULL) ? 1 : json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "Polarization Angle", &dummyValue);
	pids->POLARIZATION_ANGLE_1 = (dummyValue == NULL) ? 90 : json_object_get_double(dummyValue);
	
	json_object_object_get_ex(jsonControls, "Polarization Angle 2", &dummyValue);
	pids->POLARIZATION_ANGLE_2 = (dummyValue == NULL) ? 90 : json_object_get_double(dummyValue);

	// Electric field peak
	json_object_object_get_ex(jsonControls, "Maximum Field", &dummyValue);
	E_FIELD_MAX = (dummyValue == NULL) ? 1.0 : json_object_get_double(dummyValue);

	// Center Frequency (Gaussian Only)
	json_object_object_get_ex(jsonControls, "Center Wavelength", &dummyValue);
	CENTERFREQ = json_object_get_double(dummyValue);

	// -----------------------------------------------------------------------
	// 							MONITOR CONTROLS
	// -----------------------------------------------------------------------

	// Number of wavelengths to monitor
	json_object_object_get_ex(jsonControls, "Number of Wavelengths", &dummyValue);
	pids->NRANFREQUENCIES = (dummyValue == NULL) ? 1 : json_object_get_int(dummyValue);

	// Minimum wavelength to monitor
	json_object_object_get_ex(jsonControls, "Minimum Wavelength", &dummyValue);
	pids->STARTLAMBDA = (dummyValue == NULL) ? 100e-9 : json_object_get_double(dummyValue);

	// Maximum wavelength to monitor
	json_object_object_get_ex(jsonControls, "Maximum Wavelength", &dummyValue);
	pids->STOPLAMBDA = (dummyValue == NULL) ? 800e-9 : json_object_get_double(dummyValue);

	// -----------------------------------------------------------------------
	// 							GEOMETRY CONTROLS
	// -----------------------------------------------------------------------

	// Material index of entry layer
	json_object_object_get_ex(jsonControls, "First Medium", &dummyValue);
	pids->FIRST_MEDIUM = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Material index of middle layer
	json_object_object_get_ex(jsonControls, "Middle Medium", &dummyValue);
	pids->MIDDLE_MEDIUM = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Middle layer thickness (thin film thickness)
	json_object_object_get_ex(jsonControls, "MIDDLE_LAYER", &dummyValue);
	pids->MIDDLE_LAYER = (dummyValue == NULL) ? 0 : json_object_get_double(dummyValue);

	// Material index of final layer (exit layer)
	json_object_object_get_ex(jsonControls, "Last Medium", &dummyValue);
	pids->SECOND_MEDIUM = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Material index of nanostructure
	json_object_object_get_ex(jsonControls, "Antenna Medium", &dummyValue);
	pids->ANTENNAMATERIAL = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Flag to trigger nanostructure being drawn
	json_object_object_get_ex(jsonControls, "NANOANTENNA", &dummyValue);
	pids->NANOANTENNA = (dummyValue == NULL) ? 1 : json_object_get_int(dummyValue);

	// Nanostructure type
	json_object_object_get_ex(jsonControls, "ANTENNA_TYPE", &dummyValue);
	pids->ANTENNA_TYPE = (dummyValue == NULL) ? 1000 : json_object_get_int(dummyValue);; // read geometry from JSON file
	printf("Antenna type: %d\n", pids->ANTENNA_TYPE);

	// Various parameters for the nanostructure
	json_object_object_get_ex(jsonControls, "ANTENNA_SHIFT", &dummyValue);
	pids->ANTENNA_SHIFT = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "WIDTH", &dummyValue);
	pids->WIDTH = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "LENGTH", &dummyValue);
	pids->LENGTH = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "THICKNESS", &dummyValue);
	pids->THICKNESS = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "LENGTH2", &dummyValue);
	pids->LENGTH2 = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "GAP", &dummyValue);
	pids->GAP = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "LENGTH1_L", &dummyValue);
	pids->LENGTH1_L = json_object_get_double(dummyValue);
	json_object_object_get_ex(jsonControls, "LENGTH2_B", &dummyValue);
	pids->LENGTH2_B = json_object_get_double(dummyValue);

	// -----------------------------------------------------------------------
	// 						POST-PROCESSING CONTROLS
	// -----------------------------------------------------------------------

	// // Check for plots key in JSON (work in progress)
	// json_object_object_get_ex(jsonControls, "Plots", &dummyArray);
	// if (dummyArray != NULL)
	// {
	// 	// Check each one accordingly in a loop
	// 	int nPlotFlags = json_object_array_length(dummyArray);
	// 	for	(int plotFlagNumber = 0; plotFlagNumber < nPlotFlags; plotFlagNumber++)
	// 	{
	// 		// Read in current plot string
	// 		dummyValue = json_object_array_get_idx(dummyArray, plotFlagNumber);

	// 		// Determine which flag it sets
	// 		if ( strcmp(json_object_get_string(dummyValue), "E_DFT_3D") )
	// 		{
	// 			printf("Saving Frequency Domain Electric Fields in 3D\n");
	// 			pids->SAVE_OUTPUT = 1;
	// 		}
	// 	}
	// }

	// Which DFT cross sections to plot
	json_object_object_get_ex(jsonControls, "DFT Plot", &dummyValue);
	pids->DFT_PLOT = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Nonlinear simulation flag
	json_object_object_get_ex(jsonControls, "DFG_OPTIM", &dummyValue);
	pids->DFG_OPTIM = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// ABSOLUTELY NO IDEA
	json_object_object_get_ex(jsonControls, "DFG_THZ_FREQ", &dummyValue);
	pids->DFG_THZ_FREQ = json_object_get_double(dummyValue);

	// ABSOLUTELY NO IDEA
	json_object_object_get_ex(jsonControls, "DFG_DEEP_STEP", &dummyValue);
	pids->DFG_DEEP_STEP = json_object_get_double(dummyValue);

	// Flag for Graphene layer
	json_object_object_get_ex(jsonControls, "GRAPHENE_OPTIM", &dummyValue);
	pids->GRAPHENE_OPTIM = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "DFT Phase Normalization", &dummyValue);
	pids->DFT_PHASE_NORM = (dummyValue == NULL) ? 1 : json_object_get_int(dummyValue);

	// Time domain frames in YZ plane
	json_object_object_get_ex(jsonControls, "Movie YZ", &dummyValue);
	pids->PLOT_YZ = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Time domain frames in XZ plane
	json_object_object_get_ex(jsonControls, "Movie XZ", &dummyValue);
	pids->PLOT_XZ = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Time domain frames in XY plane
	json_object_object_get_ex(jsonControls, "Movie XY", &dummyValue);
	pids->PLOT_XY = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Time domain frames in XY plane
	json_object_object_get_ex(jsonControls, "Movie 3D", &dummyValue);
	pids->TIME_PLOT_3D = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	// Frequency domain output of entire domain E field
	json_object_object_get_ex(jsonControls, "Full Domain DFT", &dummyValue);
	pids->SAVE_OUTPUT = json_object_get_int(dummyValue);

	// absorption calculated by the sigma|E|^2 volume integration in ANTENNA Material
	json_object_object_get_ex(jsonControls, "Sigma Absorption", &dummyValue);
	pids->ABS_SIGMA = json_object_get_int(dummyValue);

	// Per-component or Staircasing
	json_object_object_get_ex(jsonControls, "Staircasing Type", &dummyValue);
	pids->STAIRCASING = (dummyValue == NULL) ? 0 : json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "Height", &dummyValue);
	pids->HEIGHT = json_object_get_double(dummyValue);

	json_object_object_get_ex(jsonControls, "Bloch Boundaries", &dummyValue);
	PBCCTW = json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "Advanced Model", &dummyValue);
	NonLoc = json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "Diffusion Coefficient", &dummyValue);
	pids->Diff_C = json_object_get_double(dummyValue);

	
	json_object_object_get_ex(jsonControls, "Pulse Width", &dummyValue);
	PULSEWIDTH = json_object_get_double(dummyValue);

	int adjoint;
	json_object_object_get_ex(jsonControls, "Adjoint", &dummyValue);
	adjoint = json_object_get_int(dummyValue);

	json_object_object_get_ex(jsonControls, "Multipoles", &dummyArray);
	if (dummyArray == NULL)
	{
		pids->MULTIPOLES = false;
	}
	else
	{
		pids->MULTIPOLES = true;
		SphereRad_MULTIPOLE = json_object_get_int(json_object_array_get_idx(dummyArray, 0));
		Ntheta_MULTIPOLE = json_object_get_int(json_object_array_get_idx(dummyArray, 1));
		Nphi_MULTIPOLE = json_object_get_int(json_object_array_get_idx(dummyArray, 2));
	}

	if (adjoint == 1)
		ADJOINT = true;
	if (adjoint == 0)
		ADJOINT = false;

	// These aren't nessecary
	pids->Fermi_V = 1.0;
	TIME_DOMAIN_PROBE = 0;

	pids->E_PEAK = E_FIELD_MAX;
	printf("E field max = %e\n", E_FIELD_MAX);

	pids->MODE_SOURCE = false;
	if (pids->SIGNAL_TYPE == 9)
		pids->MODE_SOURCE = true;

	pids->TwoTemperatureModel = false;
	SHG_HYDRO = false;
	pids->NONLOCAL = false;
	pids->Hydrodynamics = false;
	pids->GNOR = false;
	CALCULATE_POLARIZABILITY = false;
	// CALCULATE_NONLIN_SIGMA = false;
	CALCULATE_NONLIN_SIGMA = false;

	AddFreeElectronLayer = false;
	Seperate_Linear_Nonlinear = false;
	pids->SAVE_OUTPUT_TEMPERATURE = false;
	pids->Pump_Probe = false;
	pids->Pump = false;
	pids->Probe = false;

	// if(pids->ANTENNAMATERIAL == 4 || pids->ANTENNAMATERIAL == 5 || NonLoc > 0){
	// 			CALCULATE_POLARIZABILITY = true;
	// 		}
	THREE_D_FT = false;
	NoInterbandTransitions = false;
	multipoles = false;

	pids->E_PEAK = 1.0;
	pids->E_PEAK = E_FIELD_MAX;
	printf("E field max = %e\n", E_FIELD_MAX);

	SHELLMATERIAL = 2;

	// printf("%d\n",pids->NRCELLSCPMLY);
	printf("NonLoc: %d\n", NonLoc);
	if (NonLoc >= 1)
	{
		pids->NONLOCAL = true;
		pids->GNOR = false;
		pids->Hydrodynamics = false;
		AddFreeElectronLayer = false;
		EDensFT = true;
		if (NonLoc == 2)
		{
			pids->Hydrodynamics = false;
			pids->GNOR = true;
			printf("GNOR\n");
		}
		else if (NonLoc > 2)
		{
			pids->GNOR = false;
			pids->Hydrodynamics = true;
			// UpdateN = true;
			WithCritPoint = true;
			// printf("Hydrodynamics\n");
			pids->E_PEAK = E_FIELD_MAX;
			AddFreeElectronLayer = false;

			if (NonLoc == 3)
			{
				EDensFT = true;
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				pids->E_PEAK = E_FIELD_MAX;
				// printf("With Pressure\n");
				if (pids->ANTENNAMATERIAL < 17 && pids->SECOND_MEDIUM < 17 && pids->MIDDLE_MEDIUM < 17)
				{
					printf("CHOOSE HYDRODYNAMIC MATERIAL!\n");
					exit(-1);
				}
			}
			else if (NonLoc == 4)
			{
				pids->WithPressure = false;
				pids->WithConvection = true;
				pids->WithMagField = true;
				printf("Without Pressure\n");
			}
			else if (NonLoc == 5)
			{
				pids->WithPressure = false;
				pids->WithConvection = false;
				pids->WithMagField = true;
			}
			else if (NonLoc == 6)
			{
				pids->WithPressure = false;
				pids->WithConvection = false;
				pids->WithMagField = false;
			}
			else if (NonLoc == 7)
			{
				pids->WithPressure = false;
				pids->WithConvection = true;
				pids->WithMagField = false;
			}
			else if (NonLoc == 8)
			{
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = false;
			}
			else if (NonLoc == 9)
			{
				pids->WithPressure = true;
				pids->WithConvection = false;
				pids->WithMagField = true;
			}
			else if (NonLoc == 10)
			{
				pids->WithPressure = false;
				pids->WithConvection = false;
				pids->WithMagField = false;
				UpdateN = false;
				pids->E_PEAK = E_FIELD_MAX;
			}
			else if (NonLoc == 11)
			{
				pids->WithPressure = true;
				pids->WithConvection = false;
				pids->WithMagField = false;
				UpdateN = false;
			}
			else if (NonLoc == 12)
			{
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				UpdateN = false;
			}
			else if (NonLoc == 13)
			{
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				UpdateN = true;
				SHG_HYDRO = true;
				pids->E_PEAK = E_FIELD_MAX;
			}
			else if (NonLoc == 14)
			{
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				UpdateN = true;
				WithCritPoint = false;
			}
			else if (NonLoc == 15)
			{
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				UpdateN = false;
				WithCritPoint = false;
			}

			else if (NonLoc == 200)
			{
				multipoles = true;
				pids->WithPressure = false;
				pids->WithConvection = true;
				pids->WithMagField = true;
				pids->Seperate_Linear_Nonlinear = true;
				// CALCULATE_NONLIN_SIGMA = true;
				EDensFT = true;

				printf("Seperating Linear and Nonlinear \n");
			}

			else if (NonLoc == 201)
			{
				multipoles = true;
				pids->WithPressure = true;
				pids->WithConvection = true;
				pids->WithMagField = true;
				pids->Seperate_Linear_Nonlinear = true;
				printf("Seperating Linear and Nonlinear \n");
				EDensFT = true;

				// AddFreeElectronLayer = true;
				// SHELLMATERIAL = pids->ANTENNAMATERIAL + NumNonLinFreeElectron;
				// printf("SHELL WITH MATERIAL: %d\t",SHELLMATERIAL);
			}
			else if (NonLoc == 202)
			{
				pids->NONLOCAL = false;
				pids->TwoTemperatureModel = true;
				pids->Hydrodynamics = false;
				EDensFT = false;
				pids->LOAD_TEMPERATURE = false;
			}
			else if (NonLoc == 2020)
			{
				pids->NONLOCAL = false;
				pids->Seperate_Linear_Nonlinear = false;

				pids->TwoTemperatureModel = true;
				pids->Hydrodynamics = false;
				EDensFT = false;
				pids->LOAD_TEMPERATURE = false;
				if (pids->PLOT_YZ > 0 || pids->PLOT_XZ > 0 || pids->PLOT_XY > 0)
					pids->SAVE_OUTPUT_TEMPERATURE = true;
			}
			else if (NonLoc == 2021)
			{
				pids->NONLOCAL = false;
				pids->TwoTemperatureModel = false;
				pids->Hydrodynamics = false;
				EDensFT = false;
				pids->LOAD_TEMPERATURE = true;
			}

			else if (NonLoc == 2022)
			{
				pids->NONLOCAL = false;
				pids->Seperate_Linear_Nonlinear = false;
				pids->TwoTemperatureModel = true;
				pids->Hydrodynamics = false;
				EDensFT = false;
				pids->LOAD_TEMPERATURE = false;
				pids->SAVE_OUTPUT_TEMPERATURE = false;
				pids->Pump_Probe = true;
				pids->Pump = true;
				pids->Probe = false;
			}
			else if (NonLoc == 2023)
			{
				pids->NONLOCAL = false;
				pids->Seperate_Linear_Nonlinear = false;
				pids->TwoTemperatureModel = true;
				pids->Hydrodynamics = false;
				EDensFT = false;
				pids->LOAD_TEMPERATURE = false;
				pids->SAVE_OUTPUT_TEMPERATURE = false;
				pids->Pump_Probe = true;
				pids->Pump = false;
				pids->Probe = true;
			}
			else if (NonLoc == 203)
			{
				pids->NONLOCAL = true;
				pids->TwoTemperatureModel = true;
				pids->Hydrodynamics = true;
				pids->Seperate_Linear_Nonlinear = false;
				pids->LOAD_TEMPERATURE = false;
				if (pids->PLOT_YZ > 0 || pids->PLOT_XZ > 0 || pids->PLOT_XY > 0)
					pids->SAVE_OUTPUT_TEMPERATURE = true;
				EDensFT = false;
			}
		}

		else
		{
			pids->GNOR = false;
			pids->Hydrodynamics = false;
		}
		printf("NONLOCAL\n");
	}
	else
	{
		pids->NONLOCAL = false;
		pids->Hydrodynamics = false;
	}
	GRADIENT_DIVERGENCE = true;
	LAPLACIAN = false;
	if (PBCCTW == 1 || PBCCTW == 2)
	{
		PBC_CTW = true;

#ifdef DOUBLEPRECISION // Error Check
		FailureType = NO_COMPLEX_FIELDS;
		return true;
#endif
		if (pids->SIGNAL_TYPE != 6)
		{
			// printf("Switching to Modulated Gaussian for Constant Transverse Wavenumber\n");
			// pids->SIGNAL_TYPE = 6;
		}
	}

	if (TwoTemperatureModel)
	{
		string strS = BASE_PATH_INPUT + "/NumRK.txt";
		fstream pfin;
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			pids->NumRK = 1;
		}
		else
		{
			pids->NumRK = REAL(atoi(ReadFromTxtFile(pfin).c_str()));
		}

		pfin.close();
	}
	int TM_or_TE;
	if (PBC_CTW)
	{
		string strS = BASE_PATH_INPUT + "/pphinfoCTW.txt";
		fstream pfin;
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			CHIRP = 0.0;
			CHIRP2 = 0.0;

			CTWOffset = 0;
			trials = int(atoi(ReadFromTxtFile(pfin).c_str()));
			NUM_trials = int(atof(ReadFromTxtFile(pfin).c_str()));
			TM_or_TE = int(atof(ReadFromTxtFile(pfin).c_str()));
			k_max = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->STARTLAMBDA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->STOPLAMBDA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			CENTERFREQ = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			PULSEWIDTH = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			CHIRP = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			CHIRP2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));

			printf("trials = %d\n", trials);
			printf("NUM_trials = %d\n", NUM_trials);
			printf("Tend_inc = %d\n", Tend_inc);
		}

		pfin.close();

		if (TM_or_TE == 0)
		{
			TEy = true;
			TMy = false;
		}
		else if (TM_or_TE == 1)
		{
			TMy = true;
			TEy = false;
		}
		else
			TEy = TMy = true;
	}

	if (pids->Pump_Probe && pids->Probe)
	{
		string strS = BASE_PATH_INPUT + "/PumpProbeDelay.txt";
		fstream pfin;
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->PumpProbeDelay = REAL(atof(ReadFromTxtFile(pfin).c_str()));
		}
		pfin.close();
	}

	int probes;
	probe_structure->Number_of_Probes = 0;
	strS = strPath + "/Probes.txt";
	pfin.open(strS.c_str(), ios::in);
	printf("%s\t%s\n", strS.c_str(), "LOADED");
	if (!pfin.is_open())
	{
		printf("No Probes\n");
	}
	else
	{
		probe_structure->Number_of_Probes = atoi(ReadFromTxtFile(pfin).c_str());
		for (probes = 0; probes < probe_structure->Number_of_Probes; probes++)
		{
			probe_structure->Probe_X[probes] = atoi(ReadFromTxtFile(pfin).c_str());
			probe_structure->Probe_Y[probes] = atoi(ReadFromTxtFile(pfin).c_str());
			probe_structure->Probe_Z[probes] = atoi(ReadFromTxtFile(pfin).c_str());
			printf("Probes: %d\t%d\t%d\t%d\n", probe_structure->Number_of_Probes, probe_structure->Probe_X[probes], probe_structure->Probe_Y[probes], probe_structure->Probe_Z[probes]);
		}

		pfin.close();
	}

	if (pids->ANTENNAMATERIAL == NonLinFreeElectron + DispMateOffsetD + 1 || pids->FIRST_MEDIUM == NonLinFreeElectron + DispMateOffsetD + 1 || pids->MIDDLE_MEDIUM == NonLinFreeElectron + DispMateOffsetD + 1 || pids->SECOND_MEDIUM == NonLinFreeElectron + DispMateOffsetD + 1)
	{
		string strS = BASE_PATH_INPUT + "/pphREADINMATERIAL.txt";
		fstream pfin;
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->READ_IN_EPREINF = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_OME_PLASMA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_GAM_PLASMA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGOME_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGGAM_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_ACOEFF_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
		}

		pfin.close();
	}

	if (pids->ANTENNAMATERIAL == ModLorentz + DispMateOffsetD + 2 || pids->FIRST_MEDIUM == ModLorentz + DispMateOffsetD + 2 || pids->MIDDLE_MEDIUM == ModLorentz + DispMateOffsetD + 2 || pids->SECOND_MEDIUM == ModLorentz + DispMateOffsetD + 2)
	{
		string strS = BASE_PATH_INPUT + "/pphREADINMATERIAL.txt";
		fstream pfin;
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->READ_IN_EPREINF = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_OME_PLASMA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_GAM_PLASMA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGOME_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGOME_CRIPOI_2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGGAM_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_BIGGAM_CRIPOI_2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_ACOEFF_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_ACOEFF_CRIPOI_2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_PHI_CRIPOI = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_PHI_CRIPOI_2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_FERMI_VELOCITY = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->READ_IN_ELECTRON_DENSITY = REAL(atof(ReadFromTxtFile(pfin).c_str()));
		}

		pfin.close();
	}

	// To enable nonlinear simulations (polarization approach - Antonio, 2017 papers)
	/* if(pids->DFG_OPTIM == 3 || pids->DFG_OPTIM == 4 || pids->DFG_OPTIM == 5 || pids->DFG_OPTIM == 6 || (pids->ANTENNA_TYPE==58 && pids->GAP>0) )
		{
	strS = strPath + "/pphshgini.txt";
	pfin.open( strS.c_str(), ios::in );
	printf("%s\t%s\n",strS.c_str(),"LOADED");
	if ( !pfin.is_open() ) {
		FailureType = NO_READABLE_INIT_DATA;
		return true;

	}
	else {
		pids->NL_mat_gap= atoi( ReadFromTxtFile(pfin).c_str() );
		pids->CHI_NL_ANTENNA = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		pids->CHI_NL_GAP = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		pids->CHI_NL_FIRST = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		pids->CHI_NL_MIDDLE = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		pids->CHI_NL_SECOND = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		pids->E_PEAK = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );
		}

	pfin.close();
	} */

	//-----------------------------Izzat-----------------------------------------------
	if (pids->DFG_OPTIM == 5 || pids->DFG_OPTIM == 6) // SHG
	{
		strS = strPath + "/pphnlsini.txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{ // Is the external file readable/exist?
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{ // If the file exist then read...
			// string NL_ANTENNA_MAT, NL_FIRST_MAT, NL_MIDDLE_MAT, NL_SECOND_MAT;
			NL_ANTENNA_MAT = ReadFromTxtFile(pfin).c_str();					// antenna material (string)
			NL_FIRST_MAT = ReadFromTxtFile(pfin).c_str();					// first layer (superstrate) material (string)
			NL_MIDDLE_MAT = ReadFromTxtFile(pfin).c_str();					// middle layer (superstrate) material (string)
			NL_SECOND_MAT = ReadFromTxtFile(pfin).c_str();					// second layer (substrate) material (string)
			pids->PNL_APPROACH = REAL(atoi(ReadFromTxtFile(pfin).c_str())); // Nonlinear approach:  0=Electric field update, 1=Polarization update
		}
		pfin.close();

		// open chi2=2*dmn tensor elements of antenna material
		strS = strPath + "/chi2/" + NL_ANTENNA_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi2=2*dmn tensor elements of antenna material
			for (i = 0; i < 18; ++i)
			{
				pids->CHI2A[i] = REAL(2 * atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi2=2*dmn tensor elements of first layer (superstrate) material
		strS = strPath + "/chi2/" + NL_FIRST_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi2=2*dmn tensor elements of first layer (superstrate) material
			for (i = 0; i < 18; ++i)
			{
				pids->CHI2F[i] = REAL(2 * atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi2=2*dmn tensor elements of middle layer material
		strS = strPath + "/chi2/" + NL_MIDDLE_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi2=2*dmn tensor elements of middle layer material
			for (i = 0; i < 18; ++i)
			{
				pids->CHI2M[i] = REAL(2 * atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi2=2*dmn tensor elements of second layer (substrate) material
		strS = strPath + "/chi2/" + NL_SECOND_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi2=2*dmn tensor elements of second layer (substrate) material
			for (i = 0; i < 18; ++i)
			{
				pids->CHI2S[i] = REAL(2 * atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();
	}

	if (pids->DFG_OPTIM == 3 || pids->DFG_OPTIM == 4) // THG
	{
		// Read and open external THG initialization file "/pphshgini.txt"
		strS = strPath + "/pphnlsini.txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{ // Is the external file readable/exist?
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{																	// If the file exist then read...
			NL_ANTENNA_MAT = ReadFromTxtFile(pfin).c_str();					// antenna material (string)
			NL_FIRST_MAT = ReadFromTxtFile(pfin).c_str();					// first layer (superstrate) material (string)
			NL_MIDDLE_MAT = ReadFromTxtFile(pfin).c_str();					// middle layer (superstrate) material (string)
			NL_SECOND_MAT = ReadFromTxtFile(pfin).c_str();					// second layer (substrate) material (string)
			pids->PNL_APPROACH = REAL(atoi(ReadFromTxtFile(pfin).c_str())); // Nonlinear approach:  0=Electric field update, 1=Polarization update
		}
		pfin.close();

		// open chi3 tensor of antenna material
		strS = strPath + "/chi3/" + NL_ANTENNA_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi3 tensor elements of antenna material
			for (i = 0; i < 21; ++i)
			{
				pids->CHI3A[i] = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi3 tensor of first layer (superstrate) material
		strS = strPath + "/chi3/" + NL_FIRST_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi3 tensor elements of first layer (superstrate) material
			for (i = 0; i < 21; ++i)
			{
				pids->CHI3F[i] = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi3 tensor of middle layer material
		strS = strPath + "/chi3/" + NL_MIDDLE_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi3 tensor elements of middle layer material
			for (i = 0; i < 21; ++i)
			{
				pids->CHI3M[i] = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();

		// open chi3 tensor of second layer (substrate) material
		strS = strPath + "/chi3/" + NL_SECOND_MAT.c_str() + ".txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			// chi3 tensor elements of second layer (substrate) material
			for (i = 0; i < 21; ++i)
			{
				pids->CHI3S[i] = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			}
		}
		pfin.close();
	}
	//-------------------------------------------------------------------------------------

	if ((pids->ANTENNA_TYPE == 1000) || (pids->ANTENNA_TYPE == 81) || (pids->ANTENNA_TYPE == 82) || (pids->ANTENNA_TYPE == 83) || (pids->ANTENNA_TYPE == 84) || (pids->ANTENNA_TYPE == 85) || (pids->ANTENNA_TYPE == 75) || (pids->ANTENNA_TYPE == 49) || (pids->ANTENNA_TYPE == 1001) || (pids->ANTENNA_TYPE == 24) || (pids->ANTENNA_TYPE == 25) || (pids->ANTENNA_TYPE == 26) || (pids->ANTENNAMATERIAL > int(0) && pids->ANTENNAMATERIAL < DispMateOffsetD) || (pids->FIRST_MEDIUM > int(0) && pids->FIRST_MEDIUM < DispMateOffsetD) || (pids->MIDDLE_MEDIUM > int(0) && pids->MIDDLE_MEDIUM < DispMateOffsetD) || (pids->SECOND_MEDIUM > int(0) && pids->SECOND_MEDIUM < DispMateOffsetD) || (pids->NL_mat_gap > int(0) && pids->NL_mat_gap < DispMateOffsetD) || (pids->ANTENNA_TYPE == 117))
	{
		strS = strPath + "/pphmatini.txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->NREFR1 = REAL(atof(ReadFromTxtFile(pfin).c_str())); // NREFR1
			pids->NREFR2 = REAL(atof(ReadFromTxtFile(pfin).c_str())); // NREFR2
			pids->NREFR3 = REAL(atof(ReadFromTxtFile(pfin).c_str())); // NREFR3
			for (i = 4; i < 100; i++)
			{
				NREFR[i] = REAL(atof(ReadFromTxtFile(pfin).c_str()));
				printf("%e\n", NREFR[i]);
			}

			if (pids->ANTENNA_TYPE == 81 || pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 83 || pids->ANTENNA_TYPE == 84 || pids->ANTENNA_TYPE == 85)
				pids->CARRIER_DENSITY = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			if (pids->ANTENNA_TYPE == 81 || pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 83 || pids->ANTENNA_TYPE == 84 || pids->ANTENNA_TYPE == 85)
				pids->CARRIER_DENSITY_BACK = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			if (pids->ANTENNA_TYPE == 81 || pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 83 || pids->ANTENNA_TYPE == 84 || pids->ANTENNA_TYPE == 85)
				pids->GAMMA_FACT_ANTENNA = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			if (pids->ANTENNA_TYPE == 81 || pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 83 || pids->ANTENNA_TYPE == 84 || pids->ANTENNA_TYPE == 85)
				pids->GAMMA_FACT_PERT_ITO = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			if (pids->ANTENNA_TYPE == 81 || pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 83 || pids->ANTENNA_TYPE == 84 || pids->ANTENNA_TYPE == 85)
				pids->GAMMA_FACT_NONPERT_ITO = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			if (pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 84)
				pids->PERC_EXTRA_LAYER_LIDAR = REAL(atof(ReadFromTxtFile(pfin).c_str())); // percentage of MIDDLELAYER
			if (pids->ANTENNA_TYPE == 82 || pids->ANTENNA_TYPE == 84)
				pids->MAT_EXTRA_LAYER_LIDAR = atoi(ReadFromTxtFile(pfin).c_str());
			if (pids->ANTENNA_TYPE == 49)
				pids->NREFR0 = REAL(atof(ReadFromTxtFile(pfin).c_str())); // NREFR0
			if (pids->ANTENNA_TYPE == 24 || pids->ANTENNA_TYPE == 25)
				pids->STRETCH_EPS_X = REAL(atof(ReadFromTxtFile(pfin).c_str())); // flakes
			if (pids->ANTENNA_TYPE == 24 || pids->ANTENNA_TYPE == 25)
				pids->STRETCH_EPS_Z = REAL(atof(ReadFromTxtFile(pfin).c_str())); // flakes
			if (pids->ANTENNA_TYPE == 24 || pids->ANTENNA_TYPE == 26)
				pids->SHEETS_DISTANCE = REAL(atof(ReadFromTxtFile(pfin).c_str())); // flakes
																				   // if( pids->ANTENNA_TYPE==84 )pids->LENGTH2_UP = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );	//lidar_double_emb
																				   // if( pids->ANTENNA_TYPE==84 )pids->LENGTH2_DOWN = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );	//lidar_double_emb
																				   // if( pids->ANTENNA_TYPE==84 )pids->ANTENNA_DIST = REAL( atof ( ReadFromTxtFile(pfin).c_str() ) );	//lidar_double_emb
		}
		pfin.close();
	}
	if ((pids->ANTENNA_TYPE == 82) || (pids->ANTENNA_TYPE == 84))
	{
		strS = strPath + "/pphlidar.txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\t%s\n", strS.c_str(), "LOADED");
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->CONNECTOR_POS = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->CONNECTOR_WIDTH = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->GAP2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->LENGTH2_UP = REAL(atof(ReadFromTxtFile(pfin).c_str()));	// lidar_double_emb
			pids->LENGTH2_DOWN = REAL(atof(ReadFromTxtFile(pfin).c_str())); // lidar_double_emb
			pids->ANTENNA_DIST = REAL(atof(ReadFromTxtFile(pfin).c_str())); // lidar_double_emb
			pids->WIDTH2 = REAL(atof(ReadFromTxtFile(pfin).c_str()));		// lidar_double_emb
		}
		pfin.close();
	}

	if (pids->TFSF == 3)
	{
		strS = strPath + "/gaussian_beam.txt";
		pfin.open(strS.c_str(), ios::in);
		printf("%s\n", strS.c_str());
		if (!pfin.is_open())
		{
			FailureType = NO_READABLE_INIT_DATA;
			return true;
		}
		else
		{
			pids->THETA_NR_INTERVALS = REAL(atoi(ReadFromTxtFile(pfin).c_str()));
			pids->PHI_NR_INTERVALS = REAL(atoi(ReadFromTxtFile(pfin).c_str()));
			pids->THETA_MAX_DEGREES = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->PHI_MAX_DEGREES = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->BEAM_W0_COEFF = REAL(atof(ReadFromTxtFile(pfin).c_str()));
			pids->FOCAL_LENGTH_COEFF = REAL(atof(ReadFromTxtFile(pfin).c_str()));
		}

		pfin.close();
	}

	return false;
}


/**
 * @brief Read in custom DFT JSON file.
 *
 * @details Custom DFT monitors can be set (in bulk) using this function.
 * Assumes the path "./INIDEF/dfts.json", and unpacks the information NOTE: fuction has been modified to not asssume path anymore but use global input path 
 * into a global variables "dft_struct".
 *
 */
void ReadDFTPositions()
{
	fstream pfin;
	string strS;
	FILE *fp;
	char buffer[102422];
	size_t length_file;
	struct json_object *allStructures;
	struct json_object *dummyValue;
	struct json_object *jsonControls;
	struct json_object *dummyArray;

	int Size_X, Size_Y, Size_Z;
	int Center_X, Center_Y, Center_Z;

	string controlPathJSON = BASE_PATH_INPUT + "/dfts.json";
	fp = fopen(controlPathJSON.c_str(), "r");
	if (fp == NULL)
	{
		return;
	}
	else
	{
		size_t successReadCount = fread(buffer, 102422, 1, fp);
		fclose(fp);

		allStructures = json_tokener_parse(buffer);
		dft_struct->Num_DFT = json_object_array_length(allStructures);
		printf("DFTs: %d\n", dft_struct->Num_DFT);
		for (int p = 0; p < dft_struct->Num_DFT; p++)
		{
			jsonControls = json_object_array_get_idx(allStructures, p);
			json_object_object_get_ex(jsonControls, "Size", &dummyArray);
			dummyValue = json_object_array_get_idx(dummyArray, 0);
			Size_X = json_object_get_int(dummyValue);
			dummyValue = json_object_array_get_idx(dummyArray, 1);
			Size_Y = json_object_get_int(dummyValue);
			dummyValue = json_object_array_get_idx(dummyArray, 2);
			Size_Z = json_object_get_int(dummyValue);

			jsonControls = json_object_array_get_idx(allStructures, p);
			json_object_object_get_ex(jsonControls, "Center", &dummyArray);
			dummyValue = json_object_array_get_idx(dummyArray, 0);
			Center_X = json_object_get_int(dummyValue);
			dummyValue = json_object_array_get_idx(dummyArray, 1);
			Center_Y = json_object_get_int(dummyValue);
			dummyValue = json_object_array_get_idx(dummyArray, 2);
			Center_Z = json_object_get_int(dummyValue);

			dft_struct->min_X[p] = int(ceil(Center_X - 0.5 * Size_X));
			dft_struct->max_X[p] = int(ceil(Center_X + 0.5 * Size_X));
			dft_struct->min_Y[p] = int(ceil(Center_Y - 0.5 * Size_Y));
			dft_struct->max_Y[p] = int(ceil(Center_Y + 0.5 * Size_Y));
			dft_struct->min_Z[p] = int(ceil(Center_Z - 0.5 * Size_Z));
			dft_struct->max_Z[p] = int(ceil(Center_Z + 0.5 * Size_Z));
			printf("DFT %d, X: %d,%d, Y: %d,%d, Z: %d,%d \n", p, dft_struct->min_X[p], dft_struct->max_X[p], dft_struct->min_Y[p], dft_struct->max_Y[p], dft_struct->min_Z[p], dft_struct->max_Z[p]);
			printf("DFT %d, Centered: %d,%d,%d, Size: %d,%d,%d\n", p, Center_X, Center_Y, Center_Z, Size_X, Size_Y, Size_Z);
		}
	}
	return;
}
