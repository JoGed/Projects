#ifndef TOPOPT_H
#define TOPOPT_H

#include <petsc.h>
// #include <petsc-private/dmdaimpl.h>
#include "MMA.h"
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <petsc/private/dmdaimpl.h>
#include <sstream>
using namespace std;

/*
 Adapted Code for electromagnetic optimization problems, modified by Johannes Gedeon
 Original Code:
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
 Updated: June 2019, Niels Aage
Copyright (C) 2013-2019,

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

/*
 *
 * Parameter container for the topology optimization problem
 *
 * min_x fx
 * s.t. gx_j <= 0, j=1..m
 *      xmin_i <= x_i <= xmax_i, i=1..n
 *
 * with filtering and a volume constraint
 *
 */

class TopOpt {

  public:

    // Max. number of regions that will be stored
    static const PetscInt NUM_MAX_REGIONS = 10;
    // Element is 'true', if corresponding design region is defined in Json file
    static bool DRisDefined[NUM_MAX_REGIONS];
    // Element is 'true', if corresponding observation region is defined in Json file
    static bool ORisDefined[NUM_MAX_REGIONS];
    // ((xmin, xmax), (ymin, ymax), (zmin, zmax)) list per max. number of opt. regions
    static PetscInt boundsDR[NUM_MAX_REGIONS][3][2]; // Design regions boundary indices
    static PetscInt boundsOR[NUM_MAX_REGIONS][3][2]; // Objective regions boundary indices

    // Constructor/Destructor
    TopOpt();
    ~TopOpt();

    // Reads TopOpt parameters from .json file
    void ReadParamsFromFile(std::string path, int rankRead);
    // Setup DMDA and initialize vectors
    PetscErrorCode SetUp(PetscInt identifier, PetscInt numProcDR[3], PetscInt* cellsArrDR[3], PetscReal spacing[3]);   
    // Process of rankSend sends params to processes of comm. "comm"
    void SendParams(MPI_Comm comm, PetscInt rankSend);
    // Save densities into binary file
    PetscErrorCode SaveDensities(std::string prefix, bool* saveVecs);
    // Initialize densities to a sphere (for testing)
    PetscErrorCode InitializeDensityToSphere(PetscInt centerx, PetscInt centery, PetscInt centerz, 
                                             PetscReal SphereRadius, bool staircasing);

    // Method to allocate MMA with/without restarting
    PetscErrorCode AllocateMMAwithRestart(PetscInt* itr, MMA** mma);
    PetscErrorCode WriteRestartFiles(PetscInt* itr, MMA* mma);

    PetscInt     identifierDR; // Index of design Region in the simulation domain
    // Yee grid (basis for design)
    DM da_yee_lattice;
    // Optimization parameters
    PetscInt     n;            // Total number of yee cells
    PetscInt     nloc;         // Local number of yee cells?
    PetscInt     m;            // Number of constraints
    PetscScalar  fx;           // Objective value
    PetscScalar  fscale;       // Scaling factor for objective
    PetscScalar* gx;           // Array with constraint values
    PetscScalar  Densmin;      // Min. value of design variables
    PetscScalar  Densmax;      // Max. value of design variables
    PetscScalar  movlim;       // Max. change of design variables
    PetscScalar  volfrac;      // Volume fraction
    PetscInt     maxItr;       // Max iterations  ( == max. number of full simulations, NOT FDTD iterations!)
    PetscScalar  rmin;         // Filter radius (for uniform filtering between density components) in unit length
    PetscScalar  rmin_M[3][3]; // Filter radius matrix (r_xx, r_xy, r_xz, ...), generalization of rmin in unit length!
    PetscScalar  rminVec[6];   // Only needed to build radius matrix
    PetscScalar  wmin;         // Weight (for uniform filtering between density components)
    PetscScalar  wminVec[6];   // Only needed to build weight matrix
    PetscScalar  w_M[3][3];    // Weight matrix (w_xx, w_xy, w_xz, ...)
    PetscInt     filter;       // Filter type: 0=sens, 1=dens, val == no filtering
    PetscReal    beta;
    PetscReal    betaFinal;
    PetscReal    eta;
    PetscBool    projectionFilter; // Smooth heaviside projectionFilter
    PetscReal    unitLength;     // Unit lenght used for coordiantes (including dx, dy, dz) in meters
    Vec  dfdx;             // Sensitivities of objective
    Vec* dgdx;             // Sensitivities of constraints (vector array)
    Vec  densmin, densmax; // Vectors with max and min values of x
    
    Vec dens;      // Global vector containing density values of designable nanostructure on the Ex, Ey and Ez subgrid
    Vec dens_old;  // Global vector containing density values of designable nanostructure on the Ex, Ey and Ez subgrid. 
    Vec densTilde; // Global vector containing filtered density values of designable nanostructure on the Ex, Ey and Ez subgrid. 
    Vec densPhys;  // Global vector containing projected density values of designable nanostructure on the Ex, Ey and Ez subgrid
    
    PetscInt cellsDR[3]; // Total number of cells per axis
    PetscInt startIdx[3];// Startindices (globally)
    PetscInt locLength[3]; // Local chunk size (not including ghos points)
    PetscInt numGPr[3]; // Local number of "right" ghost points
    PetscReal spacing[3]; // spacing dx, dy, dz in unitLength

    PetscInt petscRank; // Rank in PETSc's communicator
    PetscInt petscSize; // Size of PETSc communicator

    PetscReal decayCond; // Minimum remainder of fractional value of objective values measured in neighbored time windows  
    PetscInt decayWindow; // Lenght of each time window (in timesteps), where objective function is computed to measure decay condition
    //----------------------------------------------------
    // Restart data for MMA:
    PetscBool   restart, flip;
    std::string restdens_1, restdens_2;
    Vec         xo1, xo2, U, L;
    //----------------------------------------------------
    inline void InputError(std::string err, int rankPrint) {
      int rankThis;
      MPI_Comm_rank(MPI_COMM_WORLD, &rankThis);
      if (rankThis == rankPrint){
        std::string str = "";
        std::string strSep = "################################\n";
        str += strSep + "Input error: " + err + "\n" + strSep;
        fprintf(stderr, str.c_str());
      }
	    MPI_Abort(MPI_COMM_WORLD, 1);
    };
    
  private:
    // Allocate and set default values
    void Init();
    PetscErrorCode SetUpMESH(PetscInt numProcDR[3], PetscInt* cellsArrDR[3]);    
    PetscErrorCode SetUpOPT();    
    // ----------------------------
    // Memory Management functions
    // ----------------------------
    void Allocate3DArray(PetscScalar*** &arr);
    void Deallocate3DArray(PetscScalar*** arr);
    // ----------------------------
    
    // Restart filenames
    std::string filename00, filename00Itr, filename01, filename01Itr;

    // File existence
    inline PetscBool fexists(const std::string& filename) {
        std::ifstream ifile(filename.c_str());
        if (ifile) {
            return PETSC_TRUE;
        }
        return PETSC_FALSE;
    }
};

#endif


