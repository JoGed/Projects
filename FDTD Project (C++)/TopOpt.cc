#include "TopOpt.h"
#include <cmath>
#include <json-c/json.h>
#include "io.h"
#include "RegionEmbedding.h"

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


TopOpt::TopOpt() {

    Init();
}

void TopOpt::Init() { 
    
    m = 1; // default (number of constraints)

    dens      = NULL;
    dens_old  = NULL;
    densTilde = NULL;
    densPhys  = NULL;
    densmin = NULL;
    densmax = NULL;

    dfdx = NULL;
    dgdx = NULL;
    gx   = NULL;
    xo1  = NULL;
    xo2  = NULL;
    U    = NULL;
    L    = NULL;
    da_yee_lattice = NULL;

    //-------------------------------------------------
    // SET DEFAULTS for optimization problems
    volfrac = 0.5;
    maxItr  = 1;
    wmin    = 1;
    filter  = 1;   // 0=sens, 1=dens, val == no filtering
    Densmin = 0.0;
    Densmax = 1.0;
    movlim  = 0.2;
    restart = PETSC_FALSE; // default True!

    //-------------------------------------------------
    // Projection filter
    projectionFilter = PETSC_TRUE; // default False!
    beta             = 0.1;//0.1;
    betaFinal        = 48;
    eta              = 0.0;  
    //-------------------------------------------------
    cellsDR[3] = {};
    decayCond = 0;
    decayWindow = 0;
}

TopOpt::~TopOpt() {

    // Delete vectors
    if (dens != NULL) {
        VecDestroy(&dens);
    }
    if (dens_old != NULL) {
        VecDestroy(&dens_old);
    }
    if (densTilde != NULL) {
        VecDestroy(&densTilde);
    }
    if (densPhys != NULL) {
        VecDestroy(&densPhys);
    }
    if (dfdx != NULL) {
        VecDestroy(&dfdx);
    }
    if (dgdx != NULL) {
        VecDestroyVecs(m, &dgdx);
    }
    if (densmin != NULL) {
        VecDestroy(&densmin);
    }
    if (densmax != NULL) {
        VecDestroy(&densmax);
    }
    if (da_yee_lattice != NULL) {
        DMDestroy(&(da_yee_lattice));
    }
    //Delete constraints
    if (gx != NULL) {
        delete[] gx;
    }

    // mma restart method
    if (xo1 != NULL) {
        VecDestroy(&xo1);
    }
    if (xo2 != NULL) {
        VecDestroy(&xo2);
    }
    if (L != NULL) {
        VecDestroy(&L);
    }
    if (U != NULL) {
        VecDestroy(&U);
    }
}

void TopOpt::ReadParamsFromFile(std::string path, int rankRead) { 
    
    // Specification of Design Region, unit_length, Filter radius (OR Filter radius matrix) expected as a minimum
    int RMatrix[6]; // 6 radius entries fine, since matrix is symmetric
    int WMatrix[6]; // 6 weight entries fine, since matrix is symmetric
    FILE *fp;
    int buff_size = 100000;
    char buffer[buff_size]; 
    struct json_object *parsed_json;
    struct json_object *params;
    struct json_object *val, *valR, *valRMatrix, *valWMatrix;
    
    fp = fopen(path.c_str(), "r");
    
    if (fp == NULL){
        printf("------------------------------------\n");
        perror (path.c_str());
        printf("------------------------------------\n");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
    else{
        printf("%s\t%s\n", path.c_str(), "LOADED");
        fread(buffer, buff_size, 1, fp);
        fclose(fp);
        params = json_object_array_get_idx(json_tokener_parse(buffer), 0); 
        // -------------------------------------------
        json_object_object_get_ex(params, "Unit length [SI]", &val);
        if(val){
            unitLength = json_object_get_double(val);   
        }
        else{
            InputError("Specifiy a 'Unit length [SI]'!", rankRead);
        }    
        // -------------------------------------------
        for (int i = 0; i < NUM_MAX_REGIONS; i++){
            std::string strParam;; 
            // -------------------------------------------
            strParam = "DR " + to_string(i) + " in X [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsDR[i][0][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsDR[i][0][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }
            // -------------------------------------------
            strParam = "DR " + to_string(i) + " in Y [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsDR[i][1][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsDR[i][1][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }
            // -------------------------------------------
            strParam = "DR " + to_string(i) + " in Z [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsDR[i][2][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsDR[i][2][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }

            // -------------------------------------------
            strParam = "OR " + to_string(i) + " in X [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsOR[i][0][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsOR[i][0][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }
            // -------------------------------------------
            strParam = "OR " + to_string(i) + " in Y [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsOR[i][1][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsOR[i][1][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }
            // -------------------------------------------
            strParam = "OR " + to_string(i) + " in Z [Idx]";
            json_object_object_get_ex(params, strParam.c_str(), &val);
            if(val){
                boundsOR[i][2][0] = json_object_get_int(json_object_array_get_idx(val, 0));
                boundsOR[i][2][1] = json_object_get_int(json_object_array_get_idx(val, 1));
            }
            else if (!val && i == 0){
                InputError("'" + strParam + "' has to be specified!", rankRead);
            }     
        }
        // -------------------------------------------
        json_object_object_get_ex(params, "Minimum Objective Change", &val);
        if(val) decayCond = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Window of Change Evaluation [timesteps]", &val);
        if(val) decayWindow = json_object_get_int(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Initial density value", &val);
        if(val) volfrac = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Number of Iterations", &val);
        if(val) maxItr = json_object_get_int(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Filter radius (default) [SI]", &valR);
        if(valR) rmin = json_object_get_double(valR);
        // -------------------------------------------
        json_object_object_get_ex(params, "Filtertype", &val);
        if(val) filter = json_object_get_int(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Density min", &val);
        if(val) Densmin = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Density max", &val);
        if(val) Densmax = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Moving Limit", &val);
        if(val) movlim = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Restart", &val);
        if(val) restart = json_object_get_boolean(val) ? PETSC_TRUE : PETSC_FALSE;
        // -------------------------------------------
        json_object_object_get_ex(params, "Projection filter", &val);
        if(val) projectionFilter = json_object_get_boolean(val) ? PETSC_TRUE : PETSC_FALSE;
        // -------------------------------------------
        json_object_object_get_ex(params, "Beta", &val);
        if(val) beta = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Beta final", &val);
        if(val) betaFinal = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "Eta", &val);
        if(val) eta = json_object_get_double(val);
        // -------------------------------------------
        json_object_object_get_ex(params, "(Rxx, Ryx, Ryy, Rzx, Rzy, Rzz) [SI]", &valRMatrix);
        // -------------------------------------------
        json_object_object_get_ex(params, "(Wxx, Wyx, Wyy, Wzx, Wzy, Wzz)", &valWMatrix);
        // -------------------------------------------
        // Initialize radius vector (in units of unitLength!)
        if (!valRMatrix){
            if(!valR){
                InputError("Specifiy a filter radius ('rmin') or a filter matrix (==symmetric)!", rankRead);
            }
            else{
                for(int i = 0; i < 6; i++){
                    // Initiliazise all radius (Rxx, Ryx, ...) to rmin      
                    rminVec[i] = rmin / unitLength;   
                }  
            }
        }
        else{
            for(int i = 0; i < 6; i++){  
                rminVec[i] = json_object_get_double(json_object_array_get_idx(valRMatrix, i)) / unitLength;   
            }  
        }
        // Initialize weight vector
        for(int i = 0; i < 6; i++){
            wminVec[i] = !valWMatrix ? wmin : json_object_get_double(json_object_array_get_idx(valWMatrix, i));
        }
        // ------------------------------------------- 
        
    }
}

void  TopOpt::SendParams(MPI_Comm comm, PetscInt rankSend) {
    // Wrong results when using PETSC_INT, PETSC_REAL...
    MPI_Bcast(&unitLength, 1, MPI_DOUBLE, rankSend, comm);
    for (int i = 0; i < NUM_MAX_REGIONS; i++){
        MPI_Bcast(&boundsDR[i][0][0], 2, MPI_INT, rankSend, comm); 
        MPI_Bcast(&boundsDR[i][1][0], 2, MPI_INT, rankSend, comm); 
        MPI_Bcast(&boundsDR[i][2][0], 2, MPI_INT, rankSend, comm); 
        MPI_Bcast(&boundsOR[i][0][0], 2, MPI_INT, rankSend, comm); 
        MPI_Bcast(&boundsOR[i][1][0], 2, MPI_INT, rankSend, comm); 
        MPI_Bcast(&boundsOR[i][2][0], 2, MPI_INT, rankSend, comm); 
    }    
    MPI_Bcast(&volfrac, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&maxItr, 1, MPI_INT, rankSend, comm); 
    MPI_Bcast(&rmin, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&filter, 1, MPI_INT, rankSend, comm); 
    MPI_Bcast(&Densmin, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&Densmax, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&movlim, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&restart, 1, MPI_C_BOOL, rankSend, comm); 
    MPI_Bcast(&projectionFilter, 1, MPI_C_BOOL, rankSend, comm); 
    MPI_Bcast(&beta, 1, MPI_DOUBLE, rankSend, comm);
    MPI_Bcast(&betaFinal, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&eta, 1, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&rminVec[0], 6, MPI_DOUBLE, rankSend, comm); 
    MPI_Bcast(&wminVec[0], 6, MPI_DOUBLE, rankSend, comm); 

};

PetscErrorCode TopOpt::SetUp(PetscInt identifier, PetscInt numProcDR[3], PetscInt* cellsArrDR[3], PetscReal spacingSI[3]) {

    PetscErrorCode ierr;
    MPI_Comm_rank(PETSC_COMM_WORLD, &petscRank);
	MPI_Comm_size(PETSC_COMM_WORLD, &petscSize);
    identifierDR = identifier;
    for (int i = 0; i < 3; i++){spacing[i] = spacingSI[i] / unitLength;}
	ierr = SetUpMESH(numProcDR, cellsArrDR); CHKERRQ(ierr);
    ierr = SetUpOPT(); CHKERRQ(ierr);
    return (ierr);
}

PetscErrorCode TopOpt::SetUpMESH(PetscInt numProcDR[3], PetscInt* cellsArrDR[3]) {

    // Setup grid object inluding domain decomposition ussing the parameters from the FDTD code
    PetscErrorCode ierr;
    DMDALocalInfo info;  
    PetscInt numnodaldof; // Stored values per cell. 
    PetscInt stencilwidth;  // Determines the number of neighbored ghost cells being accessible.  
    DMBoundaryType bx;  // Determines boundary type along x direction. 
    DMBoundaryType by;  // Determines boundary type along y direction. 
    DMBoundaryType bz;  // Determines boundary type along z direction. 
    DMDAStencilType stype;  // Stencil type. 
    stencilwidth = 1;
    numnodaldof = 3; // (dens_x, dens_y, dens_z) for each yee cell
    bx = DM_BOUNDARY_NONE;
    by = DM_BOUNDARY_NONE;
    bz = DM_BOUNDARY_NONE;
    stype = DMDA_STENCIL_BOX; // maybe star type sufficient enough?
    int coordMax[3] = {0}; // max. coordinates of a cell

    auto ColSum = [&](int idx) -> int{
        // Column-wise sum of cellsArrDr (global variable)
        cellsDR[idx] = 0; 
        for(int i = 0; i < numProcDR[idx]; i++){
            cellsDR[idx] += *(cellsArrDR[idx] + i);
        }
        coordMax[idx] = spacing[idx] * (cellsDR[idx] - 1);
    };
    ColSum(0); // X
    ColSum(1); // Y
    ColSum(2); // Z

    // Create Yee Grid object
    // -> an object that will manage the communication of a 3D (x DOFs) regular array data 
    // that is distributed across some processors 
    ierr = DMDACreate3d(PETSC_COMM_WORLD, 
                        bx, by, bz, stype, 
                        // global dimension in each direction of the array
                        cellsDR[0], cellsDR[1], cellsDR[2], 
                        // corresponding number of processors in each dimension 
                        numProcDR[0], numProcDR[1], numProcDR[2], 
                        numnodaldof, stencilwidth, 
                        // arrays containing the number of cell in each chunk along x, y and z
                        cellsArrDR[0], cellsArrDR[1], cellsArrDR[2], 
                        &(da_yee_lattice)); 
                        CHKERRQ(ierr);

    ierr = DMSetUp(da_yee_lattice); CHKERRQ(ierr);
    ierr = DMDASetElementType(da_yee_lattice, DMDA_ELEMENT_Q1); CHKERRQ(ierr); // necessary?
    // Set the coordinates (assuming the origin is at the left, lower, front corner of the Yee grid)
    ierr = DMDASetUniformCoordinates(da_yee_lattice, 0., coordMax[0], 0., coordMax[1], 0., coordMax[2]); CHKERRQ(ierr);

    DMDAGetLocalInfo(da_yee_lattice, &info);
    startIdx[0] = info.xs;
    startIdx[1] = info.ys;
    startIdx[2] = info.zs;
    locLength[0] = info.xm;
    locLength[1] = info.ym;
    locLength[2] = info.zm;
    numGPr[0] = (info.gxs + info.gxm) - (info.xs + info.xm); 
	numGPr[1] = (info.gys + info.gym) - (info.ys + info.ym); 
	numGPr[2] = (info.gzs + info.gzm) - (info.zs + info.zm);
  
    return (ierr);
}

PetscErrorCode TopOpt::SetUpOPT() {

    PetscErrorCode ierr = 0;
    ierr = DMCreateGlobalVector(da_yee_lattice, &densPhys); CHKERRQ(ierr);

    // Total number of design variables for one component-density
    VecGetSize(densPhys, &n);

    // Print optimization paramteres
    PetscPrintf(PETSC_COMM_WORLD, "------------------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "Opt. Params for DR No. %d:\n", identifierDR);
    PetscPrintf(PETSC_COMM_WORLD, "------------------------\n");
    PetscPrintf(PETSC_COMM_WORLD, "n (3x No. of cells):      %i\n", n);
    PetscPrintf(PETSC_COMM_WORLD, "m (constraints):          %i\n", m);
    PetscPrintf(PETSC_COMM_WORLD, "filter (0=sens., 1=dens): %i\n", filter);
    PetscPrintf(PETSC_COMM_WORLD, "rmin:                     %f\n", rmin);
    PetscPrintf(PETSC_COMM_WORLD, "projectionFilter (0/1):   %i\n", projectionFilter);
    PetscPrintf(PETSC_COMM_WORLD, "beta:                     %f\n", beta);
    PetscPrintf(PETSC_COMM_WORLD, "betaFinal:                %f\n", betaFinal);
    PetscPrintf(PETSC_COMM_WORLD, "eta:                      %f\n", eta);
    PetscPrintf(PETSC_COMM_WORLD, "volfrac:                  %f\n", volfrac);
    PetscPrintf(PETSC_COMM_WORLD, "maxItr:                   %i\n", maxItr);
    PetscPrintf(PETSC_COMM_WORLD, "movlim:                   %f\n", movlim);
    

    // Allocate after input
    gx = new PetscScalar[m];
    if (filter == 0) {
        Densmin = 0.001; // Prevent division by zero in filter
    }

    // Allocate the optimization vectors
    ierr = VecDuplicate(densPhys, &dens); CHKERRQ(ierr);
    ierr = VecDuplicate(densPhys, &densTilde); CHKERRQ(ierr);
    ierr = VecDuplicate(dens, &dens_old);

    ierr = VecSet(dens, volfrac); CHKERRQ(ierr);      // Initialize to volfrac!
    ierr = VecSet(densTilde, volfrac); CHKERRQ(ierr); // Initialize to volfrac!
    ierr = VecSet(densPhys, volfrac); CHKERRQ(ierr);  // Initialize to volfrac! 
    ierr = VecSet(dens_old, volfrac);                 // Initialize to volfrac! 

    // Sensitivity vectors
    ierr = VecDuplicate(dens, &dfdx); CHKERRQ(ierr); // change "dx" to dens later !!!
    ierr = VecDuplicateVecs(dens, m, &dgdx); CHKERRQ(ierr); // change "dx" to dens later !!! 

    // Bounds and extend densmin/densmax to component-wise vectors
    VecDuplicate(dens, &densmin);
    VecDuplicate(dens, &densmax);

    return (ierr);
}

PetscErrorCode TopOpt::SaveDensities(std::string prefix, bool* saveVecs){
    
    PetscErrorCode ierr;                               
    Vec densities[3] = {dens, densTilde, densPhys};
    std::string strID = to_string(identifierDR);
    std::string strVecs[3] = {"Dens" + strID, "DensTilde" + strID, "DensPhys" + strID};

    for (int s = 0; s < 3; s++){
        if (saveVecs[s]){
            Vec densLoc;
            PetscScalar ****densLocArr;
            ierr = DMCreateLocalVector(da_yee_lattice, &densLoc); CHKERRQ(ierr);
            ierr = DMGlobalToLocalBegin(da_yee_lattice, densities[s], INSERT_VALUES, densLoc); CHKERRQ(ierr);
            ierr = DMGlobalToLocalEnd(da_yee_lattice, densities[s], INSERT_VALUES, densLoc); CHKERRQ(ierr);
            ierr = DMDAVecGetArrayDOF(da_yee_lattice, densLoc, &densLocArr); CHKERRQ(ierr);
            PetscScalar*** array[3];
            for (int i = 0; i < 3; i++){
                Allocate3DArray(array[i]);
            } 
            for (int i = 0; i < locLength[0]; i++) {
                for (int j = 0; j < locLength[1]; j++) {
                    for (int k = 0; k < locLength[2]; k++) {
                        for (int c = 0; c < 3; c++){
                            array[c][i][j][k] = densLocArr[startIdx[2] + k][startIdx[1] + j][startIdx[0] + i][c];
                        }
                    }
                }
            }
            SaveAsBinary3D(prefix + strVecs[s] + ".cev", array[0], array[1], array[2], PETSC_COMM_WORLD, 
                                   startIdx, cellsDR, locLength, NULL, false, true, false);    
            for (int i = 0; i < 3; i++){
                Deallocate3DArray(array[i]);
            } 
            ierr = DMDAVecRestoreArrayDOF(da_yee_lattice, densLoc, &densLocArr); CHKERRQ(ierr);
        } 
    }
    return (ierr);
}

PetscErrorCode TopOpt::InitializeDensityToSphere(PetscInt centerx, PetscInt centery, PetscInt centerz, 
                                                 PetscReal SphereRadius, bool staircasing) {

    /* 
    Creates a radial-symmetric density distribution with 
    max. radius $R = \text{SphereRadius}$ and $\rho(r) = (R-r) / R for r <= R$ and $\rho(r) = 0$ otherwise
	Create 4D array (3 spatial dimesions + dDof) from global vector it ("local array", does NOT include 
    ghost cells, which are not necessary here anyway), 
	insert index values, restore array / update global vector 
	*/
    PetscErrorCode ierr;
    PetscReal x, xf, y, yf, z, zf, rEx, rEy, rEz, SphereRadius_ul; 
    PetscReal ****ArrInsert; 
    DMDALocalInfo info;      

    ierr = VecSet(dens, 0.); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da_yee_lattice, &info); CHKERRQ(ierr); 
    ierr = DMDAVecGetArrayDOF(da_yee_lattice, dens, &ArrInsert); CHKERRQ(ierr);  // array indices in global dimensions!
    SphereRadius_ul = SphereRadius / unitLength;

    // Run through all elements in the mesh - should not include ghosts
    for (int i = info.xs; i < info.xs + info.xm; i++){	
        x = (i - centerx) * spacing[0];
        xf = x + 0.5 * spacing[0];

        for (int j = info.ys; j < info.ys + info.ym; j++){			
            y = (j - centery) * spacing[1];
            yf = y + 0.5 * spacing[1];

            for (int k = info.zs; k < info.zs + info.zm; k++){
                z = (k - centerz) * spacing[2];
                zf = z + 0.5 * spacing[2];						
                rEx = sqrt(xf * xf + y * y + z * z);
                rEy = sqrt(x * x + yf * yf + z * z);
                rEz = sqrt(x * x + y * y + zf * zf);

                if (staircasing){
                    if (rEx <= SphereRadius_ul && rEy <= SphereRadius_ul && rEz <= SphereRadius_ul){	
                        ArrInsert[k][j][i][0] = 1; // (SphereRadius_ul - rEx) / SphereRadius_ul; 
                        ArrInsert[k][j][i][1] = 1; // (SphereRadius_ul - rEy) / SphereRadius_ul; 
                        ArrInsert[k][j][i][2] = 1; // (SphereRadius_ul - rEz) / SphereRadius_ul;                   
                    }
                }	
                else{
                    if (rEx <= SphereRadius_ul){
                        ArrInsert[k][j][i][0] = 1; // (SphereRadius_ul - rEx) / SphereRadius_ul; 
                    }
                    if (rEy <= SphereRadius_ul){
                        ArrInsert[k][j][i][1] = 1; // (SphereRadius_ul - rEy) / SphereRadius_ul;
                    }
                    if (rEz <= SphereRadius_ul){
                        ArrInsert[k][j][i][2] = 1; // (SphereRadius_ul - rEz) / SphereRadius_ul;
                    }
                }				
            }
        }
    }
    CHKMEMQ;
    ierr = DMDAVecRestoreArrayDOF(da_yee_lattice, dens, &ArrInsert); CHKERRQ(ierr);
    
	return ierr;                                         
}

void TopOpt::Allocate3DArray(PetscScalar*** &arr){

    // Store 1 more cell per axis (ghost) to use io.h's SaveToBinary()
    int lx =  locLength[0] + 1;
    int ly =  locLength[1] + 1;
    int lz =  locLength[2] + 1;

    arr = new PetscScalar**[lx]; 
    for (int i = 0; i < lx; i++) {
        arr[i] = new PetscScalar*[ly];
        for (int j = 0; j < ly; j++) {
            arr[i][j] = new PetscScalar[lz]();
        }
    }
}

void TopOpt::Deallocate3DArray(PetscScalar*** arr){

    int lx =  locLength[0] + 1;
    int ly =  locLength[1] + 1;

    for (int i = 0; i < lx; i++) {
        for (int j = 0; j < ly; j++) {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}









PetscErrorCode TopOpt::AllocateMMAwithRestart(PetscInt* itr, MMA** mma) {

    PetscErrorCode ierr = 0;

    // Set MMA parameters (for multiple load cases)
    PetscScalar aMMA[m];
    PetscScalar cMMA[m];
    PetscScalar dMMA[m];
    for (PetscInt i = 0; i < m; i++) {
        aMMA[i] = 0.0;
        dMMA[i] = 0.0;
        cMMA[i] = 1000.0;
    }

    // Check if restart is desired
    restart                  = PETSC_FALSE;  // DEFAULT USES RESTART
    flip                     = PETSC_TRUE;  // BOOL to ensure that two dump streams are kept
    PetscBool onlyLoadDesign = PETSC_FALSE; // Default restarts everything

    // Get inputs
    PetscBool flg;
    char      filenameChar[PETSC_MAX_PATH_LEN];
    PetscOptionsGetBool(NULL, NULL, "-restart", &restart, &flg);
    PetscOptionsGetBool(NULL, NULL, "-onlyLoadDesign", &onlyLoadDesign, &flg);

    if (restart) {
        ierr = VecDuplicate(dens, &xo1); 
        CHKERRQ(ierr);
        ierr = VecDuplicate(dens, &xo2); 
        CHKERRQ(ierr);
        ierr = VecDuplicate(dens, &U); 
        CHKERRQ(ierr);
        ierr = VecDuplicate(dens, &L); 
        CHKERRQ(ierr);
    }

    // Determine the right place to write the new restart files
    std::string filenameWorkdir = "./";
    PetscOptionsGetString(NULL, NULL, "-workdir", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {
        filenameWorkdir = "";
        filenameWorkdir.append(filenameChar);
    }
    filename00    = filenameWorkdir;
    filename00Itr = filenameWorkdir;
    filename01    = filenameWorkdir;
    filename01Itr = filenameWorkdir;

    filename00.append("/Restart00.dat");
    filename00Itr.append("/Restart00_itr_f0.dat");
    filename01.append("/Restart01.dat");
    filename01Itr.append("/Restart01_itr_f0.dat");

    // Where to read the restart point from
    std::string restartFileVec = ""; // NO RESTART FILE !!!!!
    std::string restartFileItr = ""; // NO RESTART FILE !!!!!

    PetscOptionsGetString(NULL, NULL, "-restartFileVec", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {
        restartFileVec.append(filenameChar);
    }
    PetscOptionsGetString(NULL, NULL, "-restartFileItr", filenameChar, sizeof(filenameChar), &flg);
    if (flg) {
        restartFileItr.append(filenameChar);
    }

    // Which solution to use for restarting
    PetscPrintf(PETSC_COMM_WORLD, "MMA: \n");
    PetscPrintf(PETSC_COMM_WORLD, "# Continue from previous iteration (-restart): %i \n", restart);
    PetscPrintf(PETSC_COMM_WORLD, "# Restart file (-restartFileVec): %s \n", restartFileVec.c_str());
    PetscPrintf(PETSC_COMM_WORLD, "# Restart file (-restartFileItr): %s \n", restartFileItr.c_str());
    PetscPrintf(PETSC_COMM_WORLD,
                "# New restart files are written to (-workdir): %s "
                "(Restart0x.dat and Restart0x_itr_f0.dat) \n",
                filenameWorkdir.c_str());

    // Check if files exist:
    PetscBool vecFile = fexists(restartFileVec);
    if (!vecFile) {
        PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileVec.c_str());
    }
    PetscBool itrFile = fexists(restartFileItr);
    if (!itrFile) {
        PetscPrintf(PETSC_COMM_WORLD, "File: %s NOT FOUND \n", restartFileItr.c_str());
    }

    // Read from restart point

    PetscInt nGlobalDesignVar;
    VecGetSize(dens,
               &nGlobalDesignVar); // ASSUMES THAT SIZE IS ALWAYS MATCHED TO CURRENT MESH
    if (restart && vecFile && itrFile) {

        PetscViewer view;
        // Open the data files
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, restartFileVec.c_str(), FILE_MODE_READ, &view);

        VecLoad(dens, view);
        VecLoad(densPhys, view);
        VecLoad(xo1, view);
        VecLoad(xo2, view);
        VecLoad(U, view);
        VecLoad(L, view);
        PetscViewerDestroy(&view);

        // Read iteration and fscale
        std::fstream itrfile(restartFileItr.c_str(), std::ios_base::in);
        itrfile >> itr[0];
        itrfile >> fscale;

        // Choose if restart is full or just an initial design guess
        if (onlyLoadDesign) {
            PetscPrintf(PETSC_COMM_WORLD, "# Loading design from file: %s \n", restartFileVec.c_str());
            *mma = new MMA(nGlobalDesignVar, m, dens, aMMA, cMMA, dMMA); 
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "# Continue optimization from file: %s \n", restartFileVec.c_str());
            *mma = new MMA(nGlobalDesignVar, m, *itr, xo1, xo2, U, L, aMMA, cMMA, dMMA);
        }

        PetscPrintf(PETSC_COMM_WORLD, "# Successful restart from file: %s and %s \n", restartFileVec.c_str(),
                    restartFileItr.c_str());
    } else {
        *mma = new MMA(nGlobalDesignVar, m, dens, aMMA, cMMA, dMMA); 
    }

    return ierr;
}

PetscErrorCode TopOpt::WriteRestartFiles(PetscInt* itr, MMA* mma) {

    PetscErrorCode ierr = 0;
    // Only dump data if correct allocater has been used
    if (!restart) {
        return -1;
    }

    // Get restart vectors
    mma->Restart(xo1, xo2, U, L);

    // Choose previous set of restart files
    if (flip) {
        flip = PETSC_FALSE;
    } else {
        flip = PETSC_TRUE;
    }

    // Write file with iteration number of f0 scaling
    // and a file with the MMA-required vectors, in the following order:
    // : x,xPhys,xold1,xold2,U,L
    PetscViewer view;         // vectors
    PetscViewer restartItrF0; // scalars

    PetscViewerCreate(PETSC_COMM_WORLD, &restartItrF0);
    PetscViewerSetType(restartItrF0, PETSCVIEWERASCII);
    PetscViewerFileSetMode(restartItrF0, FILE_MODE_WRITE);

    // Open viewers for writing
    if (!flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename00.c_str(), FILE_MODE_WRITE, &view);
        PetscViewerFileSetName(restartItrF0, filename00Itr.c_str());
    } else if (flip) {
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename01.c_str(), FILE_MODE_WRITE, &view);
        PetscViewerFileSetName(restartItrF0, filename01Itr.c_str());
    }

    // Write iteration and fscale
    PetscViewerASCIIPrintf(restartItrF0, "%d ", itr[0]);
    PetscViewerASCIIPrintf(restartItrF0, " %e", fscale);
    PetscViewerASCIIPrintf(restartItrF0, "\n");

    // Write vectors
    VecView(dens, view); 
    VecView(densPhys, view); 
    VecView(xo1, view);
    VecView(xo2, view);
    VecView(U, view);
    VecView(L, view);

    // Clean up
    PetscViewerDestroy(&view);
    PetscViewerDestroy(&restartItrF0);

    // PetscPrintf(PETSC_COMM_WORLD,"DONE WRITING DATA\n");
    return ierr;
}
