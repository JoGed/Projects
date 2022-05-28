#include "externs.h" // contains Filter.h, TopOpt.h, MMA.h
#include "OptFunctions.h"
#include "RegionEmbedding.h"
#include "io.h"
#include <algorithm>

/*
This file provides multiple functions needed to setup and embed Topology Optimization into the existing FDTD code
(by using its global variables), e.g. decomposition of Design and Observation regions based in the FDTD domain 
decomposition and assigning "optimizer-ranks" which handle the optimization, as well as the computation of 
the objective function, the gradients, injecting the adjoint source, monitoring change of objectives, 
computing mean frequency spectrum for OR's and more...
*/


// -------------------------------------------------------------------------------------------------------
/* 
Here, TopOpt object will be initalized by ALL processes first (even if they do not handle
the optimization later). ROOTRANK will read and send the TopOpt parameters (from a json file) 
to all the other ranks 
*/
void ReadAndSendParams(){
	// Here, ALL processes initialize TopOpt object. Will be sorted out later
	opt = new TopOpt(); 
	// Root process reads TopOpt parameters from file and sends them to the others
	if (RankThis == ROOTPROCESS)
	{	
		opt->ReadParamsFromFile("./INIDEF/paramsTopOpt.json", ROOTPROCESS);
		obsDecayWindow = opt->decayWindow;
		obsDecayMinimum = opt->decayCond;
	}
	opt->SendParams(MPI_COMM_WORLD, ROOTPROCESS); 
	MPI_Bcast(&obsDecayWindow, 1, MPI_INT, ROOTPROCESS, MPI_COMM_WORLD);
	MPI_Bcast(&obsDecayMinimum, 1, MPIREAL, ROOTPROCESS, MPI_COMM_WORLD);
	// Can not intialize strings right after the following ROOTPRROCESS routine...?	
	strRegions[0] = "Observation Region";
  strRegions[1] = "Design Region"; 
	// Get number of max TopOpt iterations
	maxOptItr = opt->maxItr;
	if(RankThis == ROOTPROCESS){
		printf("------------------------------------\n");
		printf("Program is running in optimization mode!\n");
		printf("------------------------------------\n");
	}
}

// -------------------------------------------------------------------------------------------------------
/*
Check if Design- and Observation regions are located inside Total Field Region
*/
void CheckInput(){
	// Get the FDTD decomposition info from global variables
	// -------------------------------------------
	// Global FDTD boundary indices in x (indluding PMLs)
  boundsChunkFDTD[0][0] = NRCPMLPADNEAX;
  boundsChunkFDTD[0][1] = NRCPMLPADNEAX + NRCELLSXthis - 1;
  boundsChunkFDTD[1][0] = NRCPMLPADNEAY;
  boundsChunkFDTD[1][1] = NRCPMLPADNEAY + NRCELLSYthis - 1;
  boundsChunkFDTD[2][0] = NRCPMLPADNEAZ;
  boundsChunkFDTD[2][1] = NRCPMLPADNEAZ + NRCELLSZthis - 1;

	// Left  Total Field boundary idx in x: pids->NRCELLSCPMLX + pids->NRINFLATEX - 1;
	// Right Total Field boundary idx in x: pids->NRCELLSTOTX - pids->NRCELLSCPMLX - pids->NRINFLATEX - 1;
	int boundsTF[3][2] = {{pids->NRCELLSCPMLX + pids->NRINFLATEX - 1, 
                         pids->NRCELLSTOTX - pids->NRCELLSCPMLX - pids->NRINFLATEX - 1}, 
                        {pids->NRCELLSCPMLY + pids->NRINFLATEY - 1, 
                         pids->NRCELLSTOTY - pids->NRCELLSCPMLY - pids->NRINFLATEY - 1}, 
	                      {pids->NRCELLSCPMLZ + pids->NRINFLATEZ - 1, 
                         pids->NRCELLSTOTZ - pids->NRCELLSCPMLZ - pids->NRINFLATEZ - 1}}; 
	// Length array of local FDTD subdomain														 
	lengthChunkFDTD[0] = NRCELLSXthis;
  lengthChunkFDTD[1] = NRCELLSYthis;
  lengthChunkFDTD[2] = NRCELLSZthis;
	// Spacing in SI units
	spacing[0] = DX;
  spacing[1] = DY; 
  spacing[2] = DZ;
	// -------------------------------------------
	// 1.) Check if regions boundaries are valid
	// -------------------------------------------
	for (int i = 0; i < NRMAXREGIONS; i++){
		TopOpt::DRisDefined[i] = (TopOpt::boundsDR[i][0][1] - TopOpt::boundsDR[i][0][0]) > 0
		                      && (TopOpt::boundsDR[i][1][1] - TopOpt::boundsDR[i][1][0]) > 0
													&& (TopOpt::boundsDR[i][2][1] - TopOpt::boundsDR[i][2][0]) > 0;
		// -------------------------------------------
		// 1.a) Check if design region intersects observation region
		// ------------------------------------------- 

		for (int j = 0; j < NRMAXREGIONS; j++){
			if (i == 0){
				TopOpt::ORisDefined[j] = (TopOpt::boundsOR[j][0][1] - TopOpt::boundsOR[j][0][0]) > 0
		                          && (TopOpt::boundsOR[j][1][1] - TopOpt::boundsOR[j][1][0]) > 0
													    && (TopOpt::boundsOR[j][2][1] - TopOpt::boundsOR[j][2][0]) > 0;									
			}
			// if (Intersection(TopOpt::boundsOR[j], TopOpt::boundsDR[i]) && TopOpt::ORisDefined[j] && TopOpt::DRisDefined[i]){
			// 	opt->InputError("Observation Region intersects Design Region!", ROOTPROCESS);
			// }
		}
		// -------------------------------------------
		// 1.b) Check if both design-and observation region are located inside Total Field Region
		// ------------------------------------------- 
		// if (!(Inside(TopOpt::boundsDR[i], boundsTF)) && TopOpt::DRisDefined[i]){
		// 	opt->InputError("Design Region must be located Inside the Total Field Region!", ROOTPROCESS);
		// }
		if (!(Inside(TopOpt::boundsOR[i], boundsTF)) && TopOpt::ORisDefined[i]){
			opt->InputError("Observation Region must be located inside the Total Field Region!", ROOTPROCESS);
		} 
    // printf("\n------------------\n%d, %d, %d, %d, %d, %d \n %d, %d, %d, %d, %d, %d\n", 
		// 																		TopOpt::boundsOR[0][0][0], TopOpt::boundsOR[0][0][1],
	  //                                     TopOpt::boundsOR[0][1][0], TopOpt::boundsOR[0][1][1], 
		// 																		TopOpt::boundsOR[0][2][0], TopOpt::boundsOR[0][2][1],
		// 																		boundsTF[0][0], boundsTF[0][1],
		// 																		boundsTF[1][0], boundsTF[1][1],
		// 																		boundsTF[2][0], boundsTF[2][1]);
		// -------------------------------------------
		//  Check if design region intersects other design region
		// ------------------------------------------- 
		// -> will be checked later while PETSc communicator is being initialized!
	}
}

// -------------------------------------------------------------------------------------------------------
/*
Setup Observation Regions and determine common rank "SendRankObjUnion" which 
will gather objective values and send it to all optimizer (DR) ranks later
*/
void InitializeOR(){
  // ------------------------
	// OBSERVATION REGIONS
	// ------------------------
	// Create Communicator for each Observation Region
	int ranksNumOR;
	noIntrsORAtAll = true; // global variable
	for (int r = 0; r < NRMAXREGIONS; r++){
		if (TopOpt::ORisDefined[r]){
			intrsOR[r] = 1 ? (Intersection(TopOpt::boundsOR[r], boundsChunkFDTD)) : 0;
		}
		MPI_Comm_split(MPI_COMM_WORLD, intrsOR[r], 0, &OBJ_COMM[r]);
		if(intrsOR[r]){
			noIntrsORAtAll = false;
			MPI_Comm_rank(OBJ_COMM[r], &ObjCommRank[r]);
			MPI_Comm_size(OBJ_COMM[r], &ObjCommSize[r]);
		}
	}
	// -------------------------------------------
	// Declare/Initallize collective communicator group/comm
	OBJ_COMM_UNION = MPI_COMM_NULL; // glob. var
	MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP); // glob. var
	MPI_Group OBJ_COMM_UNION_GROUP; 
	// -------------------------------------------
	// Get decomposition variables
	for (int r = 0; r < NRMAXREGIONS; r++){
		if (TopOpt::ORisDefined[r]){	
			// ------------------------
			// Get intersecting ranks
			int ranksNumOR;
			MPI_Allreduce(&intrsOR[r], &ranksNumOR, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			int* processesOR = new int[ranksNumOR];
			GetRanks(MPI_COMM_WORLD, intrsOR[r], &(processesOR[0]));
			// ------------------------
			// Update collective communicator group
			MPI_Group tmpGroupOR;	
			MPI_Group_incl(WORLD_GROUP, ranksNumOR, processesOR, &tmpGroupOR);
			if (r == 0) OBJ_COMM_UNION_GROUP = tmpGroupOR; // If only 1 region is defined
			MPI_Group_union(OBJ_COMM_UNION_GROUP, tmpGroupOR, &OBJ_COMM_UNION_GROUP);
			MPI_Group_free(&tmpGroupOR);
			// ------------------------
			// Get additional decomp. info and print all embedding parameters for verification
			int numProcOR[3] = {};
			int ObjRoot = 0;
			if (intrsOR[r]){ 
				GetShifts(TopOpt::boundsOR[r], boundsChunkFDTD, &shiftOR[r][0]);
				GetLocalLength(TopOpt::boundsOR[r], boundsChunkFDTD, lengthChunkFDTD, &locLengthOR[r][0]);
				GetStartIndices(TopOpt::boundsOR[r], boundsChunkFDTD, &startIdxOR[r][0]); 
				CountRanksPerAxis(OBJ_COMM[r], startIdxOR[r], &numProcOR[0]);  
				int* cellsArrOR[3];
				for (int i = 0; i < 3; i++){ 
					cellsArrOR[i] = new int[numProcOR[i]]; 
				}			
				GetCellArrayPerAxis(OBJ_COMM[r], locLengthOR[r], ranksNumOR, numProcOR, &cellsArrOR[0]);
				if (ObjCommRank[r] == ObjRoot){
					PrintDecomposition(ObjCommSize[r], GlobalProcsSize, processesOR, numProcOR, cellsArrOR, "Observation Region " + to_string(r));
				} 
				// Clean up
				for (int i = 0; i < 3; i++){ delete [] cellsArrOR[i]; }
			} 
			MPI_Comm_free(&OBJ_COMM[r]);
			// Save 1st rank of 1st OR (must be always defined anyways) as SendRankObjUnion
			if (r == 0) SendRankObjUnion = processesOR[0];
			// Clean up
			delete [] processesOR; 
		}
	}
	// -------------------------------------------
	// Create collective observation communicator
	MPI_Comm_create(MPI_COMM_WORLD, OBJ_COMM_UNION_GROUP, &OBJ_COMM_UNION);
	MPI_Group_free(&OBJ_COMM_UNION_GROUP);
	// Verify union
	if (OBJ_COMM_UNION != MPI_COMM_NULL){	
		MPI_Comm_rank(OBJ_COMM_UNION, &ObjCommUnionRank);
		MPI_Comm_size(OBJ_COMM_UNION, &ObjCommUnionSize);
		// Verify union stuff
		if (ObjCommUnionRank == 0){
			int* processesUnion = new int[ObjCommUnionSize];
			MPI_Gather(&RankThis, 1, MPI_INT, &processesUnion[0], 1, MPI_INT, 0, OBJ_COMM_UNION);
			sort(processesUnion, processesUnion + ObjCommUnionSize);
			printf("\n------------------------------------\n");
			printf("Ranks sharing Observation Regions in total (%d ranks):\n{", ObjCommUnionSize); 
  		for (int i = 0; i < ObjCommUnionSize; i++){ 
    		printf("%i, ", processesUnion[i]); 
			}
			printf("}\n\n");
			delete[] processesUnion;
		}
		else{
			MPI_Gather(&RankThis, 1, MPI_INT, NULL, 0, MPI_INT, 0, OBJ_COMM_UNION);
		}
	}
	else{
		ObjCommUnionRank = -1;
		ObjCommUnionSize = 0;
	}
  
	// -------------------------------------------
	// Get rank index of SendRankObjUnion in OBJ_COMM_UNION
	if (RankThis == SendRankObjUnion) ROOT_OBJ_COMM_UNION = ObjCommUnionRank;
	MPI_Bcast(&ROOT_OBJ_COMM_UNION, 1, MPI_INT, SendRankObjUnion, MPI_COMM_WORLD);
	// -------------------------------------------
	//  Number of right ghost points of local subdomain of observation regions (needed to update adjoint field)
	if (!noIntrsORAtAll){
		for (int r = 0; r < NRMAXREGIONS; r++){ 
			if (intrsOR[r]){
				for (int i = 0; i < 3; i++){
					// Needed in Efield-Updating function
					globBoundsOR[r][i][0] = TopOpt::boundsOR[r][i][0];
					globBoundsOR[r][i][1] = TopOpt::boundsOR[r][i][1];
					numGPrOR[r][i] = 1 ? ((startIdxOR[r][i] + locLengthOR[r][i]) <= TopOpt::boundsOR[r][i][1]) : 0;
				}
			}
		}
	} 
}

// -------------------------------------------------------------------------------------------------------
/*
Setup Design Region(s) and initialize TopOpt-, Filter- and MMA-objects
*/
void InitializeDR(){
  // -------------------
	// DESIGN REGIONS 
	// -------------------
	// Communicator will be defined differently than before, PETSC_COMM_WORLD is unique (no array of comms used)
	// Intersection of FDTD chunk with > 1 DR must be avoided!
	// -------------------------------------------
	// Check if design region intersects other design region
	// -------------------------------------------
	int intrsDRTot = 0;
	int colorDR = -1; // Assign integer to DR's for communicato splitting
	for (int r = 0; r < NRMAXREGIONS; r++){
		if (TopOpt::DRisDefined[r]){
			intrsDR[r] = 1 ? (Intersection(TopOpt::boundsDR[r], boundsChunkFDTD)) : 0;
			intrsDRTot += intrsDR[r];
			if (intrsDRTot > 1){
				opt->InputError("FDTD subdomain contains cells from different design Regions!", RankThis);
			}
			if (intrsDR[r]){
				colorDR = r;
			}
		}
	}
	MPI_Comm_split(MPI_COMM_WORLD, colorDR, 0, &PETSC_COMM_WORLD);
	//++++++++++++++++++++++++++++++++++++
	// Initialize PETSc communicator here!
	PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);  
	//++++++++++++++++++++++++++++++++++++ 
	MPI_Comm_rank(PETSC_COMM_WORLD, &PETScRank);
	MPI_Comm_size(PETSC_COMM_WORLD, &PETScCommSize);
	//printf("WORLD RANK/SIZE: %d/%d \t PETSC RANK/SIZE: %d/%d\n", RankThis, GlobalProcsSize, PETScRank, PETScCommSize); 	
	// -------------------------------------------
	OBJ_DR_COMM = MPI_COMM_NULL; // glob. var
	MPI_Comm_group(MPI_COMM_WORLD, &WORLD_GROUP); // glob. var
	processesDR = new int*[NRMAXREGIONS];
	ranksNumDR = new int[NRMAXREGIONS];
	DR_defined = 0;
	for (int r = 0; r < NRMAXREGIONS; r++){
		processesDR[r] = NULL;
		ranksNumDR[r] = 0;
		if (TopOpt::DRisDefined[r]){
			DR_defined++;
			MPI_Allreduce(&intrsDR[r], &ranksNumDR[r], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			processesDR[r] = new int[ranksNumDR[r]];
			GetRanks(MPI_COMM_WORLD, intrsDR[r], &(processesDR[r][0]));
			if (intrsDR[r]){ 
				int numProcDR[3] = {};	 
				int PETScRoot = 0;
				GetShifts(TopOpt::boundsDR[r], boundsChunkFDTD, &shiftDR[r][0]);
				GetLocalLength(TopOpt::boundsDR[r], boundsChunkFDTD, lengthChunkFDTD, &locLengthDR[r][0]);
				GetStartIndices(TopOpt::boundsDR[r], boundsChunkFDTD, &startIdxDR[r][0]);
				CountRanksPerAxis(PETSC_COMM_WORLD, startIdxDR[r], &numProcDR[0]); 
				int* cellsArrDR[3];
				for (int i = 0; i < 3; i++){ cellsArrDR[i] = new int[numProcDR[i]]{}; }
				GetCellArrayPerAxis(PETSC_COMM_WORLD, locLengthDR[r], ranksNumDR[r], numProcDR, &cellsArrDR[0]);
				if (PETScRank == PETScRoot){
					PrintDecomposition(PETScCommSize, GlobalProcsSize, processesDR[r], numProcDR, cellsArrDR, "Design Region " + to_string(r));
				} 
				// -------------------------------------------
				// SETUP OPTIMIZATION PARAMETERS, DATA AND LATTICE
				opt->SetUp(r, numProcDR, cellsArrDR, spacing);
				// -------------------------------------------
				// INITIALIZE DENSITY (+ STORE DENSITIES INTO BINARY)
				bool saveVecs[3] = {true, false, false}; // {SaveDens, SaveDensTilde, SaveDensPhys};
				// Define length of normal vectors perpendicular to XY, XZ, YZ planes
				int normalVecIdx[3] = {(TopOpt::boundsDR[r][2][1] - TopOpt::boundsDR[r][2][0] + 1) / 2,
															 (TopOpt::boundsDR[r][1][1] - TopOpt::boundsDR[r][1][0] + 1) / 2, 
															 (TopOpt::boundsDR[r][0][1] - TopOpt::boundsDR[r][0][0] + 1) / 2};

				// opt->InitializeDensityToSphere(normalVecIdx[2], normalVecIdx[1], normalVecIdx[0], 
				// 															min(DZ * normalVecIdx[2], min(DY * normalVecIdx[1], DX *normalVecIdx[0])), 
				// 															STAIRCASING);
				// -------------------------------------------
				// THE FILTER
				filter = new Filter(opt->da_yee_lattice, opt->dens, opt->filter, opt->rminVec, opt->wminVec);
				// -------------------------------------------
				// THE OPTIMIZER MMA
				opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart ! restart == false
				// -------------------------------------------
				// FILTER THE INTIAL DENSITY (without Observation Punch out)
			 	filter->FilterProject(opt->dens, opt->densTilde, opt->densPhys, opt->projectionFilter, opt->beta, opt->eta);
				// -------------------------------------------
				// CLEAN UP
				for (int i = 0; i < 3; i++){ delete [] cellsArrDR[i]; }
			}
		}
	}
  // -------------------------------------------
	// Save noIntrsDRAtAll flag, DR bounds + rGPs as global vars -> needed for E update
	noIntrsDRAtAll = true; // global variable
	for (int r = 0; r < NRMAXREGIONS; r++){ 	
		if (intrsDR[r]){
			for (int i = 0; i < 3; i++){
				// Needed in Efield-Updating function
				globBoundsDR[r][i][0] = TopOpt::boundsDR[r][i][0];
				globBoundsDR[r][i][1] = TopOpt::boundsDR[r][i][1];
				numGPrDR[r][i] = opt->numGPr[i];
			}
			noIntrsDRAtAll = false;
			break;
		}
		else{
			for (int i = 0; i < 3; i++){
				// Set "No defined values" for ranks not being involved in optimizing
				numGPrDR[r][i] = startIdxDR[r][i] = shiftDR[r][i] = locLengthDR[r][i] = -1;
			}
		}
	}
}

// -------------------------------------------------------------------------------------------------------
/*
Create communicator unifying root rank of united DR Comm. and PETSc Comm.
That is needed for sending the (collected) value of objectiv function to DR ranks
*/
void SetupORDRCommunication(){
	// Get all DR processes
	int ranksNumDR_sum = 0; 
	for (int i = 0; i < DR_defined; i++) ranksNumDR_sum += ranksNumDR[i];
	processesDR_total = new int[ranksNumDR_sum];
	for (int i = 0; i < DR_defined; i++){
		// processesDR's are distinct by definition
		for (int j = 0; j < ranksNumDR[i]; j++){
			int newRegIdx = 0 ? i == 0 : ranksNumDR[i - 1];
			processesDR_total[newRegIdx + j] = processesDR[i][j];
		}
	}
	// Built communicator
	MPI_Group ObjSendingGroup; // will consist only 1 rank
	MPI_Group DRGroup;	
	MPI_Group OBJ_DR_COMM_GROUP;
	MPI_Group_incl(WORLD_GROUP, 1, &SendRankObjUnion, &ObjSendingGroup);
	MPI_Group_incl(WORLD_GROUP, ranksNumDR_sum, processesDR_total, &DRGroup);
	MPI_Group_union(ObjSendingGroup, DRGroup, &OBJ_DR_COMM_GROUP);
	MPI_Comm_create(MPI_COMM_WORLD, OBJ_DR_COMM_GROUP, &OBJ_DR_COMM);
	MPI_Group_free(&ObjSendingGroup);
	MPI_Group_free(&DRGroup);
	MPI_Group_free(&OBJ_DR_COMM_GROUP);
	if (OBJ_DR_COMM != MPI_COMM_NULL){
		MPI_Comm_rank(OBJ_DR_COMM, &ObjPetscUnitCommRank);
		MPI_Comm_size(OBJ_DR_COMM, &ObjPetscUnitCommSize);
	}
	else{
		ObjPetscUnitCommRank = -1;
		ObjPetscUnitCommSize = 0;
	}
	// -------------------------------------------
	// Get rank index of SendRankObjUnion in OBJ_DR_COMM
	if (RankThis == SendRankObjUnion) ROOT_OBJ_DR_COMM = ObjPetscUnitCommRank;
	MPI_Bcast(&ROOT_OBJ_DR_COMM, 1, MPI_INT, SendRankObjUnion, MPI_COMM_WORLD);
	// -------------------------------------------
	// Verify
	if (OBJ_DR_COMM != MPI_COMM_NULL){
			int* processesUnion = new int[ObjPetscUnitCommSize];
			MPI_Gather(&RankThis, 1, MPI_INT, &processesUnion[0], 1, MPI_INT, 0, OBJ_DR_COMM);
			sort(processesUnion, processesUnion + ObjPetscUnitCommSize);
			if (ObjPetscUnitCommRank == 0){
				printf("\n------------------------------------\n");
				printf("Ranks sharing DR Regions in total + Observation sending rank (%d ranks):\n{", ObjPetscUnitCommSize); 
				for (int i = 0; i < ObjPetscUnitCommSize; i++){ 
				printf("%i, ", processesUnion[i]); 
				}
				printf("}\n\n\n");
			}
			delete[] processesUnion;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

// -------------------------------------------------------------------------------------------------------
/* 
Delete unneeded arrays and finalize PETCs in the case 
the rank does not intersects DR at all
*/
void CleanUpDRInit(){
	for (int i = 0; i < NRMAXREGIONS; i++) if (processesDR[i] != NULL) delete processesDR[i];
	delete [] processesDR_total;
	delete [] processesDR;
	delete [] ranksNumDR;
	if(noIntrsDRAtAll){
		delete opt;
		PETScCommSize = 0;
		PETScRank = -1;
		PetscFinalize(); // -> PETSC_COMM_WORLD = MPI_COMM_NULL
	}
}

// -------------------------------------------------------------------------------------------------------
/* 
Set density and density_phys values to 0 in observation regions
*/
void ZeroDensInOR(int r){
  if (intrsDR[r]){
    for (int o = 0; o < NRMAXREGIONS; o++){
      if(intrsOR[o]){												
        PetscReal ****dens_arr;	
        PetscReal ****densPhys_arr;	
        DMDAVecGetArrayDOF(opt->da_yee_lattice, opt->dens, &dens_arr); // array indices in global dimensions!
        DMDAVecGetArrayDOF(opt->da_yee_lattice, opt->densPhys, &densPhys_arr); // array indices in global dimensions!
        int petscStartIdx[3] = {startIdxDR[r][0] - globBoundsDR[r][0][0], 
                                startIdxDR[r][1] - globBoundsDR[r][1][0], 
                                startIdxDR[r][2] - globBoundsDR[r][2][0]};		
          int relDist[3] = {startIdxOR[o][0] - startIdxDR[r][0],
                            startIdxOR[o][1] - startIdxDR[r][1],
                            startIdxOR[o][2] - startIdxDR[r][2]};						
                            
          for (int u = 0; u < locLengthDR[r][0]; u++){
            if (u >= relDist[0] && u < (relDist[0] + locLengthOR[o][0])){
              for (int v = 0; v < locLengthDR[r][1]; v++){
                if (v >= relDist[1] && v < (relDist[1] + locLengthOR[o][1])){
                  for (int w = 0; w < locLengthDR[r][2]; w++){	
                    if (w >= relDist[2] && w < (relDist[2] + locLengthOR[o][2])){	
                      dens_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][0] = 0;
                      dens_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][1] = 0;
                      dens_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][2] = 0;
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][0] = 0;
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][1] = 0;
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][2] = 0;
                    }
                  }
                }
              }
            }
          }
        DMDAVecRestoreArrayDOF(opt->da_yee_lattice, opt->dens, &dens_arr);
        DMDAVecRestoreArrayDOF(opt->da_yee_lattice, opt->densPhys, &densPhys_arr);
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------------
/*
Initialize global simulation parameters
*/
void InitializeSimParams(){
	adjoint_mode = false;
	itr = 0;
	num_sweep = 2 * maxOptItr; // Number of FDTD simulations (== 2 x opt. Iterations)
	preSim = true;
	if (!noIntrsDRAtAll){
		objConstant = 10;
	}
	PetscReal objFunLoc_window = 0;
	integrationTime = NRMAXTIMESTEPS; 
	abortSim = false;
	measureCount = 1;
	fx_window_old = 0;
	optStartTime = MPI_Wtime();
  signFac = 1;
}

// -------------------------------------------------------------------------------------------------------
/*
Save design, observation region(s) and material indices into a single binary file (for verification)
*/
void SaveRegions(){
  int*** regionsIdx; // DR == 5 , OR == 6, 0 == vacuum, else: materialIndex
			MemAllocat3DArray(regionsIdx, NRCELLSXthis, NRCELLSYthis, NRCELLSZthis);
			// BACKGROUND (MATERIAL INDEX)
			if (!pids->CPMLONLY[RankThis]){
				for(int i = 0; i < NRXinnerthis; i++){
					for(int j = 0; j < NRYinnerthis; j++){
						for(int k = 0; k < NRZinnerthis; k++){
							// Using Z comp. of  material here!
							regionsIdx[i + NRSHIFTX][j + NRSHIFTY][k + NRSHIFTZ] = matindz[i][j][k];
						}
					}	
				}
			}
			for (int r = 0; r < NRMAXREGIONS; r++){
				// DESIGN REGION
				if(intrsDR[r]){
					for(int i = 0; i < locLengthDR[r][0]; i++){
						for(int j = 0; j < locLengthDR[r][1]; j++){
							for(int k = 0; k < locLengthDR[r][2]; k++){
								regionsIdx[i + shiftDR[r][0]][j + shiftDR[r][1]][k + shiftDR[r][2]] = 5;
							}
						}	
					}
				} 
				// OBSERVATION REGION
				if(intrsOR[r]){
					for(int i = 0; i < locLengthOR[r][0]; i++){
						for(int j = 0; j < locLengthOR[r][1]; j++){
							for(int k = 0; k < locLengthOR[r][2]; k++){
									regionsIdx[i + shiftOR[r][0]][j + shiftOR[r][1]][k + shiftOR[r][2]] = 6;
							}
						}	
					}
				}		
			} 
			int totSize[3] = {NRTOTX, NRTOTY, NRTOTZ};
			int totIdx[3] = {NRCPMLPADNEAX, NRCPMLPADNEAY, NRCPMLPADNEAZ}; 
			int locLength[3] = {NRCELLSXthis, NRCELLSYthis, NRCELLSZthis};
			//Save regions map to binary for verification
			SaveAsBinary3D("Regions.cev", regionsIdx, regionsIdx, regionsIdx, MPI_COMM_WORLD, 
																		totIdx, totSize, locLength, NULL, false, true, false);
			MemRelease3DArray(regionsIdx, NRCELLSXthis, NRCELLSYthis, NRCELLSZthis);
}

// -------------------------------------------------------------------------------------------------------
/*
Injection of adjoint source based on enhancement objective function into r-th OR
*/
void AdjSource(int r){
  // Example x-component:
  // grad{Ex} = -2 * eps_x * Ex_forward(T-t) on the volume of interest, Ex_forward = forward x field		

  if (intrsOR[r]){						
    // Get materialindex startindices
    int mIdx[3] = {- NRSHIFTX + shiftOR[r][0],
                    - NRSHIFTY + shiftOR[r][1],
                    - NRSHIFTZ + shiftOR[r][2]};	
    REAL prefac_x, prefac_y, prefac_z;	
    // Ghost points needed here 
    for (int u = 0; u < locLengthOR[r][0] + numGPrOR[r][0]; u++){
      for (int v = 0; v < locLengthOR[r][1] + numGPrOR[r][1]; v++){
        for (int w = 0; w < locLengthOR[r][2] + numGPrOR[r][2]; w++){
          prefac_x = EP0 * EPRE[matindx[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]] 
                    * CE2[matindx[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]];
          prefac_y = EP0 * EPRE[matindy[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]] 
                    * CE2[matindy[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]];
          prefac_z = EP0 * EPRE[matindz[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]] 
                    * CE2[matindz[u + mIdx[0]][v + mIdx[1]][w + mIdx[2]]];
          Ex[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] += 
            signFac * prefac_x * Ex_timedep_OR[r][maxTimeItr - timestep][u][v][w];
          Ey[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] += 
            signFac * prefac_y * Ey_timedep_OR[r][maxTimeItr - timestep][u][v][w];
          Ez[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] += 
            signFac * prefac_z * Ez_timedep_OR[r][maxTimeItr - timestep][u][v][w];						
        }
      }
    }  
  }
}

// -------------------------------------------------------------------------------------------------------
/*
Update gradients (in r-th DR) based on enhancement objective function
Formulas for non-dispersive materials used here
See paper "Topology optimization of dispersive plasmonic nanostructures in the time-domain", Eq. (25b)
*/
void UpdateGradients(int r){
  // Ghost points NOT needed here!
  int petscStartIdx[3] = {startIdxDR[r][0] - globBoundsDR[r][0][0], 
                          startIdxDR[r][1] - globBoundsDR[r][1][0], 
                          startIdxDR[r][2] - globBoundsDR[r][2][0]};	
                                                      
  for (int u = 0; u < locLengthDR[r][0]; u++){
    for (int v = 0; v < locLengthDR[r][1]; v++){
      for (int w = 0; w < locLengthDR[r][2]; w++){					
        // divided by DT, since we set DX*DY*DZ*DT = 1 in objective		
        dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][0] +=	
          -1 * signFac * (- EP0 * Ex_timedep_DR[r][maxTimeItr - timestep][u][v][w])
                  * ( Ex[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]] - Ex_adj_DR_old[r][u][v][w]) / DT;
        dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][1] +=	
          -1 * signFac * (- EP0 * Ey_timedep_DR[r][maxTimeItr - timestep][u][v][w])
                  * ( Ey[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]] - Ey_adj_DR_old[r][u][v][w]) / DT;
        dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][2] +=	
          -1 * signFac * (- EP0 * Ez_timedep_DR[r][maxTimeItr - timestep][u][v][w])
                  * ( Ez[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]] - Ez_adj_DR_old[r][u][v][w]) / DT;
        Ex_adj_DR_old[r][u][v][w] = Ex[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]];
        Ey_adj_DR_old[r][u][v][w] = Ey[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]];
        Ez_adj_DR_old[r][u][v][w] = Ez[u + shiftDR[r][0]][v + shiftDR[r][1]][w  + shiftDR[r][2]];
      }
    }
  }
  // Zero the gradients (locally) in the case there is a observation region located in DR subdomain
  for (int o = 0; o < NRMAXREGIONS; o++){
    int relDist[3];
    if (intrsOR[o]){
      for (int i = 0; i < 3; i++) relDist[i] = startIdxOR[o][i] - startIdxDR[r][i];
      for (int u = 0; u < locLengthDR[r][0]; u++){
        if (u >= relDist[0] && u < (relDist[0] + locLengthOR[o][0])){
          for (int v = 0; v < locLengthDR[r][1]; v++){
            if (v >= relDist[1] && v < (relDist[1] + locLengthOR[o][1])){
              for (int w = 0; w < locLengthDR[r][2]; w++){
                if (w >= relDist[2] && w < (relDist[2] + locLengthOR[o][2])){							
                  dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][0] = 0;
                  dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][1] = 0;
                  dfdx_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][2] = 0;
                }
              }
            }
          }
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------------
/*
Store electric field in r-th DR of the forward simulation, which will be used to update the gradients
*/
void StoreEfieldDR(int r){
  // Ghost points not needed
  for (int u = 0; u < locLengthDR[r][0]; u++){
    for (int v = 0; v < locLengthDR[r][1]; v++){
      for (int w = 0; w < locLengthDR[r][2]; w++){	
        // Saving E(t,x,y,z)
        //! -1 because, timestep >= 1, when called here (initial value = 1)
        Ex_timedep_DR[r][timestep - 1][u][v][w] = Ex[u + shiftDR[r][0]][v + shiftDR[r][1]][w + shiftDR[r][2]];
        Ey_timedep_DR[r][timestep - 1][u][v][w] = Ey[u + shiftDR[r][0]][v + shiftDR[r][1]][w + shiftDR[r][2]];
        Ez_timedep_DR[r][timestep - 1][u][v][w] = Ez[u + shiftDR[r][0]][v + shiftDR[r][1]][w + shiftDR[r][2]];
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------------
/*
Store electric field in OR(s) of the forward simulation, which will be used to build the adjoint source
and compute objective function (locally). Both subroutine will be done during one single function call
for the sake of effiency, since they share the same nested triple loop
*/
void StoreEfieldORAndComputeObj(int r){
  int matStartIdx[3] = {- NRSHIFTX + shiftOR[r][0],
                        - NRSHIFTY + shiftOR[r][1],
                        - NRSHIFTZ + shiftOR[r][2]};
  REAL objFunLoc_tmp = 0; // Local objective for current timestep and r
  for (int u = 0; u < locLengthOR[r][0] + numGPrOR[r][0]; u++){
    for (int v = 0; v < locLengthOR[r][1] + numGPrOR[r][1]; v++){
      for (int w = 0; w < locLengthOR[r][2] + numGPrOR[r][2]; w++){	
        // Saving E(t,x,y,z) -> ghost points needed!!!
        //! -1 because, timestep >= 1, when arriving here (initial value = 1)
        if (!preSim){
          Ex_timedep_OR[r][timestep - 1][u][v][w] = Ex[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]];
          Ey_timedep_OR[r][timestep - 1][u][v][w] = Ey[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]];
          Ez_timedep_OR[r][timestep - 1][u][v][w] = Ez[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]];
        }
        
        // Compute objective function 0.5 * integ {eps_0 * eps_inf * E^2} dxdydz, assuming "eps" is the background material
        // No ghost points included!
        if (u != locLengthOR[r][0] && v != locLengthOR[r][1] && w != locLengthOR[r][2]){
          objFunLoc_tmp = signFac * 0.5 * EP0 *
                          (  
                            EPRE[matindx[u + matStartIdx[0]][v + matStartIdx[1]][w + matStartIdx[2]]] 
                            * Ex[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] 
                            * Ex[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]]
                          + EPRE[matindy[u + matStartIdx[0]][v + matStartIdx[1]][w + matStartIdx[2]]]  
                            * Ey[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] 
                            * Ey[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]]
                          + EPRE[matindz[u + matStartIdx[0]][v + matStartIdx[1]][w + matStartIdx[2]]] 
                            * Ez[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]] 
                            * Ez[u + shiftOR[r][0]][v + shiftOR[r][1]][w + shiftOR[r][2]]
                          );
        }
      }
    }
  } 
  // Update objective
  objFunLoc[r] += objFunLoc_tmp;
  // --------------------------
  // Needed to evaluate change of fx
  objFunLoc_window += objFunLoc_tmp;	
}

// -------------------------------------------------------------------------------------------------------
/*
Evaluate change in fx to possibly abort forward sim
*/
void CheckObjChange(){
  if (OBJ_COMM_UNION != MPI_COMM_NULL){
    REAL fx_window = 0.;	
    MPI_Reduce(&objFunLoc_window, &fx_window, 1, MPIREAL, MPI_SUM, ROOT_OBJ_COMM_UNION, OBJ_COMM_UNION);
    if (abs(fx_window_old) > 1E-40){
      REAL frac = abs(1 - abs(fx_window / fx_window_old));
      if (RankThis == SendRankObjUnion) printf("Obj.change (t=%d): %.3e\n", timestep, frac);
      if (frac < obsDecayMinimum){
        printf("\n------------------------------------\n");
        printf("\nSimulation aborted after %d timesteps, since decay condition is satisfied (%.3e < %.3e)!\n", timestep, frac, obsDecayMinimum);
        printf("\n------------------------------------\n");		
        integrationTime = timestep;
        fx_window = 0; // Reset fx_window_old below
        abortSim = true; // Will abort sim
      }
    }
    else{
      if (RankThis == SendRankObjUnion) printf("Obj.change (t=%5d): not defined yet\n", timestep);
    }
    fx_window_old = fx_window;
    objFunLoc_window = 0;
  }
  // Communicate with all other processes 
  measureCount = 0;
  MPI_Bcast(&abortSim, 1, MPI_C_BOOL, SendRankObjUnion, MPI_COMM_WORLD);
}

// -------------------------------------------------------------------------------------------------------
/*
Sum up local objective values and send global objective value to optimizer ranks,
returns objective value 
*/
REAL SendObjToOptRanks(){
  // Call MPI allreduce to get the global objective value, stored by 1 process "SendRankObjUnion"
  REAL fx_tmp = 0;
  if (OBJ_COMM_UNION != MPI_COMM_NULL){
    REAL objFunLoc_union = 0.;
    for (int r = 0; r < NRMAXREGIONS; r++){
      if(intrsOR[r]) objFunLoc_union += objFunLoc[r];
    }					
    MPI_Reduce(&objFunLoc_union, &fx_tmp, 1, MPIREAL, MPI_SUM, ROOT_OBJ_COMM_UNION, OBJ_COMM_UNION);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  // Send the reduced objective value to all ranks sharing the Design Regions (all together!)
  if (OBJ_DR_COMM != MPI_COMM_NULL){
    MPI_Bcast(&fx_tmp, 1, MPIREAL, ROOT_OBJ_DR_COMM, OBJ_DR_COMM);
  }
  return fx_tmp;
}

// -------------------------------------------------------------------------------------------------------
/*
Zero objective function, save densities of current iteration and get densPhys_loc_arr which
will be used to update the E field based on a certain density <-> permittivity mapping
*/
void SetupForwardSim(){
  PetscErrorCode ierr;
  // Memory Allocation (For storing adjoint fields, ect...)
  AllocTopOptFields(); 
  // -------------------------------------------
  // Reset local objective function values
  for (int r = 0; r < NRMAXREGIONS; r++)	objFunLoc[r] = 0; 
  // -------------------------------------------
  if (!noIntrsDRAtAll){
    opt->fx = 0; // reset summed objective values
    bool saveVecs[3] = {true, true, true}; 	// MOVE PLOT OPTIONS TO JSON !!!	
    opt->SaveDensities(BASE_PATH_OUTPUT + RELATIVE_PATH_OPT_OUTPUT + "/Itr" + to_string(itr + 1) + "_", saveVecs);
    // ------------------------------------	
    // Get density array, which will be mapped to material later
    DMCreateLocalVector(opt->da_yee_lattice, &densPhys_loc);
    // Use local vector & create 4D (3 dims + 3 dofs) array from it (ghost cells accessible, which ARE necessary here!), 
    DMGlobalToLocalBegin(opt->da_yee_lattice, opt->densPhys, INSERT_VALUES, densPhys_loc); // delete that later
    DMGlobalToLocalEnd(opt->da_yee_lattice, opt->densPhys, INSERT_VALUES, densPhys_loc); // delete that later
    DMDAVecGetArrayDOF(opt->da_yee_lattice, densPhys_loc, &densPhys_loc_arr); // array indices in global dimensions!	 
    // ------------------------------------
  }
  adjoint_mode = false;
  itr++;
}

// -------------------------------------------------------------------------------------------------------
/*
Send global objective value to optimizer ranks, zero gradients, get gradient array which is
used in UpdateGradients()
*/				
void SetupAdjointSim(){
  // Sum up local objective values and send global objective value to optimizer ranks
  REAL fx_tmp = SendObjToOptRanks();

  // Zero gradients and constraints, get dfdx_arr
  if (OBJ_DR_COMM != MPI_COMM_NULL){
    // Assign fx_tmp to opt->fx for ranks sharing DRs
    if (!noIntrsDRAtAll){ // or PETSC_COMM_WORLD != MPI_COMM_NULL (==equivalent)
      // Compute volume constraint gx[0]
      opt->gx[0] = 0; // == NO CONSTRAINTS!
      VecSet(opt->dgdx[0], 0); 
      /* 
      // ------------------------------------
      // Volume constraint (used in mechanical problems usually)
      PetscInt ncelltot;
      VecGetSize(opt->densPhys, &ncelltot);
      VecSum(opt->densPhys, &(opt->gx[0]));
      opt->gx[0] = opt->gx[0] / (((PetscScalar)ncelltot)) - opt->volfrac;
      */
      // ------------------------------------
      // Assign objective value
      opt->fx = fx_tmp;

      // Scale objective and sens
      opt->fx = opt->fx * opt->fscale;

      // ------------------------------------
      // Reset dfdx and get array
      // -> will be computed in the following backward simulation
      VecSet(opt->dfdx, 0); 
      DMDAVecGetArrayDOF(opt->da_yee_lattice, opt->dfdx, &dfdx_arr); // array indices in global dimensions!
    }		
  }
  // ------------------------------------
  // Activate adjoint mode
  adjoint_mode = true;
}

// -------------------------------------------------------------------------------------------------------
/*
Update density, zero densPhys in observation regions again, print optimization variables
*/			
void UpdateDesign(){
  PetscErrorCode ierr;
  // ----------------------------------------
  // UPDATE DENSITY
  // ----------------------------------------
  DMDAVecRestoreArrayDOF(opt->da_yee_lattice, opt->dfdx, &dfdx_arr);
  VecScale(opt->dfdx, opt->fscale);

  // Filter sensitivities (chainrule)
  filter->Gradients(opt->dens, opt->densTilde, opt->dfdx, opt->m, opt->dgdx, opt->projectionFilter, 
                          opt->beta, opt->eta); 
  
  // Sets outer movelimits on design variables
  mma->SetOuterMovelimit(opt->Densmin, opt->Densmax, opt->movlim, opt->dens, opt->densmin, opt->densmax); 

  // Update design by MMA
  mma->Update(opt->dens, opt->dfdx, opt->gx, opt->dgdx, opt->densmin, opt->densmax);

  // Inf norm on the design change
  PetscScalar ch = mma->DesignChange(opt->dens, opt->dens_old);

  // Increase beta if needed
  PetscBool changeBeta = PETSC_FALSE;
  if (opt->projectionFilter) {
      changeBeta = filter->IncreaseBeta(&(opt->beta), opt->betaFinal, opt->gx[0], itr, ch);
  }
  
  // Filter design field
  filter->FilterProject(opt->dens, opt->densTilde, opt->densPhys, opt->projectionFilter, opt->beta, opt->eta);

  // Project densPhys in observation regions onto 0
  for (int r = 0; r < NRMAXREGIONS; r++){
    if (intrsDR[r]){
      for (int o = 0; o < NRMAXREGIONS; o++){
        if(intrsOR[o]){												
          PetscReal ****densPhys_arr;	
          DMDAVecGetArrayDOF(opt->da_yee_lattice, opt->densPhys, &densPhys_arr); // array indices in global dimensions!
          int petscStartIdx[3] = {startIdxDR[r][0] - globBoundsDR[r][0][0], 
                                  startIdxDR[r][1] - globBoundsDR[r][1][0], 
                                  startIdxDR[r][2] - globBoundsDR[r][2][0]};		
          int relDist[3] = {startIdxOR[o][0] - startIdxDR[r][0],
                            startIdxOR[o][1] - startIdxDR[r][1],
                            startIdxOR[o][2] - startIdxDR[r][2]};						
                            
          for (int u = 0; u < locLengthDR[r][0]; u++){
            if (u >= relDist[0] && u < (relDist[0] + locLengthOR[o][0])){
              for (int v = 0; v < locLengthDR[r][1]; v++){
                if (v >= relDist[1] && v < (relDist[1] + locLengthOR[o][1])){
                  for (int w = 0; w < locLengthDR[r][2]; w++){	
                    if (w >= relDist[2] && w < (relDist[2] + locLengthOR[o][2])){	
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][0] = 0;
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][1] = 0;
                      densPhys_arr[w + petscStartIdx[2]][v + petscStartIdx[1]][u + petscStartIdx[0]][2] = 0;
                    }
                  }
                }
              }
            }
          }
          DMDAVecRestoreArrayDOF(opt->da_yee_lattice, opt->densPhys, &densPhys_arr);
        }
      }
    }
  }

  // Discreteness measure
  PetscScalar mnd = filter->GetMND(opt->densPhys);

  // Print to screen
  PetscPrintf(PETSC_COMM_WORLD, "DR: %d: It.: %i,\nRuntime (s): %f,\nTime steps: %d,\nTrue fx: %24.21f,\nScaled fx: %f,\ngx[0]: %f,\nmnd.: %f,\nDens.change: %f,\nbeta: %f,\n",
                                  opt->identifierDR, itr, MPI_Wtime() - optStartTime, timestep - 1,  
                                  opt->fx / opt->fscale, opt->fx, opt->gx[0], mnd, ch, opt->beta);

  // Release TopOpt arrays
  DMDAVecRestoreArrayDOF(opt->da_yee_lattice, densPhys_loc, &densPhys_loc_arr);
  VecDestroy(&densPhys_loc);
  // Release field arrays
  ReleaseTopOptFields();
}

// -------------------------------------------------------------------------------------------------------
/*
Compute mean DFT_Ampl(freq) values in each OR respectively to evaluate broadband response
*/
void MeanDFT(int r){
  string MeanFileStr = BASE_PATH_OUTPUT + RELATIVE_PATH_OPT_OUTPUT + "/DFTMean_ObsReg" + to_string(r) + ".txt";
  FILE *string_DFT = fopen(MeanFileStr.c_str(), "w+");
  fprintf(string_DFT, "LAMBDA[nm]\t X\t Y\t Z\n");
  REAL volume = REAL((1 + TopOpt::boundsOR[r][0][1] - TopOpt::boundsOR[r][0][0])	
                    * (1 + TopOpt::boundsOR[r][1][1] - TopOpt::boundsOR[r][1][0])	
                    * (1 + TopOpt::boundsOR[r][2][1] - TopOpt::boundsOR[r][2][0]));
  for (int m = NRANFREQUENCIES - 1; m >= 0; m--){ 
    REAL meanLoc[3] = {0}; // x, y, z component
    string lambda_str = to_string(int((PI2 * C0 / OMEGA[m]) * REAL(1.e9)));
    int DFTStartIdx[3] = {- NRSHIFTX + shiftOR[r][0],
                          - NRSHIFTY + shiftOR[r][1],
                          - NRSHIFTZ + shiftOR[r][2]};

    for (int u = 0; u < locLengthOR[r][0]; u++){
      for (int v = 0; v < locLengthOR[r][1]; v++){
        for (int w = 0; w < locLengthOR[r][2]; w++){	
          meanLoc[0] += sqrt(  pow(ReTFEx[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2)
                            + pow(ImTFEx[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2));
          meanLoc[1] += sqrt(  pow(ReTFEy[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2)
                            + pow(ImTFEy[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2));
          meanLoc[2] += sqrt(  pow(ReTFEz[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2)
                              + pow(ImTFEz[m][u + DFTStartIdx[0]][v + DFTStartIdx[1]][w + DFTStartIdx[2]], 2));
        }
      }
    }
    REAL mean[3] = {0};		
    for (int i = 0; i < 3; i++){
      MPI_Reduce(&meanLoc[i], &mean[i], 1, MPIREAL, MPI_SUM, ROOT_OBJ_COMM_UNION, OBJ_COMM_UNION);
      mean[i] = mean[i] / volume;
    }
    if (RankThis == SendRankObjUnion){
      fprintf(string_DFT, "%e\t%e\t%e\t%e\n", (PI2 * C0 / OMEGA[m]) * REAL(1.e9), mean[0], mean[1], mean[2]);			
    }
  }
  fclose(string_DFT);			
}

// -------------------------------------------------------------------------------------------------------
/*
Finalize PETCs & delete optimization objects
*/
void FinalizeOpt(){
  delete opt;
  delete filter;
  delete mma;
  PetscFinalize(); 
}

// -------------------------------------------------------------------------------------------------------
/*
Allocate fields needed for the optimization procedure (computed during a forward/backward simulation)
*/
void AllocTopOptFields(){
	for (int r = 0; r < NRMAXREGIONS; r++){ 
		if (intrsDR[r]){
			// Ghost points not needed
			int spatialLenght[3] = {locLengthDR[r][0] - 1, 
															locLengthDR[r][1] - 1, 
															locLengthDR[r][2] - 1};
			MemAllocat4DArray(Ex_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			MemAllocat4DArray(Ey_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			MemAllocat4DArray(Ez_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			MemAllocat3DArray(Ex_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			MemAllocat3DArray(Ey_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			MemAllocat3DArray(Ez_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			// "E_adj_DR_NEW" would be the E field computed by backward simulation (but with adjoint source)
		}	
		if (intrsOR[r]){
			// Ghost points needed (because E field will updated via adjoint source)
				int spatialLenght[3] = {locLengthOR[r][0] - 1 + numGPrOR[r][0], 
																locLengthOR[r][1] - 1 + numGPrOR[r][1], 
																locLengthOR[r][2] - 1 + numGPrOR[r][2]};
			MemAllocat4DArray(Ex_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			MemAllocat4DArray(Ey_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			MemAllocat4DArray(Ez_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
		}
	}
}
// -------------------------------------------------------------------------------------------------------
/*
Release fields needed for the optimization procedure (computed during a forward/backward simulation)
*/
void ReleaseTopOptFields(){
	for (int r = 0; r < NRMAXREGIONS; r++){ 
		if (intrsDR[r]){
			// Ghost points not needed
			int spatialLenght[3] = {locLengthDR[r][0] - 1, 
															locLengthDR[r][1] - 1, 
															locLengthDR[r][2] - 1};
			if (Ex_timedep_DR[r] != NULL) MemRelease4DArray(Ex_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			if (Ey_timedep_DR[r] != NULL) MemRelease4DArray(Ey_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			if (Ez_timedep_DR[r] != NULL) MemRelease4DArray(Ez_timedep_DR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			if (Ex_adj_DR_old[r] != NULL) MemRelease3DArray(Ex_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			if (Ey_adj_DR_old[r] != NULL) MemRelease3DArray(Ey_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			if (Ez_adj_DR_old[r] != NULL) MemRelease3DArray(Ez_adj_DR_old[r], spatialLenght[0], spatialLenght[1], spatialLenght[2]);
			// "E_adj_DR_NEW" would be the E field computed by backward simulation (but with adjoint source)
		}
		if (intrsOR[r]){
			// Ghost points needed (because E field will updated via adjoint source)
			int spatialLenght[3] = {locLengthOR[r][0] - 1 + numGPrOR[r][0], 
															locLengthOR[r][1] - 1 + numGPrOR[r][1], 
															locLengthOR[r][2] - 1 + numGPrOR[r][2]};
			if (Ex_timedep_OR[r] != NULL) MemRelease4DArray(Ex_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			if (Ey_timedep_OR[r] != NULL) MemRelease4DArray(Ey_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
			if (Ez_timedep_OR[r] != NULL) MemRelease4DArray(Ez_timedep_OR[r], spatialLenght[0], spatialLenght[1], spatialLenght[2], NRMAXTIMESTEPS); 
		}
	}
}