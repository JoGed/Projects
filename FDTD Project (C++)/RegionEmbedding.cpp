#include <mpi.h>
#include <algorithm>
#include <string>
#include "RegionEmbedding.h"

using namespace std;
const char coordChar[3] = {'X', 'Y', 'Z'};

//---------------------------------------------------------------------------------------
bool Intersection(const int boundsA[3][2], const int boundsB[3][2]){
  bool AintrsB[3];
  for (int i = 0; i < 3; i++){
    AintrsB[i] = !(boundsA[i][1] < boundsB[i][0] || boundsA[i][0] > boundsB[i][1]);
  }
	return (AintrsB[0] && AintrsB[1] && AintrsB[2]);
};
//---------------------------------------------------------------------------------------
bool Inside(const int boundsA[3][2], const int boundsB[3][2]){
  bool AintrsB[3];
  for (int i = 0; i < 3; i++){
    AintrsB[i] = boundsA[i][0] > boundsB[i][0] && boundsA[i][1] < boundsB[i][1];
  }
  return (AintrsB[0] && AintrsB[1] && AintrsB[2]);
};
//---------------------------------------------------------------------------------------
void SumUpRanks(const int intrs, int* ranksNum, MPI_Comm comm){
  MPI_Allreduce(&intrs, ranksNum, 1, MPI_INT, MPI_SUM, comm); 
};
//---------------------------------------------------------------------------------------
void GetShifts(const int boundsA[3][2], const int boundsB[3][2], int* shift){
  for(int i = 0; i < 3; i++){
    shift[i] = max(0, boundsA[i][0] - boundsB[i][0]);
  } 
};
//---------------------------------------------------------------------------------------
void GetStartIndices(const int boundsA[3][2], const int boundsB[3][2], int* startIdx){
  for(int i = 0; i < 3; i++){
    startIdx[i] = boundsB[i][0] + max(0, boundsA[i][0] - boundsB[i][0]);
  } 
};
//---------------------------------------------------------------------------------------
void GetLocalLength(const int boundsA[3][2], const int boundsB[3][2], const int lengthChunkFDTD[3], 
                    int* localLength){
  for(int i = 0; i < 3; i++){
    localLength[i] = lengthChunkFDTD[i] - max(0, boundsA[i][0] - boundsB[i][0]) - max(0, boundsB[i][1] - boundsA[i][1]);
  } 
};
//---------------------------------------------------------------------------------------
void GetRanks(const MPI_Comm commParent, const int intrs, int* processes){

  int rankFDTD, globalProcSize;                                
  MPI_Comm_rank(commParent, &rankFDTD);
	MPI_Comm_size(commParent, &globalProcSize);                               
  int* intrsAll = new int[globalProcSize];
	MPI_Allgather(&intrs, 1, MPI_INT, intrsAll, 1, MPI_INT, commParent);

	int r = 0;
  for (int i = 0; i < globalProcSize; i++){		
      if(intrsAll[i]){	
        processes[r] = i;
        r++;	
      }
  } 
	delete [] intrsAll;    
};
//---------------------------------------------------------------------------------------
void CountRanksPerAxis(const MPI_Comm comm, const int* startIdx, int* numProc){

      fill_n(numProc, 3, 1); // Initialize to 1 (assuming 1 process per axis)
      int size;
	    MPI_Comm_size(comm, &size);
      int* coordsArrAll = new int[3 * size];
      MPI_Allgather(&startIdx[0], 3, MPI_INT, coordsArrAll, 3, MPI_INT, comm);
      int coordsCell[3] = {coordsArrAll[0], coordsArrAll[1], coordsArrAll[2]};

      // X is the fastest index for the arrangement of subdomains!
      auto count = [&](const int scale, const int coord) -> int{
        int max_it =  size / scale;
        // Start to count at next neighbour
        for (int d = 1; d < max_it + 1; d++){
          bool isSameCoord = coordsCell[coord] == coordsArrAll[3 * d * scale + coord];
          if (d == max_it){ 
            break; 
          }
          // <- Avoid memory corruption
          if (isSameCoord){
            break;  
          }
          else {
            numProc[coord] += 1;
          }
        }    
      };

      if (size > 1){
          count(1, 0);
          count(numProc[0], 1);
          count(numProc[0] * numProc[1], 2);
      }
};
//---------------------------------------------------------------------------------------
void GetCellArrayPerAxis(const MPI_Comm comm, const int* localLength, const int ranksNum, 
                                const int* numProc, int* cellsArr[3]){
  
  // Arrays containing number of cells of each process in X
  int* cellsArrAll_x = new int[ranksNum]();
  int* cellsArrAll_y = new int[ranksNum]();
  int* cellsArrAll_z = new int[ranksNum]();

  MPI_Allgather(&localLength[0], 1, MPI_INT, cellsArrAll_x, 1, MPI_INT, comm);
  MPI_Allgather(&localLength[1], 1, MPI_INT, cellsArrAll_y, 1, MPI_INT, comm);
  MPI_Allgather(&localLength[2], 1, MPI_INT, cellsArrAll_z, 1, MPI_INT, comm); 

  for(int i = 0; i < numProc[0]; i++){ *(cellsArr[0] + i) = cellsArrAll_x[i]; }
  for(int j = 0; j < numProc[1]; j++){ *(cellsArr[1] + j) = cellsArrAll_y[j * numProc[0]]; }
  for(int k = 0; k < numProc[2]; k++){ *(cellsArr[2] + k) = cellsArrAll_z[k * numProc[0] * numProc[1]]; }

  delete [] cellsArrAll_x;
  delete [] cellsArrAll_y;
  delete [] cellsArrAll_z; 
}
//---------------------------------------------------------------------------------------
void PrintDecomposition(const int commSize, const int globalProcSize, const int* processes, 
                        const int* numProc, const int* cellsArr[3], const string regionName){
  printf("\n\n");
  printf("------------------------------------\n");
  printf(" %s\n", regionName);
  printf("------------------------------------\n");
  printf("Ranks sharing '%s':\n{", regionName); 
  for (int i = 0; i < commSize; i++){ 
    printf("%i, ", processes[i]); 
  }
  printf("}\n--------\n"); 
  printf("--> %i / %i ranks \n", commSize, globalProcSize); 
  printf("------------------------------------\n");
  for(int i = 0; i < 3; i++)
  {
    printf("Num. of ranks sharing '%s' in %c: %i \n", regionName, coordChar[i], numProc[i]); 				
  }
  if (cellsArr != NULL){
    for (int i = 0; i < 3; i++)
    {
      printf("------------------------------------\n");
      printf("Cells of '%s' in %c:\n{", regionName, coordChar[i]);
      for (int p = 0; p < numProc[i]; p++)
      {
        printf("%i, ", *(cellsArr[i] + p));
      }
      printf("}\n");
    } 
  }
};
