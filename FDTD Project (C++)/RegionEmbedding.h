#ifndef EMBEDDINGOPT_H
#define EMBEDDINGOPT_H
#include <mpi.h>
#include <string>

bool Intersection(const int boundsA[3][2], const int boundsB[3][2]);

bool Inside(const int boundsA[3][2], const int boundsB[3][2]);

void SumUpRanks(const int intrs, int* ranksNum, MPI_Comm comm);

void GetShifts(const int boundsA[3][2], const int boundsB[3][2], int* shift);

void GetStartIndices(const int boundsA[3][2], const int boundsB[3][2], int* startIdx);

void GetLocalLength(const int boundsA[3][2], const int boundsB[3][2], const int lengthChunkFDTD[3], int* localLength);

void GetRanks(const MPI_Comm commParent, const int intrs, int* processes);

void CountRanksPerAxis(const MPI_Comm comm, const int* startIdx, int* numProc);

void GetCellArrayPerAxis(const MPI_Comm comm, const int* localLength, const int ranksNum, 
                                const int* numProc, int* cellsArr[3]);

void PrintDecomposition(const int commSize, const int globalProcSize, const int* processes, 
                               const int* numProc, const int* cellsArr[3], const std::string regionName);

#endif