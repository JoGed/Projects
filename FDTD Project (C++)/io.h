#ifndef IO_H
#define IO_H

#include <string>
#include "mpi.h"

using namespace std;

template <typename T> void MemAllocat2DArray(T **&, int const &, int const &);
template <typename T> void MemRelease2DArray(T **&, int const &, int const &);
template <typename T> void SaveAsBinary2D(string const &, T ***&, T ***&, T ***&, MPI_Comm const &, int const &, 
int* const &, int* const &, int* const &, int* const &, int const &, int const &, bool const &, bool const &, bool const&);

template <typename T> void MemAllocat3DArray(T ***&, int const &, int const &, int const &);
template <typename T> void MemRelease3DArray(T ***&, int const &, int const &, int const &);
template <typename T> void SaveAsBinary3D(string const &, T ***&, T ***&, T ***&, MPI_Comm const &, int* const &, int* const &, 
int* const &, int* const &, bool const &, bool const &, bool const&);

template <typename T> void MemAllocat4DArray(T ****&, int const &, int const &, int const &, int const &);
template <typename T> void MemRelease4DArray(T ****&, int const &, int const &, int const &, int const &);

bool EnsureDirectory(std::string &dirPath); 

#endif  