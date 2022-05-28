/**
 * @file io.cpp
 * @author FDTD-Lab
 * @brief Functions to write out binary files and check output paths
 * @version 0.1
 * @date 2022-04-12
 * 
 * @copyright Copyright (c) 2022
 */
#include <string>
#include "mpi.h"
#include <typeinfo>
#include "definitions.h"
#include <type_traits>
//#include <filesystem>

using namespace std;
//namespace fs = std::filesystem;


/**
 * @brief Checks to see if path exists and is directory
 * If doesn't exist, will create the path
 * 
 * @param dirPath [in] path string to directory you wish to check
 * @return true if error eccured
 * @return false if there were no issues
 */
// bool EnsureDirectory(std::string &dirPath)
// {
// 	fs::path tmpPath{dirPath};

// 	if ( fs::exists(tmpPath) )
// 	{
// 		if ( !fs::is_directory(dirPath) )
// 		{
// 			printf("\nERROR: Path exists but is not a directory\n\n");
// 			return true;
// 		}
// 		else
// 			return false;
// 	}
// 	return !fs::create_directory(dirPath);
// }

// -------------------------------
// MEMORY ALLOCATION 2D 
// -------------------------------
template <typename T>
void MemAllocat2DArray(T **&p, int const &NX, int const &NY)
{
	int nx = NX + 1, ny = NY + 1;
	int i, j, shift;

	p = new T *[nx];

	p[0] = new T[nx * ny];

	for (i = 0; i < nx; ++i)
	{
		shift = i * ny;
		if (i != 0)
		{
			p[i] = p[0] + shift;
		}
	}

	for (i = 0; i < nx; ++i)
	{
		for (j = 0; j < ny; ++j)
		{
			p[i][j] = static_cast<T>(0);
		}
	}

	return;
}
template void MemAllocat2DArray<int>(int **&p, int const &NX, int const &NY);
template void MemAllocat2DArray<REAL>(REAL **&p, int const &NX, int const &NY);
template void MemAllocat2DArray<CMPLX>(CMPLX **&p, int const &NX, int const &NY);

// -------------------------------
// MEMORY RELEASE 2D
// -------------------------------
template <typename T>
void MemRelease2DArray(T **&p, int const &NX, int const &NY)
{
	delete[] p[0];
	delete[] p;

	return;
}

template void MemRelease2DArray<int>(int **&p, int const &NX, int const &NY);
template void MemRelease2DArray<REAL>(REAL **&p, int const &NX, int const &NY);
template void MemRelease2DArray<CMPLX>(CMPLX **&p, int const &NX, int const &NY);

// -------------------------------
// MEMORY ALLOCATION 3D
// -------------------------------
template <typename T>
void MemAllocat3DArray(T ***&p, int const &NX, int const &NY, int const &NZ)
{
	int nx = NX + 1, ny = NY + 1, nz = NZ + 1;
	int i, j, k, shift1, shift2;
	int dim2 = ny * nz, dim3 = nx * ny * nz;

	p = new T **[nx];

	for (i = 0; i < nx; ++i)
		p[i] = new T *[ny];

	p[0][0] = new T[dim3];

	for (i = 0; i < nx; ++i)
	{
		shift1 = i * dim2;

		for (j = 0; j < ny; ++j)
		{
			shift2 = j * nz;

			if ((i != 0) || (j != 0))
				p[i][j] = p[0][0] + shift1 + shift2;
		}
	}

	for (i = 0; i < nx; ++i)
	{
		for (j = 0; j < ny; ++j)
		{
			for (k = 0; k < nz; ++k)
			{
				p[i][j][k] = static_cast<T>(0);
			}
		}
	}

	return;
}
template void MemAllocat3DArray<int>(int ***&p, int const &NX, int const &NY, int const &NZ);
template void MemAllocat3DArray<REAL>(REAL ***&p, int const &NX, int const &NY, int const &NZ);
template void MemAllocat3DArray<CMPLX>(CMPLX ***&p, int const &NX, int const &NY, int const &NZ);


// -------------------------------
// MEMORY RELEASE 3D 
// -------------------------------
template <typename T>
void MemRelease3DArray(T ***&p, int const &NX, int const &NY, int const &NZ)
{
	int i, nx = NX + 1;

	delete[] p[0][0];

	for (i = 0; i < nx; ++i)
		delete[] p[i];

	delete[] p;

	return;
}

template void MemRelease3DArray<int>(int ***&p, int const &NX, int const &NY, int const &NZ);
template void MemRelease3DArray<REAL>(REAL ***&p, int const &NX, int const &NY, int const &NZ);
template void MemRelease3DArray<CMPLX>(CMPLX ***&p, int const &NX, int const &NY, int const &NZ);

// -------------------------------
// MEMORY ALLOCATION 4D 
// -------------------------------
template <typename T>
void MemAllocat4DArray(T ****&arr, int const &NX, int const &NY, int const &NZ, int const &t_max)
{	
	int nx = NX + 1, ny = NY + 1, nz = NZ + 1;
	int i, j, k, t, shift1, shift2, shift3;

	arr = new T ***[t_max];
	
	for (t = 0; t < t_max; t++)
	{
		arr[t] = new T **[nx];
		for (i = 0; i < nx; i++)
		{
			arr[t][i] = new T *[ny];
		}
		
	}
	
	arr[0][0][0] = new T[t_max * nx * ny * nz];

	for (t = 0; t < t_max; t++)
	{
		shift1 = t * nx * ny * nz;

		for (i = 0; i < nx; i++)
		{
			shift2 = i *  ny * nz;
			
			for (j = 0; j < ny; j++)
			{
				shift3 = j * nz;

				if ((t != 0) || (i != 0) || (j != 0))
					arr[t][i][j] = arr[0][0][0] + shift1 + shift2 + shift3;
			}
		}
	}
	
	for (t = 0; t < t_max; t++)
	{
		for (i = 0; i < nx; i++)
		{
			for (j = 0; j < ny; j++)
			{
				for (k = 0; k < nz; k++)
				{ 	
					arr[t][i][j][k] = static_cast<T>(0);	
				}
			}
		}
	}
}

template void MemAllocat4DArray<int>(int ****&p, int const &NX, int const &NY, int const &NZ, int const &time);
template void MemAllocat4DArray<REAL>(REAL ****&p, int const &NX, int const &NY, int const &NZ, int const &time);
template void MemAllocat4DArray<CMPLX>(CMPLX ****&p, int const &NX, int const &NY, int const &NZ, int const &time);

template <typename T>
void MemRelease4DArray(T ****&arr, int const &NX, int const &NY, int const &NZ, int const &t_max)
{
	int nx = NX + 1, ny = NY + 1, nz = NZ + 1; // +1 extra point == right ghost point
	int t, i, j, k;
	
	delete[] arr[0][0][0];

	for (t = 0; t < t_max; ++t)
	{
		for (i = 0; i < nx; ++i)
		{
			delete[] arr[t][i];
		}
		delete[] arr[t];	
	}
	delete[] arr;

	return;
}
template void MemRelease4DArray<int>(int ****&p, int const &NX, int const &NY, int const &NZ, int const &time);
template void MemRelease4DArray<REAL>(REAL ****&p, int const &NX, int const &NY, int const &NZ, int const &time);
template void MemRelease4DArray<CMPLX>(CMPLX ****&p, int const &NX, int const &NY, int const &NZ, int const &time);

// -------------------------------
// OUPUT BINARY FILE FOR 2D PLOT
// -------------------------------
/**
 * @brief Stores the values of all 3 components of a given 3D field for a specific plane (2D) into a binary file.
 * The first 3 stored values are integers specifying the number of points of the
 * simulation domain (without PMLs) in X, Y, Z direction respectively. The 2D array are
 * stored in ij order, where j is the inner/fastest loop.
 * 
 * @tparam Location to put pointer of type T to the 3D array. Supports Integer, Float and Double type
 * @param filename Name of the binary file
 * @param arr_x x component of the 3D array
 * @param arr_x y component of the 3D array
 * @param arr_x z component of the 3D array
 * @param outputCommPlotPlane Communicator holding the group of processes whose domain chunks intersect the 2D plane and are non-pure CPML domains
 * @param normalToPlaneIdx Index (or lenght in cells) of vector normal to the plane
 * @param startindices Global (regarding TFSF region without CPML) smallest indices whose corresponding cells being part of the domain
 * @param nrInnerThis Local number of cells lying inside the physical domain (without CPMLs) 
 * @param nrShift Local Indices shift from total simulation domain (including CPMLs) to "physical" domain (not including CPMLs)
 * @param idx_coord_A Indices of positions of a list {"X", "Y", "Z"}, where idx_coord_A must be smaller than idx_coord_B
 * @param idx_coord_B Indices of positions of a list {"X", "Y", "Z"}, where idx_coord_B must be bigger than idx_coord_A
 * @param use_global_coords Set "true" if arrays are indexed globally (regarding domain without CPMLs). 
 * Set "false" if arrays are indexed locally (ALL starting indices == [0][0][0])
 * @param ijk_order Set 'true', if 3D arrays are ordered as arr[i][j][k], 
 * otherwise set 'false' it they are ordered as arr[k][j][i]
 * @param with_PMLs Set 'true', if 3D arrays (in a global view) are defined on the physical domain + CPMLs, otherwise choose 'false'
 */
template <typename T>
void SaveAsBinary2D(string const &filename, T ***&arr_x, T ***&arr_y, T ***&arr_z, MPI_Comm const &outputCommPlotPlane, int const &normalToPlaneIdx, 
int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, int const &idx_coord_A, int const &idx_coord_B, 
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs)
{	
	int localarraysizePlane;
	int idxShift[3] = {}; // Initialize to 0
	T **arr_x_Plane;
	T **arr_y_Plane;
	T **arr_z_Plane;
	MPI_File mpiFile;
	MPI_Datatype Vector_dim;
	MPI_Datatype VectorPlane_pre;
	MPI_Datatype VectorPlane;
	MPI_Datatype MPITYPE;
	int MPITYPE_size;
	MPI_Status status;
	//-------------------------
	int dim = 3; // x,y,z component
	localarraysizePlane = nrInnerThis[idx_coord_A] * nrInnerThis[idx_coord_B];
	//-------------------------
	// Adapt indice shift for proper array access
	if (use_global_coords)
	{
		idxShift[0] = startindices[0];
		idxShift[1] = startindices[1];
		idxShift[2] = startindices[2];
	}
	if (with_PMLs)
	{
		idxShift[0] += nrShift[0];
		idxShift[1] += nrShift[1];
		idxShift[2] += nrShift[2];
	}
	//-------------------------
	// Check T type
	if (is_same<T, int>::value)
	{ 
		MPITYPE = MPI_INT; 
	}
	else if (is_same<T, float>::value)
	{ 
		MPITYPE = MPI_FLOAT; 
	}
	else if (is_same<T, double>::value)
	{ 
		MPITYPE = MPI_DOUBLE; 
	}
	else if (is_same<T, complex<double>>::value)
	{ 
		MPITYPE = MPI_DOUBLE_COMPLEX; 
	}
	else{
		int printRank;
		MPI_Comm_rank(outputCommPlotPlane, &printRank);
		if(printRank == 0){
			printf("------------------------------------\n");
			printf("Warning: Datatype of given arrays in %s are not supported!\n", __func__);
			printf("------------------------------------\n");
		}
	}	
	MPI_Type_size(MPITYPE,  &MPITYPE_size);		
	//-------------------------
	// Allocate memory
	MemAllocat2DArray(arr_x_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);
	MemAllocat2DArray(arr_y_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);
	MemAllocat2DArray(arr_z_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);
	//-------------------------
	// initialization of 2D arrays (plane)
	for (int i = 0; i < nrInnerThis[idx_coord_A]; i++)
	{
		for (int j = 0; j < nrInnerThis[idx_coord_B]; j++)
		{
			if ((idx_coord_A == 0) && (idx_coord_B == 1)) // XY
			{
				if (ijk_order)
				{	
					arr_x_Plane[i][j] = arr_x[i + idxShift[0]][j + idxShift[1]][normalToPlaneIdx + idxShift[2]];
					arr_y_Plane[i][j] = arr_y[i + idxShift[0]][j + idxShift[1]][normalToPlaneIdx + idxShift[2]];
					arr_z_Plane[i][j] = arr_z[i + idxShift[0]][j + idxShift[1]][normalToPlaneIdx + idxShift[2]];
					
				}
				else
				{
					arr_x_Plane[i][j] = arr_x[normalToPlaneIdx + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
					arr_y_Plane[i][j] = arr_y[normalToPlaneIdx + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
					arr_z_Plane[i][j] = arr_z[normalToPlaneIdx + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
				}
			}
			else if ((idx_coord_A == 0) && (idx_coord_B == 2)) // XZ
			{
				if (ijk_order)
				{
					arr_x_Plane[i][j] = arr_x[i + idxShift[0]][normalToPlaneIdx + idxShift[1]][j + idxShift[2]];
					arr_y_Plane[i][j] = arr_y[i + idxShift[0]][normalToPlaneIdx + idxShift[1]][j + idxShift[2]];
					arr_z_Plane[i][j] = arr_z[i + idxShift[0]][normalToPlaneIdx + idxShift[1]][j + idxShift[2]];
				}
				else
				{
					arr_x_Plane[i][j] = arr_x[j + idxShift[2]][normalToPlaneIdx + idxShift[1]][i + idxShift[0]];
					arr_y_Plane[i][j] = arr_y[j + idxShift[2]][normalToPlaneIdx + idxShift[1]][i + idxShift[0]];
					arr_z_Plane[i][j] = arr_z[j + idxShift[2]][normalToPlaneIdx + idxShift[1]][i + idxShift[0]];
				}
			}
			else if ((idx_coord_A == 1) && (idx_coord_B == 2)) // YZ
			{
				if (ijk_order)
				{
					arr_x_Plane[i][j] = arr_x[normalToPlaneIdx + idxShift[0]][i + idxShift[1]][j + idxShift[2]];
					arr_y_Plane[i][j] = arr_y[normalToPlaneIdx + idxShift[0]][i + idxShift[1]][j + idxShift[2]];
					arr_z_Plane[i][j] = arr_z[normalToPlaneIdx + idxShift[0]][i + idxShift[1]][j + idxShift[2]];
				}
				else
				{
					arr_x_Plane[i][j] = arr_x[j + idxShift[2]][i + idxShift[1]][normalToPlaneIdx + idxShift[0]];
					arr_y_Plane[i][j] = arr_y[j + idxShift[2]][i + idxShift[1]][normalToPlaneIdx + idxShift[0]];
					arr_z_Plane[i][j] = arr_z[j + idxShift[2]][i + idxShift[1]][normalToPlaneIdx + idxShift[0]];
				}
			}
		}
	}
	//-------------------------
	// Writing 2D arrays to binary file
	MPI_Type_create_hvector(nrInnerThis[idx_coord_B], 1, dim * MPITYPE_size, MPITYPE, &VectorPlane_pre);
	MPI_Type_commit(&VectorPlane_pre);
	MPI_Type_create_hvector(nrInnerThis[idx_coord_A], 1, dim * MPITYPE_size * globalInnerSize[idx_coord_B], VectorPlane_pre, &VectorPlane);
	MPI_Type_commit(&VectorPlane);
	MPI_Type_create_hvector(3, 1, MPI_INT_size, MPI_INT, &Vector_dim);
	MPI_Type_commit(&Vector_dim);
	//-------------------------
	int disp_origin = 0;
	int disp;
	MPI_File_open(outputCommPlotPlane, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile);
	MPI_File_set_view(mpiFile, disp_origin, MPI_INT, Vector_dim, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &globalInnerSize[0], 3, MPI_INT, &status); //  inner_dim <-> globalInnerSize[0], why does just "globalInnerSize" not work
	disp = 3 * MPI_INT_size + dim * MPITYPE_size * (startindices[idx_coord_B] + startindices[idx_coord_A] * globalInnerSize[idx_coord_B]);
	MPI_File_set_view(mpiFile, disp, MPITYPE, VectorPlane, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_x_Plane[0][0], localarraysizePlane, MPITYPE, &status);
	MPI_File_set_view(mpiFile, disp + 1 * MPITYPE_size, MPITYPE, VectorPlane, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_y_Plane[0][0], localarraysizePlane, MPITYPE, &status);
	MPI_File_set_view(mpiFile, disp + 2 * MPITYPE_size, MPITYPE, VectorPlane, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_z_Plane[0][0], localarraysizePlane, MPITYPE, &status);
	MPI_File_close(&mpiFile);
	//-------------------------
	// Release memory
	MemRelease2DArray(arr_x_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);
	MemRelease2DArray(arr_y_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);
	MemRelease2DArray(arr_z_Plane, nrInnerThis[idx_coord_A] - 1, nrInnerThis[idx_coord_B] - 1);	
	//-------------------------
	MPI_Barrier(outputCommPlotPlane);
}

template void SaveAsBinary2D<int>(string const &filename, int ***&arr_x, int ***&arr_y, int ***&arr_z, MPI_Comm const &outputCommPlotPlane, 
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, 
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void SaveAsBinary2D<REAL>(string const &filename, REAL ***&arr_x, REAL ***&arr_y, REAL ***&arr_z, MPI_Comm const &outputCommPlotPlane, 
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,  
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);
template void SaveAsBinary2D<CMPLX>(string const &filename, CMPLX ***&arr_x, CMPLX ***&arr_y, CMPLX ***&arr_z, MPI_Comm const &outputCommPlotPlane, 
int const &normalToPlaneIdx, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,  
int const &idx_coord_A, int const &idx_coord_B, bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);


// -------------------------------
// OUPUT BINARY FILE FOR 3D PLOT
// -------------------------------
/**
 * @brief Stores the values of all 3 components of a given 3D field for a specific plane into a binary file.
 * The first 3 stored values are integers specifying the number of points of the
 * simulation domain (without PMLs) in X, Y, Z direction respectively. The 3D array are
 * stored in ijk order, where k is the inner/fastest loop.
 * 
 * @tparam Location to put pointer of type T to the 3D array. Supports Integer, Float and Double type
 * @param filename Name of the binary file
 * @param arr_x x component of the 3D array
 * @param arr_x y component of the 3D array
 * @param arr_x z component of the 3D array
 * @param outputCommPlot Communicator holding the group of processes whose domain chunks are non-pure CPML domains
 * @param startindices Global (regarding TFSF region without CPML) smallest indices whose corresponding cells being part of the domain
 * @param nrInnerThis Local number of cells lying inside the physical domain (without CPMLs) 
 * @param nrShift Local Indices shift from total simulation domain (including CPMLs) to "physical" domain (not including CPMLs)
 * @param use_global_coords Set "true" if arrays are indexed globally (regarding domain without CPMLs). 
 * Set "false" if arrays are indexed locally (ALL starting indices == [0][0][0])
 * @param ijk_order Set 'true', if 3D arrays are ordered as arr[i][j][k], 
 * otherwise set 'false' it they are ordered as arr[k][j][i]
 * @param with_PMLs Set 'true', if 3D arrays (in a global view) are defined on the physical domain + CPMLs, otherwise choose 'false'
 */
template <typename T>
void SaveAsBinary3D(string const &filename, T ***&arr_x, T ***&arr_y, T ***&arr_z,  MPI_Comm const &outputCommPlot,
int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift,
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs)
{	
	int localarraysizeVolume;
	int idxShift[3] = {0, 0, 0}; // Initialize to 0
	T ***arr_x_Vol;
	T ***arr_y_Vol;
	T ***arr_z_Vol;
	MPI_File mpiFile;
	MPI_Datatype Vector_dim;
	MPI_Datatype VectorPlane_pre;
	MPI_Datatype VectorPlane;
	MPI_Datatype VectorVolume;
	MPI_Datatype MPITYPE;
	int MPITYPE_size;
	MPI_Status status;
	//-------------------------
	int dim = 3; // x,y,z component
	localarraysizeVolume = nrInnerThis[0] * nrInnerThis[1] * nrInnerThis[2];
	//-------------------------
	// Adapt indice shift for proper array access
	if (use_global_coords)
	{
		idxShift[0] = startindices[0];
		idxShift[1] = startindices[1];
		idxShift[2] = startindices[2];
	}
	if (with_PMLs) 
	{
		idxShift[0] += nrShift[0];
		idxShift[1] += nrShift[1];
		idxShift[2] += nrShift[2];
	}
	//-------------------------
	// Check T type
	if (is_same<T, int>::value)
	{ 
		MPITYPE = MPI_INT; 
	}
	else if (is_same<T, float>::value)
	{ 
		MPITYPE = MPI_FLOAT; 
	}
	else if (is_same<T, double>::value)
	{ 
		MPITYPE = MPI_DOUBLE; 
	}
	else if (is_same<T, complex<double>>::value)
	{ 
		MPITYPE = MPI_DOUBLE_COMPLEX; 
	}
	else{
		int printRank;
		MPI_Comm_rank(outputCommPlot, &printRank);
		if(printRank == 0){
			printf("------------------------------------\n");
			printf("Warning: Datatype of given arrays in %s are not supported!\n", __func__);
			printf("------------------------------------\n");
		}
	}
	MPI_Type_size(MPITYPE,  &MPITYPE_size);
	//-------------------------
	// Allocate memory 
	MemAllocat3DArray(arr_x_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1); 
	MemAllocat3DArray(arr_y_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1);
	MemAllocat3DArray(arr_z_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1);
	//-------------------------

	for (int i = 0; i < nrInnerThis[0]; i++)
	{
		for (int j = 0; j < nrInnerThis[1]; j++)
		{   
			for (int k = 0; k < nrInnerThis[2]; k++)
			{   
				if (ijk_order)
				{	
					arr_x_Vol[i][j][k] = arr_x[i + idxShift[0]][j + idxShift[1]][k + idxShift[2]];
					arr_y_Vol[i][j][k] = arr_y[i + idxShift[0]][j + idxShift[1]][k + idxShift[2]];
					arr_z_Vol[i][j][k] = arr_z[i + idxShift[0]][j + idxShift[1]][k + idxShift[2]];	
				}	
				else
				{
					arr_x_Vol[i][j][k] = arr_x[k + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
					arr_y_Vol[i][j][k] = arr_y[k + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
					arr_z_Vol[i][j][k] = arr_z[k + idxShift[2]][j + idxShift[1]][i + idxShift[0]];
				}
			}
		}
	}	
	//-------------------------
	// Writing 2D arrays to binary file
	MPI_Type_create_hvector(nrInnerThis[2], 1, dim * MPITYPE_size, MPITYPE, &VectorPlane_pre);
	MPI_Type_commit(&VectorPlane_pre);
	MPI_Type_create_hvector(nrInnerThis[1], 1, dim * MPITYPE_size * globalInnerSize[2], VectorPlane_pre, &VectorPlane);
	MPI_Type_commit(&VectorPlane);
	MPI_Type_create_hvector(nrInnerThis[0], 1, dim * MPITYPE_size * globalInnerSize[2] * globalInnerSize[1], VectorPlane, &VectorVolume);
	MPI_Type_commit(&VectorVolume);
	MPI_Type_create_hvector(3, 1, MPI_INT_size, MPI_INT, &Vector_dim); // used to save info about number of cells in x, y, z of domain without PMLs
	MPI_Type_commit(&Vector_dim);
	//-------------------------
	int disp_origin = 0;
	int disp;
	MPI_File_open(outputCommPlot, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpiFile);
	MPI_File_set_view(mpiFile, disp_origin, MPI_INT, Vector_dim, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &globalInnerSize[0], 3, MPI_INT, &status);
	disp = 3 * MPI_INT_size + dim * MPITYPE_size * (startindices[2] + startindices[1] * globalInnerSize[2] + startindices[0] * globalInnerSize[1] * globalInnerSize[2]);
	MPI_File_set_view(mpiFile, disp + 0 * MPITYPE_size, MPITYPE, VectorVolume, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_x_Vol[0][0][0], localarraysizeVolume, MPITYPE, &status);
	MPI_File_set_view(mpiFile, disp + 1 * MPITYPE_size, MPITYPE, VectorVolume, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_y_Vol[0][0][0], localarraysizeVolume, MPITYPE, &status);
	MPI_File_set_view(mpiFile, disp + 2 * MPITYPE_size, MPITYPE, VectorVolume, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpiFile, &arr_z_Vol[0][0][0], localarraysizeVolume, MPITYPE, &status);
	MPI_File_close(&mpiFile);
	// Release memory
	MemRelease3DArray(arr_x_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1); 
	MemRelease3DArray(arr_y_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1); 
	MemRelease3DArray(arr_z_Vol, nrInnerThis[0] - 1, nrInnerThis[1] - 1, nrInnerThis[2] - 1); 	
	//-------------------------
	MPI_Barrier(outputCommPlot);
}

template void SaveAsBinary3D<int>(string const &filename, int ***&arr_x, int ***&arr_y, int ***&arr_z, 
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, 
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void SaveAsBinary3D<REAL>(string const &filename, REAL ***&arr_x, REAL ***&arr_y, REAL ***&arr_z, 
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, 
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);

template void SaveAsBinary3D<CMPLX>(string const &filename, CMPLX ***&arr_x, CMPLX ***&arr_y, CMPLX ***&arr_z, 
MPI_Comm const &outputComm, int *const &startindices, int *const &globalInnerSize, int *const &nrInnerThis, int *const &nrShift, 
bool const &use_global_coords, bool const &ijk_order, bool const &with_PMLs);