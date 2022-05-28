#include "Filter.h"

/* -----------------------------------------------------------------------------
Adapted Code for electromagnetic optimization problems, modified by Johannes Gedeon
Original Code:
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
Copyright (C) 2013-2014,

This Filter implementation is licensed under Version 2.1 of the GNU
Lesser General Public License.

This MMA implementation is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This Module is distributed in the hope that it will be useful,implementation
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
-------------------------------------------------------------------------- */

Filter::Filter(DM da_yee_lattice, Vec x, PetscInt filterT, PetscScalar radiusVec[6], PetscScalar weightVec[6]) {
    // Set all pointers to NULL
    H       = NULL;
    Hs      = NULL;
  
    // Create radius and weight matrix
    for (PetscInt i = 0; i < 3; i++){
        for (PetscInt j = 0; j <= i; j++){
            R_M[i][j] = radiusVec[i + j];
            W_M[i][j] = weightVec[i + j];
            R_M[j][i] = R_M[i][j];
            W_M[j][i] = W_M[i][j];
        }
    }
    filterType = filterT;
    
    // Call the setup method
    SetUp(da_yee_lattice, x);
}

Filter::~Filter() {
    // Deallocate data
    if (Hs != NULL) {
        VecDestroy(&Hs);
    }
    if (H != NULL) {
        MatDestroy(&H);
    }
    if (dx != NULL) {
        VecDestroy(&dx);
    }
}

// Filter design variables
PetscErrorCode Filter::FilterProject(Vec x, Vec xTilde, Vec xPhys, PetscBool projectionFilter, PetscScalar beta,
                                     PetscScalar eta) {
    PetscErrorCode ierr;

    // Filter the design variables or copy to xPhys
    // STANDARD FILTER
    if (filterType == 1) {  
        // Filter the densitities
        ierr = MatMult(H, x, xTilde);
        CHKERRQ(ierr);
        VecPointwiseDivide(xTilde, xTilde, Hs);
    }
    // COPY IN CASE OF SENSITIVITY FILTER
    else {
        ierr = VecCopy(x, xTilde);
        CHKERRQ(ierr);
    }

    // Check for projection
    if (projectionFilter) {
        HeavisideFilter(xPhys, xTilde, beta, eta);
    } else {
        VecCopy(xTilde, xPhys);
    }

    return ierr;
}

// Filter the sensitivities
PetscErrorCode Filter::Gradients(Vec x, Vec xTilde, Vec dfdx, PetscInt m, Vec* dgdx, PetscBool projectionFilter,
                                 PetscScalar beta, PetscScalar eta) {

    PetscErrorCode ierr;
    // Cheinrule for projection filtering
    if (projectionFilter) {
          
        // Get correction
        ChainruleHeavisideFilter(dx, xTilde, beta, eta);
              
        PetscScalar *xt, *dg, *df, *dxp;
        PetscInt     locsiz;

        ierr = VecGetLocalSize(xTilde, &locsiz);
        CHKERRQ(ierr);
        ierr = VecGetArray(xTilde, &xt);
        CHKERRQ(ierr);
        ierr = VecGetArray(dx, &dxp);
        CHKERRQ(ierr);
        // Objective function
        ierr = VecGetArray(dfdx, &df);
        CHKERRQ(ierr);
        for (PetscInt j = 0; j < locsiz; j++) {
            df[j] = df[j] * dxp[j];
        }
        ierr = VecRestoreArray(dfdx, &df);
        CHKERRQ(ierr);
        // Run through all constraints
        for (PetscInt i = 0; i < m; i++) {
            ierr = VecGetArray(dgdx[i], &dg);
            CHKERRQ(ierr);
            // The eta item corresponding to the correct realization
            for (PetscInt j = 0; j < locsiz; j++) {
                dg[j] = dg[j] * dxp[j];
            }
            ierr = VecRestoreArray(dgdx[i], &dg);
            CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(dx, &dxp);
        CHKERRQ(ierr);
        ierr = VecRestoreArray(dgdx[0], &dg);
        CHKERRQ(ierr);
        ierr = VecRestoreArray(xTilde, &xt);
        CHKERRQ(ierr);
    }

    // Chainrule/Filter for the sensitivities
    if (filterType == 0)
    // Filter the sensitivities, df,dg
    {
        Vec xtmp;
        ierr = VecDuplicate(xTilde, &xtmp);
        CHKERRQ(ierr);
        VecPointwiseMult(xtmp, dfdx, x);
        MatMult(H, xtmp, dfdx);
        VecPointwiseDivide(xtmp, dfdx, Hs);
        VecPointwiseDivide(dfdx, xtmp, x);
        VecDestroy(&xtmp);
    } else if (filterType == 1) {
        // Filter the densities, df,dg: STANDARD FILTER
        Vec xtmp;
        ierr = VecDuplicate(x, &xtmp);
        CHKERRQ(ierr);
        // dfdx
        VecPointwiseDivide(xtmp, dfdx, Hs);
        MatMult(H, xtmp, dfdx);
        // dgdx
        for (PetscInt i = 0; i < m; i++) {
            VecPointwiseDivide(xtmp, dgdx[i], Hs);
            MatMult(H, xtmp, dgdx[i]);
        }
        // tidy up
        VecDestroy(&xtmp);
    } 
    
    return ierr;
}

PetscScalar Filter::GetMND(Vec x) {

    PetscScalar mnd, mndloc = 0.0;

    PetscScalar* xv;
    PetscInt     nelloc, nelglob;
    VecGetLocalSize(x, &nelloc);
    VecGetSize(x, &nelglob);

    // Compute power sum
    VecGetArray(x, &xv);
    for (PetscInt i = 0; i < nelloc; i++) {
        mndloc += 4 * xv[i] * (1.0 - xv[i]);
    }
    // Collect from procs
    MPI_Allreduce(&mndloc, &mnd, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);
    mnd = mnd / ((PetscScalar)nelglob);

    return mnd;
}

PetscErrorCode Filter::HeavisideFilter(Vec y, Vec x, PetscReal beta, PetscReal eta) {
    PetscErrorCode ierr;

    PetscScalar *yp, *xp;
    PetscInt     nelloc;
    VecGetLocalSize(x, &nelloc);
    ierr = VecGetArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &yp);
    CHKERRQ(ierr);

    for (PetscInt i = 0; i < nelloc; i++) {
        yp[i] = SmoothProjection(xp[i], beta, eta);
    }
    ierr = VecRestoreArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &yp);
    CHKERRQ(ierr);
}

PetscErrorCode Filter::ChainruleHeavisideFilter(Vec y, Vec x, PetscReal beta, PetscReal eta) {
    PetscErrorCode ierr;
    
    PetscScalar *yp, *xp;
    PetscInt     nelloc;
    VecGetLocalSize(x, &nelloc);
    ierr = VecGetArray(x, &xp);
    
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &yp);
    CHKERRQ(ierr);
    
    for (PetscInt i = 0; i < nelloc; i++) {
        yp[i] = ChainruleSmoothProjection(xp[i], beta, eta);
    }
    ierr = VecRestoreArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &yp);
    CHKERRQ(ierr);
}

// Continuation function
PetscBool Filter::IncreaseBeta(PetscReal* beta, PetscReal betaFinal, PetscScalar gx, PetscInt itr, PetscReal ch) {

    PetscBool changeBeta = PETSC_FALSE;

    // Increase beta when fitting
    if ((ch < 0.01 || itr % 10 == 0) && beta[0] < betaFinal && gx < 0.000001) {
        changeBeta = PETSC_TRUE;
        if (beta[0] < 7) {
            beta[0] = beta[0] + 1;
        } else {
            beta[0] = beta[0] * 1.2;
        }
        if (beta[0] > betaFinal) {
            beta[0]    = betaFinal;
            changeBeta = PETSC_FALSE;
        }
        PetscPrintf(PETSC_COMM_WORLD, "Beta has been increased to: %f\n", beta[0]);
    }

    return changeBeta;
}

PetscErrorCode Filter::SetUp(DM da_yee_lattice, Vec dens) {

    PetscErrorCode ierr;

    VecDuplicate(dens, &dx); //dx =  Projection filter chainrule correction here!!!
     
    VecSet(dx, 1.0);
      	
    if (filterType == 0 || filterType == 1) {
        
        PetscInt        rank;                                
        DM              da_tmp;
        PetscInt        M, N, P, md, nd, pd;
        DMBoundaryType  bx, by, bz;
        DMDAStencilType stype;
        PetscScalar     dx, dy, dz;    
        DMDALocalInfo   info;
        Vec             lcoor; // coordinate vector of edges of the Yee lattice
        PetscScalar*    lcoorp;
        PetscScalar*    lcoorp2;

        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        // Extract information from the yee lattice
        ierr = DMDAGetInfo(da_yee_lattice, NULL, &M, &N, &P, &md, &nd, &pd, NULL, NULL, &bx, &by, &bz, &stype);
        ierr = DMGetCoordinatesLocal(da_yee_lattice, &lcoor); CHKERRQ(ierr);  
        ierr = VecGetArray(lcoor, &lcoorp); CHKERRQ(ierr); 
        ierr = DMDAGetLocalInfo(da_yee_lattice, &info); CHKERRQ(ierr); 
        dx = lcoorp[0 + 3                      ] - lcoorp[0]; 
        dy = lcoorp[1 + 3 * info.gxm           ] - lcoorp[1]; 
        dz = lcoorp[2 + 3 * info.gxm * info.gym] - lcoorp[2]; 
        VecRestoreArray(lcoor, &lcoorp);

        // Create the yee cells connectivity shit
        PetscInt el_tmp = 0;
        PetscInt CellsConn_max = 0;
        PetscPrintf(PETSC_COMM_WORLD, "------------------------------------\n");
        PetscPrintf(PETSC_COMM_WORLD, "Radius matrix and weight matrix (both symmetric), 0=x, 1=y, 2=z:\n");
        for (PetscInt i = 0; i < 3; i++){
            for (PetscInt j = 0; j <= i; j++){
                
                PetscInt CellsConn = 0;
                // Check dx,dy,dz and find max conn for a given rmin
                CellsConn = (PetscInt)PetscMax(ceil(R_M[i][j] / dx) - 1, PetscMax(ceil(R_M[i][j] / dy) - 1, ceil(R_M[i][j] / dz) - 1));
                CellsConn = PetscMin(CellsConn, PetscMin((M - 1) / 2, PetscMin((N - 1) / 2, (P - 1) / 2)));
                CellsConn = PetscMax(0, CellsConn);

                // The following is needed due to roundoff errors
                PetscInt tmp;
                MPI_Allreduce(&CellsConn, &tmp, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);
                CellsConn = tmp;
                el_tmp =  CellsConn; 
                if (el_tmp > CellsConn_max){
                    CellsConn_max = el_tmp;
                }
      
                // Print to screen: mesh overlap!
                PetscPrintf(PETSC_COMM_WORLD, "W_%i%i = %5.2f | R_%i%i = %5.2f [UL] -> stencil of %i yee cells \n", 
                            i, j, W_M[i][j], i, j, R_M[i][j], CellsConn);
            }
        }
        PetscPrintf(PETSC_COMM_WORLD, "-----------------\n");
        PetscPrintf(PETSC_COMM_WORLD, "To build the filter matrix a DMDA with stencil width = %i will be chosen!\n", CellsConn_max);
        PetscPrintf(PETSC_COMM_WORLD, "Building matrix... ");

        // Create spacing vector
        PetscReal sp_v[3];
        sp_v[0] = dx;
        sp_v[1] = dy;
        sp_v[2] = dz;

        // Create spacing matrix
        PetscReal sp_M[3][3];
        for (PetscInt i = 0; i < 3; i++){
            for (PetscInt j = 0; j < 3; j++){
                if (i==j){
                    sp_M[i][j] = sp_v[i];
                }
                else{
                    sp_M[i][j] = 0;
                }
            }
        }
        
        // Get number of nodes for each partition
        const PetscInt *Lx, *Ly, *Lz;
        ierr = DMDAGetOwnershipRanges(da_yee_lattice, &Lx, &Ly, &Lz); CHKERRQ(ierr);    
        
        // Check if length of subdomain >= stencil width
        const PetscInt cellsArr[3] = {md, nd, pd};
        const PetscInt* LArr[3] = {Lx, Ly, Lz};
        for (int c = 0; c < 3; c++){
            for (int i = 0; i < cellsArr[c]; i++){
                if (LArr[c][i] < CellsConn_max){
                    PetscErrorPrintf("Local width of domain is smaller than stencil width of filter! "
                                     "Choose a smaller filter size or increase/relocate your DR.\n");
                    MPI_Abort(PETSC_COMM_WORLD, 2);
                }
            }
        }    
        // Create DM again, but wit stencil width = CellsConn_max, to build the filtering matrix:
        ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, stype, M, N, P, md, nd, pd, 3, CellsConn_max, Lx, Ly, Lz,
                     &da_tmp); CHKERRQ(ierr);
        // Initialize
        ierr = DMSetFromOptions(da_tmp); CHKERRQ(ierr);
        ierr = DMSetUp(da_tmp); CHKERRQ(ierr);
        ierr = DMDASetUniformCoordinates(da_tmp, 0., (M - 1) * dx, 0., (N - 1) * dy, 0., (P - 1) * dz); CHKERRQ(ierr);
        ierr = DMGetCoordinatesLocal(da_tmp, &lcoor); CHKERRQ(ierr); 
        ierr = VecGetArray(lcoor, &lcoorp2); CHKERRQ(ierr); 
        ierr = DMDAGetLocalInfo(da_tmp, &info); CHKERRQ(ierr); 
        // Filter matrix - Allocate and assemble
        ierr = DMCreateMatrix(da_tmp, &H); CHKERRQ(ierr); 
        ierr = DMCreateGlobalVector(da_tmp, &Hs); CHKERRQ(ierr); 

        // The variables from info that are used are described below:
        // -------------------------------------------------------------------------
        // sw = Stencil width
        // mx, my, mz = Global number of "elements" in each direction
        // xs, ys, zs = Starting point of this processor, excluding ghosts
        // xm, ym, zm = Number of grid points on this processor, excluding ghosts
        // gxs, gys, gzs = Starting point of this processor, including ghosts
        // gxm, gym, gzm = Number of grid points on this processor, including ghosts
        // -------------------------------------------------------------------------
        // Outer loop is local part = find row
        // What is done here, is:
        //
        // 1. Run through all elements in the mesh - should not include ghosts
        for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                    // The row number of the yee cells we are considering:
                    PetscInt cell_idx1 =
                        (i - info.gxs) + (j - info.gys) * (info.gxm) + (k - info.gzs) * (info.gxm) * (info.gym);

                    // 2. Loop over nodes (including ghosts) within a cubic domain with
                    // center at (i,j,k)
                    //    For each element, run through all elements in a box of size
                    //    stencilWidth * stencilWidth * stencilWidth Remark, we want to
                    //    make sure we are not running "out of the domain", therefore k2
                    //    etc. are limited to the max global index (info.mz-1 etc.)
                    
                    for (PetscInt comp1 = 0;  comp1 < 3; comp1++){
                        PetscInt row = 3 * cell_idx1 + comp1; 
                        for (PetscInt k2 = PetscMax(k - info.sw, 0); k2 <= PetscMin(k + info.sw, info.mz - 1); k2++) {
                            for (PetscInt j2 = PetscMax(j - info.sw, 0); j2 <= PetscMin(j + info.sw, info.my - 1); j2++) {
                                for (PetscInt i2 = PetscMax(i - info.sw, 0); i2 <= PetscMin(i + info.sw, info.mx - 1); i2++) {
                                    // The col number of the yee cells we are considering:
                                    PetscInt cell_idx2 = (i2 - info.gxs) + (j2 - info.gys) * (info.gxm) +
                                                    (k2 - info.gzs) * (info.gxm) * (info.gym);
                                    
                                    for (PetscInt comp2 = 0;  comp2 < 3; comp2++){
                                        PetscInt col = 3 * cell_idx2 + comp2; 
                                        PetscScalar dist = 0.0; // distance
                             
                                        // Compute the distance from the dens_comp1 in the "row" Yee Cell
                                        // to the dens_comp2 in the "col" Yee Cell 
                                        for (PetscInt kk = 0; kk < 3; kk++) {
                                            dist = dist + PetscPowScalar((lcoorp2[3 * cell_idx1 + kk] + 0.5 * sp_M[comp1][kk]) 
                                                                       - (lcoorp2[3 * cell_idx2 + kk] + 0.5 * sp_M[comp2][kk]), 2.0);
                                        }        
                                        dist = PetscSqrtScalar(dist); // why memory error in MatSetValuesLocal, if we neglect that operation???????
                                        if (dist < R_M[comp1][comp2]) {
                                            // Longer distances should have less weight                           
                                            dist = W_M[comp1][comp2] * (R_M[comp1][comp2] - dist); // weight value 
                                            ierr = MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES); CHKERRQ(ierr); //????
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
  
        // Assemble H:
        ierr = MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
        ierr = MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
        // Compute the Hs, i.e. sum the rows
        Vec dummy;
        ierr = VecDuplicate(Hs, &dummy); CHKERRQ(ierr); 
        ierr = VecSet(dummy, 1.0); CHKERRQ(ierr); 
        ierr = MatMult(H, dummy, Hs); CHKERRQ(ierr); 

        // Clean up
        VecRestoreArray(lcoor, &lcoorp);
        VecDestroy(&dummy);
        DMDestroy(&da_tmp);

        PetscPrintf(PETSC_COMM_WORLD, "done!\n");
        PetscPrintf(PETSC_COMM_WORLD, "------------------------------------\n");
    } 

    return ierr;
}
