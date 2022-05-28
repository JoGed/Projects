import numpy as np
import cubetools


# code to extract the spin densities from the wf ouput of octopus for 3 electrons

'''HOW TO USE:
Compute wave functions for N=3 using inp file - e.g.:

__________________inp file____________________________
ExperimentalFeatures = yes
FromScratch = yes

CalculationMode = gs
ConvAbsDens = 1e-08
Dimensions = 3
TheoryLevel = independent_particles
ExtraStates = 2

BoxShape = parallelepiped
Lsize = 11.5
Spacing = 0.1

%Species
  "mpot" | species_user_defined | potential_formula |"1/(1+(y-x)^2)^(1/2) + 1/(1+(z-x)^2)^(1/2) + 1/(1+(z-y)^2)^(1/2) - 1/(1+(-3.198824952949975-x)^2)^(1/2) - 1/(1+(-1.1725911208414486-x)^2)^(1/2) - 1/(1+(-3.220188170087039-x)^2)^(1/2) - 1/(1+(-3.198824952949975-y)^2)^(1/2) - 1/(1+(-1.1725911208414486-y)^2)^(1/2) - 1/(1+(-3.220188170087039-y)^2)^(1/2) - 1/(1+(-3.198824952949975-z)^2)^(1/2) - 1/(1+(-1.1725911208414486-z)^2)^(1/2) - 1/(1+(-3.220188170087039-z)^2)^(1/2)"| valence | 1 $
%

%Coordinates
  "mpot" | 0 | 0 | 0 

%

NDimModelmb = 1
NParticleModelmb = 3
%DescribeParticlesModelmb
   "electron" | 1 | 1. | 1. | fermion
   "electron" | 1 | 1. | 1. | fermion
   "electron" | 1 | 1. | 1. | fermion

%
%DensitytoCalc
   "electron" | 1 | -1
%

Output = mmb_wfs
OutputFormat = cube
____________________________________________________________________________

Use this code to extract the spin densities via:
wf_YP1, meta1 = cubetools.read_cube(r'static/modelmb/wf-st0002.cube') # projected wf onto 1st young diagram (you can choose 2nd diagram as well)
Dens_data_total, Dens_data_up , Dens_data_down  = SpinDensity.SpinDens(wf_YP1, dx=0.1) # Here, a spacing of dx=0.1 has been assumed!
'''

def Perm(wf, Permutation):
	
	if Permutation[0] == Permutation[1]:
		return wf
		
	fixed_var = np.delete([1,2,3], np.array(Permutation) - 1)
	wf_new = wf.copy()
	
	for i in range(wf.shape[0]):
		if fixed_var == 1:	
			wf_new[i,:,:] = wf[i,:,:].T
		elif fixed_var == 2:
			wf_new[:,i,:] = wf[:,i,:].T
		elif fixed_var == 3:
			wf_new[:,:,i] = wf[:,:,i].T
	return wf_new


def SpinDens(wf, dx=1):
	
	Theta_1 =  ( 2 * wf - 2 * Perm(wf, [1,2]) + Perm(wf, [1,3]) + Perm(wf, [2,3]) - Perm(Perm(wf, [2,3]), [1,2]) - Perm(Perm(wf, [1,2]), [2,3]))
	Theta_2 =  + ( Perm(wf, [1,3]) - Perm(wf, [2,3]) - Perm(Perm(wf, [2,3]), [1,2]) + Perm(Perm(wf, [1,2]), [2,3]))
	
	M   = np.zeros((2,2, *wf.shape), dtype= float)

	M[0][0] = 1 *  np.multiply(Theta_1, Theta_1)  
	M[0][1] = 3 *  np.multiply(Theta_1, Theta_2)  
	M[1][0] = 3 *  np.multiply(Theta_2, Theta_1)  
	M[1][1] = 9 *  np.multiply(Theta_2, Theta_2)  

	M_x = np.sum(M.copy(), axis = (3,4))
	M_y = np.sum(M.copy(), axis = (2,4))
	M_z = np.sum(M.copy(), axis = (2,3))

	x_up_coeffs   = np.array([[ 5,-1],
				   [-1, 1]
				   ])

	x_down_coeffs = np.array([[ 1, 1],
				   [ 1, 1]
				   ])

	y_up_coeffs   = np.array([[ 5, 1],
				   [ 1, 1]
				   ])

	y_down_coeffs = np.array([[ 1,-1],
				   [-1, 1]
				   ])

	z_up_coeffs   = np.array([[ 2, 0],
				   [ 0, 2]
				   ])

	z_down_coeffs = np.array([[ 4, 0],
				   [ 0, 0]
				   ])

	dens_x_up   = np.sum(np.multiply(np.expand_dims(x_up_coeffs,   axis = -1), M_x), axis = (0,1))
	dens_x_down = np.sum(np.multiply(np.expand_dims(x_down_coeffs, axis = -1), M_x), axis = (0,1))

	dens_y_up   = np.sum(np.multiply(np.expand_dims(y_up_coeffs,   axis = -1), M_y), axis = (0,1))
	dens_y_down = np.sum(np.multiply(np.expand_dims(y_down_coeffs, axis = -1), M_y), axis = (0,1))

	dens_z_up   = np.sum(np.multiply(np.expand_dims(z_up_coeffs,   axis = -1), M_z), axis = (0,1))
	dens_z_down = np.sum(np.multiply(np.expand_dims(z_down_coeffs, axis = -1), M_z), axis = (0,1))

	dens_up     = dens_x_up   + dens_y_up   + dens_z_up 
	dens_down   = dens_x_down + dens_y_down + dens_z_down 
	

	Norm_fac    = (dens_up.sum() + dens_down.sum()) * dx

	dens_up   *= 3. / Norm_fac
	dens_down *= 3. / Norm_fac
	dens_total = dens_up + dens_down

	return dens_total, dens_up, dens_down
	

