import numpy as np
from matplotlib import pyplot as plt
import SpinDens3Particle as SpinDensity
import sys
import re
import os
import cubetools


ConvAbsDens      = 1e-8
octExecute       = "mpirun -n 2 octopus > out.log 2>&1" #on theo109
#octExecute       = "mpirun ~ahhsw/software/progs/octopus/2020.05.19/bin/octopus > out.log 2>&1" #on cluster
#octExecute_4D    = "mpirun ~ahhsw/software/progs/octopus/2020.05.19/bin/octopus-maxdim4 > out.log 2>&1"  
nucs_Z, nucs_pos = None, None

def Create_InpFile(P, v_ext_return, plot, extra_states = 2):
    nucs_Z   = P['nucs_array'][0, :]
    nucs_pos = P['nucs_array'][1, :]

    '''_________POTENTIAL STRING_______________________________________________________________________________________'''
    def pot_str(nucs_Z, nucs_pos):
        sp_var = ['x', 'y', 'z', 'w']
        sp_var = sp_var[0:P['dim']]
        v_ee  = []
        v_ext = []
        for i in range(P['dim']):
            for j in range(0, i):
                v_ee.append('1/(' + str(P['a']) + '+(' + sp_var[i] + '-' + sp_var[j] + ')^2)^(1/2)')
            for k in range(P['nucs']):

                v_ext.append(nucs_Z[k] + '/(' + str(P['a']) + '+(' + nucs_pos[k] + '-' + sp_var[i] + ')^2)^(1/2)')

        v_ee  = " + ".join(v_ee)
        v_ext = " - ".join(v_ext)
        v = '"' + v_ee + ' - ' +  v_ext + '"'
        #print('potential = ', v)
        return v

    deleteFiles()
    '''#_________CREATE INP FILE__________#
    '''
    inp     = open('inp', 'w')
    inp.write("ExperimentalFeatures = yes\n")
    inp.write("FromScratch = yes\n\n")

    inp.write("CalculationMode = gs\n")
    inp.write("ConvAbsDens = " + str(ConvAbsDens) + "\n")
    inp.write("Dimensions = " + str(P['dim']) + "\n")
    inp.write("TheoryLevel = independent_particles\n")
    inp.write("ExtraStates = " + str(extra_states) + "\n\n")
    if P['dim'] < 4:
        inp.write("BoxShape = parallelepiped\n")
    else:
        inp.write("BoxShape = hypercube\n")
    inp.write("Lsize = " + str(P['radius']) + "\n")
    inp.write("Spacing = " + str(P['dx']) + "\n\n")
    inp.write("%Species\n")
    inp.write('  "mpot" | species_user_defined | potential_formula |'+pot_str(nucs_Z, nucs_pos) + '| valence | 1 $\n')
    inp.write("%\n\n%Coordinates\n")
    inp.write('  "mpot" '+P['dim']*"| 0 " + "\n")
    inp.write("\n%\n\n")
    inp.write("NDimModelmb = 1\n")
    inp.write("NParticleModelmb = " + str(P['dim']) + "\n")
    inp.write("%DescribeParticlesModelmb\n")
    inp.write(P['dim']*'   "electron" | 1 | 1. | 1. | fermion\n')
    inp.write("\n%\n%DensitytoCalc\n")
    inp.write('   "electron" | 1 | -1\n%\n\n')
    if P['dim'] < 3:
        inp.write("Output = mmb_den\n")
        inp.write("OutputFormat = axis_x")
    else:
        inp.write("Output = mmb_wfs\n")
        inp.write("OutputFormat = cube")

    inp.close()

    if v_ext_return:
        nucs_array = np.array([nucs_Z, nucs_pos], dtype=float)
        #print('nucs_array = ', nucs_array)
        points = np.arange(-P['radius'], P['radius'] + P['dx'], P['dx'])
        v_arr = np.zeros(len(points))
        for k in range(nucs_array.shape[1]):
            v_arr += - nucs_array[0,k] / np.sqrt(P['a'] +(nucs_array[1,k] - points)**2)
        if plot:
            plt.plot(points, v_arr)
            plt.show()
        return v_arr

def deleteFiles():
    os.system("rm -r exec")
    os.system("rm -r output_iter")
    os.system("rm -r restart")
    os.system("rm -r static")
    os.system("rm -r inp")
    os.system("rm -r out.log")

def return_Dens(P, extra_states, exc_mode):
    output_num = None
    nspindown  = []
    energies   = []
    exS_index  = []
    Energy_2dExState = None

    '''____________OCTOPUS CALCULATION____________________________________________________________________________'''

    if P['dim'] < 4:
        os.system(octExecute)
    else:
        os.system(octExecute_4D)

    if 'FATAL ERROR' in open('out.log').read():
        raise Exception('FATAL ERROR')
    if 'SCF *not* converged after' in open('out.log').read():
        raise Exception('Not converged! System will be dropped.')

    '''___________READ ENERGY FROM OUTPUT_________________________________________________________________________'''
    with open('static/modelmb/youngprojections', 'r') as youngfile:
        GS_found = False
        youngs = youngfile.read()
        sign   = 1
        for i in range(extra_states + 1):
            line_pos = re.findall("(?<= {6}" + str(i + 1) +" {7}).*",  youngs) #for positiv energies
            line_neg = re.findall("(?<= {6}" + str(i + 1) +" {6}-).*", youngs) #for negativ energies
            if not line_pos:
                line_split = line_neg[0].split()
                sign = -1
            else:
                line_split = line_pos[0].split()

            if line_split[-1] != "diagram":
                if not GS_found:
                    output_num = i + 1
                    GS_found = True
                energies.append(sign * float(line_split[0]))
                nspindown.append(int(line_split[-2]))

                if exc_mode == False and P['dim'] != 2: # additional loop for Spin_Koopmans Theorem needed
                    break
                else:
                    exS_index.append(output_num)

        if len(energies) == 0:
            raise Exception('No associated Young diagram found!')

    '''___________READ DENSITY FROM OUTPUT________________________________________________________________________'''
    if exc_mode == False:
        if P['dim'] < 3:
            points, Dens_data_total = np.loadtxt(r'static/modelmb/density_ip001_imb0' + str(output_num), unpack=True)
            if P['dim'] == 1:
                Dens_data_up   = Dens_data_total
                Dens_data_down = np.zeros(len(points))
            if P['dim'] == 2:
                Dens_data_up   = Dens_data_total  / 2.
                Dens_data_down = Dens_data_total  / 2.

        else:
            wf_YP1, meta1 = cubetools.read_cube(r'static/modelmb/wf-st0002.cube') # projected wf onto 1st young diagram
            wf_YP2, meta2 = cubetools.read_cube(r'static/modelmb/wf-st0003.cube') # projected wf onto 2nd young diagram
            dens_total_wf_YP1, dens_up_wf_YP1, dens_down_wf_YP1 = SpinDensity.SpinDens(wf_YP1, dx=P['dx'])
            dens_total_wf_YP2, dens_up_wf_YP2, dens_down_wf_YP2 = SpinDensity.SpinDens(wf_YP2, dx=P['dx'])

            Dens_data_total       = dens_total_wf_YP1
            Dens_data_up   = dens_up_wf_YP1
            Dens_data_down = dens_down_wf_YP1


        excited_states_no    = [str(output_num) + " (" + "?" + ")"]
        excited_nspindown    = [nspindown]
    else:
        raise Exception('SPINDENS CREATION ONLY AVAILABLE ONLY FOR GROUND STATE!')

    if P['dim'] != 3:
        return Dens_data_total , Dens_data_up, Dens_data_down, \
               energies, excited_states_no, excited_nspindown,\
               None, None, None, \
               None, None, None
    if P['dim'] == 3:
        return Dens_data_total , Dens_data_up, Dens_data_down, \
               energies, excited_states_no, excited_nspindown, \
               dens_total_wf_YP1, dens_up_wf_YP1, dens_down_wf_YP1, \
               dens_total_wf_YP2, dens_up_wf_YP2, dens_down_wf_YP2

























