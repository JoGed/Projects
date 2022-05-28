import numpy as np
import csv
import sys
import KohnShamSpin
import datetime
import pandas as pd
from prettytable import PrettyTable
from tqdm import tqdm

spin = True
alpha = np.array([1.]) # fractional densities

if len(sys.argv) < 5:
    print("_________HOPPLA!__________________________________________________")
    print("Arguments are missing! Type 'python V_xc_Creator.py <index> <float1> <float2> <string>'', \n"
          "where 'index' corresponds to the DensData_<index>.csv file.\n"
          "Set conditions for convergences: MSE < 'float1' & maxDensDiff < 'float2'.\n"
          "'string' specifies the calculation mode, e.g.: g - groundstate | e - excited mode.")
    print("__________________________________________________________________")
    sys.exit()

def  GS_dens_splitter(particles):
    if particles < 1:
        raise Exception("particles < 1!")
    rounded, append = int(particles), particles - int(particles)
    if rounded % 2 == 0:
        s = int(rounded / 2.)
        up_occ     = np.ones(s + 1)
        up_occ[-1] = append
        down_occ = np.ones(s + 1)
        down_occ[-1] = 0
    else:
        s = int((rounded - 1) / 2.)
        up_occ = np.ones(s + 1)
        down_occ = np.ones(s + 1)
        down_occ[-1] = append

    up_occ   = np.concatenate((up_occ, [0]))
    down_occ = np.concatenate((down_occ, [0]))
    
    return  up_occ + down_occ, up_occ, down_occ

def Consistency_verifier(arr, max):
    if np.mean(arr) != np.mean(np.arange(1, max + 1, 1)):
        print("Inconcistency found!")
        sys.exit()

# Conditions for convergences: MSE < cond_MSE & maxDensDiff < cond_maxDiff
cond_MSE     = np.prod(np.array(sys.argv[2].replace("*", " ").split()).astype(float))
cond_maxDiff = np.prod(np.array(sys.argv[3].replace("*", " ").split()).astype(float))
csv_index    = str(sys.argv[1])

if sys.argv[4] == 'g':
    csv_input    = r"DensSets/DensData_" + str(csv_index) + ".csv"
    csv_ouput    = r"TrainingsSets/TrainingsSet_" + str(csv_index) + ".csv"
    dim_collumn  = np.array(pd.read_csv(csv_input)['dim'])
    nuc_collumn  = np.array(pd.read_csv(csv_input)['nucs_array'])
    max_dim      = np.amax(dim_collumn)
else:
    raise Exception('Works only for groundstate mode yet')

Consistency_verifier(dim_collumn.astype(int), max_dim)
data_lenght  = len(dim_collumn)
data_lenght2 = (len(dim_collumn) / 3) * (1 + len(alpha) * 2)


systems      = []
test_x = []
test_y = []
def floatArray(string):
    return np.array(string.replace('[', '').replace(']', '').split()).astype(float)

with open(csv_input) as csvfile:
    reader       = csv.DictReader(csvfile)
    start_time   = datetime.datetime.now().replace(microsecond=0)
    counter      = 0
    failure      = False
    sample       = (1 + len(alpha) * (max_dim - 1)) * [None]
    Energies_integerDim = max_dim * [None]
    dropped      = 0
    failed_dims  = []
    failed_nucs  = []

    error_table = PrettyTable()
    error_table.field_names = ["row", "at dim", "nucs_array", "Error"]

    last_int_particle_rows = None
    sample_idx = 0
    idx_newCSV = 0
    Energy_2dExState = None

    for i, row in enumerate(reader):

        counter += 1 # integer particle number

        if (counter > max_dim):
            print("_______________________")
            if (None in sample) == False:
                for s in sample:
                    systems.append(s)
            sample       = (1 + len(alpha) * (max_dim - 1)) * [None]
            Energies_integerDim = max_dim * [None]
            counter      = 1
            Energy_2dExState = None
            sample_idx   = 0
            failure      = False
            print("\n\n")

        if failure == False:

            #print("#########################################################################")
            for a in range(len(alpha)):
                sample_idx += 1
                idx_newCSV += 1

                '''________________________PARAMETERS_________________________________________________________________________________________________________________'''
                p = {
                    # ______SYSTEM PROPERTIES______________________________________________________________________________________________________________
                    'dim':             int(row['dim']),              # number of particles
                    'radius':          row['radius'],                # 2 * radius = 'box_length' in 1D
                    'dx':              float(row['dx']),             # grid spacing
                    'N':               int(row['N']),                # number of spatial points
                    'nucs_array':      row['nucs_array'],            # array containg [charge_arr , position_arr] (dtype = str) of each nuclei
                    'total_charge':    row['total_charge'],          # total charge of the nuclei
                    'nuc_pos_radius':  row['nuc_pos_radius'],        # max_range of positions of nuclei
                    'a':               row['a'],                     # Coulomb parameter
                    'nucs':            row['nucs'],                  # number of nuclei
                    'Energy':          float(row['Energy']),         # energy of the groundstate
                    'state':           row['state'],                 # state index
                    'nspindown':       row['nspindown'],             # number of spin downs
                    'points':          floatArray(row['points']),    # spatial points
                    'v_ext':           floatArray(row['v_ext']),     # external potential
                    'Dens_data_total': floatArray(row['Dens_data_total']), # exact Density (calculated by octopus)
                    'Dens_data_up':    floatArray(row['Dens_data_up']),  # exact Density_up (calculated by octopus)
                    'Dens_data_down':  floatArray(row['Dens_data_down']),  # exact Density_down (calculated by octopus)

                    'dens3_total_wf_YP1': floatArray(row['dens3_total_wf_YP1']), # total dens (down + up), projected wf on 1st young diagram as primitive spatial part used
                    'dens3_up_wf_YP1':    floatArray(row['dens3_up_wf_YP1']),   # dens_up, projected wf on 1st young diagram as primitive spatial part used
                    'dens3_down_wf_YP1':  floatArray(row['dens3_down_wf_YP1']), # dens_down, projected wf on 1st young diagram as primitive spatial part used

                    'dens3_total_wf_YP2': floatArray(row['dens3_total_wf_YP2']),# total dens (down + up), projected wf on 2nd young diagram as primitive spatial part used
                    'dens3_up_wf_YP2':    floatArray(row['dens3_up_wf_YP2']),# dens_up, projected wf on 2nd young diagram as primitive spatial part used
                    'dens3_down_wf_YP2':  floatArray(row['dens3_down_wf_YP2']),# dens_down, projected wf on 2nd young diagram as primitive spatial part used

                    'Energy_2dExState': None if len(row['Energy_2dExState']) == 0 else float(row['Energy_2dExState']),
                    # Energy of 1st excited state for 2 particles (needed for Koopman's Theorem for open Shells)

                    #______INVERSION PARAMETERS______________________________________________________________________________________________________________
                    'LODF':            1,                           # occupation factor (last orbital)
                    'laplOrder':       4,                           # accuracy of laplace operator
                    'dirichlet':       0.,                          # fixes v_xc to intial values for dens <= dirichlet
                    'weight':          1,                           # weight_function for MSE Loss
                    'v_xc_init':       'FA',                        # initial v_xc. Choose 'FA' for Fermi–Amaldi potential that satisfies corrrect -1/r decay,
                                                                    # otherwise use an arbitrary lamba function (P['Dens_data'] will be input)

                    # ______TO CALCULATE _____________________________________________________________________________________________________________________
                    'occ':             [],                          # occupation of KS orbitals, e.g. [1,1] means, that the 1st and 2nd orbital is occupied once!
                    'evals':           None,                        # Kohn Sham eigenvalues
                    'e_h':             None,                        # highest occupied Kohn Sham eigenvalue
                    'Dens_inv':        None,                        # density (by inversion)
                    'Psi':             None,                        # KohnSham orbitals (by inversion)
                    'v_xc':            None,                        # v_xc (by inversion)
                    'v_H':             None,                        # Hartree potential
                    'MSE':             None,                        # value of Loss function (MSE) at minimum
                    'E_kin':           None,                        # kinetic Energy
                    'E_ext':           None,                        # external Energy
                    'E_H':             None,                        # Hartree Energy
                    'E_xc_total(kin)':      None,                        # Exchange Correlation Energy
                    'E_xc_total(evals)': None,              # Exchange Correlation Energy


                    'spin':            spin,                         # Spin Calculation mode
                    'E_xc_spin(kin)':   None,  # Exchange Correlation Energy
                    'E_xc_spin(evals)': None,  # Exchange Correlation Energy

                    'up_occ':          [],
                    'evals_up':        None,
                    'e_h_up':          None,
                    'Dens_inv_up':     None,
                    'v_xc_up':         None,
                    'MSE_up':          None,
                    'E_kin_up':        None,
                    'Psi_up':          None,

                    'down_occ':        [],
                    'evals_down':      None,
                    'e_h_down':        None,
                    'Dens_inv_down':   None,
                    'v_xc_down':       None,
                    'MSE_down':        None,
                    'E_kin_down':      None,
                    'Psi_down':        None,

                }
                Energies_integerDim[counter - 1] = p['Energy']
                points = p['points']


                if p['dim'] == 2:
                    Energy_2dExState = float(row['Energy_2dExState'])

                if counter > 1:
                    p["dim"]             = int(last_int_particle_rows["dim"]) + alpha[a]
                    p['Energy']          = (1 - alpha[a]) * float(last_int_particle_rows["Energy"])               + alpha[a] * p["Energy"]
                    p['Dens_data_total'] = (1 - alpha[a]) * floatArray(last_int_particle_rows['Dens_data_total']) + alpha[a] * p['Dens_data_total']
                    p['Dens_data_up']    = (1 - alpha[a]) * floatArray(last_int_particle_rows['Dens_data_up'])    + alpha[a] * p['Dens_data_up']
                    p['Dens_data_down']  = (1 - alpha[a]) * floatArray(last_int_particle_rows['Dens_data_down'])  + alpha[a] * p['Dens_data_down']

                #test
                test_x.append(p["dim"])
                test_y.append(p['Energy'])

                # test
                p['occ'], p['up_occ'], p['down_occ'] = GS_dens_splitter(p['dim'])
                print("________________________________________________________________________")
                print("dim = ", p["dim"], " | i = ", i, " | counter = ", counter, " | a = ", a, " | alpha = ", alpha[a] if counter > 1 else 0, " | sample_idx", sample_idx)
                print("GS_occ: ", p['occ'], " | Up_occ: ", p['up_occ'], " | Down_occ: ", p['down_occ'])
                print("nucs_array: ", p['nucs_array'])
                print(".........................................................................")
                ''' _________________________________________________________________________________________________________________________________________'''

                #if (p['Dens_data'][0] > p['dirichlet']) or  (p['Dens_data'][-1] > p['dirichlet']):
                #    print('-----------------------------------------------------------------------------')
                #    print('**** WARNING ****')
                #    print('Density does not vanish regarding dirichlet conditions (', p['dirichlet'], ')')
                #    print('-----------------------------------------------------------------------------')

                if p['v_xc_init'] == 'FA':
                    p['v_xc_init'] = -KohnShamSpin.V_Hartree(p, p['Dens_data_total']) / p['dim']
                else:
                    p['v_xc_init'] = row['v_xc_init'](p['Dens_data_total'])

                '''_________INVERSION___________________________________________________________________________________________________________________________________'''
                #try:
                # if p['Dens_data'][0] > 1e-6 or p['Dens_data'][-1] > 1e-6:
                #    raise Exception(
                #        'Density > 1e-6 at boundaries! (' + str(p['Dens_data'][0]) + ", " + str(p['Dens_data'][-1]) + ")")
                try:
                    '''_________________EXCHANGE CORRELATION POTENTIAL CALCULATION______________________________________________'''
                    print("Inverting (spin_indep)...")
                    evals_inv, Psi_inv, Dens_inv, v_xc_inv, v_H, MSE, maxErr \
                        = KohnShamSpin.Inversion(p, occ = p['occ'], Dens_spin = p['Dens_data_total'])

                    for j in range(len(evals_inv)):
                        if p['occ'][-1 - j] != 0:
                            p['e_h']  = evals_inv[-1 - j]
                            break

                    if float(MSE) >= cond_MSE or float(maxErr) >= cond_maxDiff:
                        raise Exception('Not converged! (MSE = ' + str(MSE) + ", maxDensDiff = " + str(maxErr) + ")")

                    if spin:
                        print("Inverting (up) ...")
                        #print("p['up_occ']: ", p['up_occ'], "p['down_occ']", p['down_occ'])
                        evals_inv_up, Psi_inv_up, Dens_inv_up, v_xc_inv_up, v_H_up, MSE_up, maxErr_up \
                            = KohnShamSpin.Inversion(p, occ = p['up_occ'], Dens_spin = p['Dens_data_up'])

                        if float(MSE_up) >= cond_MSE or float(maxErr_up) >= cond_maxDiff:
                            raise Exception('Not converged! (MSE_up = ' + str(MSE_up) + ", maxDensDiff_up = " + str(maxErr_up) + ")")

                        for j in range(len(evals_inv_up)):
                            if p['up_occ'][-1 - j] != 0:
                                p['e_h_up'] = evals_inv_up[-1 - j]
                                break

                        if p['down_occ'].sum() < 1e-3:
                            evals_inv_down, Psi_inv_down, Dens_inv_down, v_xc_inv_down, v_H_down, MSE_down, maxErr_down \
                            = np.zeros(1), np.zeros((len(p["points"]),1)), np.zeros((len(p["points"]),)), np.zeros((len(p["points"]),)), np.zeros((len(p["points"]),)), np.zeros(1), np.zeros(1)

                        else:
                            print("Inverting (down) ...")
                            evals_inv_down, Psi_inv_down, Dens_inv_down, v_xc_inv_down, v_H_down, MSE_down, maxErr_down \
                            = KohnShamSpin.Inversion(p, occ = p['down_occ'], Dens_spin = p['Dens_data_down'])

                        if float(MSE_down) >= cond_MSE or float(maxErr_down) >= cond_maxDiff:
                            raise Exception('Not converged! (MSE_down = ' + str(MSE_down) + ", maxDensDiff_down = " + str(maxErr_down) + ")")

                        for j in range(len(evals_inv_down)):
                            if p['down_occ'][-1 - j] != 0:
                                p['e_h_down'] = evals_inv_down[-1 - j]
                                break

                    '''_________________KOOPMAN CORRECTION______________________________________________________________________'''
                    shift = None
                    if counter == 1:
                        shift = Energies_integerDim[0] - p['e_h']
                    else:
                        shift = Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2] - p['e_h']

                    v_xc_inv  += shift
                    evals_inv += shift


                    if spin:
                        shift_up   = None
                        shift_down = None

                        if counter == 1:
                            shift_up   = Energies_integerDim[0] - p['e_h_up']
                            shift_down = 0 # because V_xc = 0

                        if counter == 2:
                            shift_up   = Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2] - p['e_h_up']
                            shift_down = Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2] - p['e_h_down']

                        if counter == 3:
                            shift_up   = Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2] - p['e_h_up']
                            shift_down = Energies_integerDim[counter - 1] - Energy_2dExState          - p['e_h_down']

                        v_xc_inv_up    += shift_up
                        v_xc_inv_down  += shift_down
                        evals_inv_up   += shift_up
                        evals_inv_down += shift_down

                    '''________________EXCHANGE CORRELATION ENERGY CALCULATION__________________________________________________'''
                    v_H_diag = KohnShamSpin.V_Hartree(p, p['Dens_data_total'])
                    v_ext_diag = p['v_ext']

                    E_kin      = KohnShamSpin.E_kinetic(p, Psi_inv, occ=p['occ'])
                    E_ext, E_H = KohnShamSpin.Energies(p, p['Dens_data_total'], v_ext_diag, v_H_diag, occ = p['occ'])
                    E_xc_total_kin       = p['Energy'] - E_ext - E_H - E_kin

                    if spin:
                      E_kin_up = KohnShamSpin.E_kinetic(p, Psi_inv_up, p['up_occ'])

                      if p['down_occ'].sum() < 1e-3:
                        E_kin_down = 0
                      else:
                        E_kin_down = KohnShamSpin.E_kinetic(p, Psi_inv_down, p['down_occ'])
                      E_xc_spin_kin  = p['Energy'] - E_ext - E_H - (E_kin_up + E_kin_down)

                    if sys.argv[-1] == "test":
                        print(">>>>>>>>>>>> TEST <<<<<<<<<<<<")
                        v_eff_diag = v_H_diag + v_ext_diag + v_xc_inv
                        evals, Psi, D_Matrix = KohnShamSpin.Orbitals(P=p, v_eff_diag=v_eff_diag, occ=p["occ"])
                        v_eff_diag_up = v_H_diag + v_ext_diag + v_xc_inv_up
                        evals_up, Psi_up, D_Matrix_up = KohnShamSpin.Orbitals(P=p, v_eff_diag=v_eff_diag_up, occ=p["up_occ"])
                        v_eff_diag_down = v_H_diag + v_ext_diag + v_xc_inv_down
                        evals_down, Psi_down, D_Matrix_down = KohnShamSpin.Orbitals(P=p, v_eff_diag=v_eff_diag_down, occ=p["down_occ"])
                        if counter == 1:
                            print("evals:      ", evals,      " | Energies_integerDim[0]: ", Energies_integerDim[0])
                            print("evals_up:   ", evals_up,   " | Energies_integerDim[0]: ", Energies_integerDim[0])
                        if counter == 2:
                            print("evals:      ", evals,      " | Energies_integerDim[1] - Energies_integerDim[0]: ", Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2])
                            print("evals_up:   ", evals_up,   " | Energies_integerDim[1] - Energies_integerDim[0]: ", Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2])
                            print("evals_down: ", evals_down, " | Energies_integerDim[1] - Energies_integerDim[0]: ", Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2])
                        if counter == 3:
                            print("evals:      ", evals,      " | Energies_integerDim[2] - Energies_integerDim[1]: ", Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2])
                            print("evals_up:   ", evals_up,   " | Energies_integerDim[2] - Energies_integerDim[1]: ", Energies_integerDim[counter - 1] - Energies_integerDim[counter - 2])
                            print("evals_down: ", evals_down, " | Energies_integerDim[2] - Energy_2dExState:       ", Energies_integerDim[counter - 1] - Energy_2dExState)
                    # ____TEST________________

                    E_eigen, E_vxc, E_H = KohnShamSpin.Energies2(p, p['Dens_data_total'], evals_inv, v_xc_inv, v_H_diag, occ = p['occ'])
                    E_eigen_up, E_vxc_up, E_H_up_WRONG = KohnShamSpin.Energies2(p, p['Dens_data_up'], evals_inv_up, v_xc_inv_up, v_H_diag, occ=p['up_occ'])
                    E_eigen_down, E_vxc_down, E_H_down_WRONG = KohnShamSpin.Energies2(p, p['Dens_data_down'], evals_inv_down, v_xc_inv_down, v_H_diag, occ=p['down_occ'])
                    E_xc_total_evals = p['Energy'] - E_eigen + E_vxc + E_H
                    E_xc_spin_evals  = p['Energy'] - (E_eigen_up + E_eigen_down) + (E_vxc_up + E_vxc_down) + E_H
                    print("Energy                 : ", Energies_integerDim[counter - 1])
                    print("E_xc_total (with E_kin): ", E_xc_total_kin)
                    print("E_xc_total (with evals): ", E_xc_total_evals)
                    print("E_xc_spin  (with E_kin): ", E_xc_spin_kin)
                    print("E_xc_spin  (with evals): ", E_xc_spin_evals)


                    '''___________________________SAVING PARAMETERS______________________________________________________________'''
                    p['evals']             = evals_inv
                    p['Dens_inv']          = Dens_inv
                    p['v_xc']              = v_xc_inv
                    p['v_H']               = v_H
                    p['MSE']               = MSE
                    p['maxDensDiff']       = maxErr
                    p['E_kin']             = E_kin
                    p['E_ext']             = E_ext
                    p['E_H']               = E_H
                    p['E_xc_total(kin)']   = E_xc_total_kin
                    p['E_xc_total(evals)'] = E_xc_total_evals
                    p['Psi']               = Psi_inv


                    if spin:
                        p['E_xc_spin(kin)']    = E_xc_spin_kin
                        p['E_xc_spin(evals)']  = E_xc_spin_evals

                        p['Dens_inv_up']       = Dens_inv_up
                        p['evals_up']          = evals_inv_up
                        p['v_xc_up']           = v_xc_inv_up
                        p['MSE_up']            = MSE_up
                        p['E_kin_up']          = E_kin_up
                        p['Psi_up']            = Psi_inv_up

                        p['Dens_inv_down']     = Dens_inv_down
                        p['evals_down']        = evals_inv_down
                        p['v_xc_down']         = v_xc_inv_down
                        p['MSE_down']          = MSE_down
                        p['E_kin_down']        = E_kin_down
                        p['Psi_down']          = Psi_inv_down

                    sample[sample_idx - 1] =     p.copy()

                    if sys.argv[-1] == "test":
                        print("E_xc_total(kin):   ", p["E_xc_total(kin)"])
                        print("E_xc_total(evals): ", p["E_xc_total(evals)"])
                        print("E_xc_spin(kin):    ", p["E_xc_spin(kin)"])
                        print("E_xc_spin(evals):  ", p["E_xc_spin(evals)"])
                        
                except Exception as e:
                    print("E: ", e)
                    failed_dims.append(p["dim"])
                    failed_nucs.append(nuc_collumn[i])
                    error_table.add_row([str(i), str(failed_dims[-1]), str(failed_nucs[-1]).replace('\n',' ').replace("'", ''), str(e)])
                    dropped += 2*len(alpha) + 1
                    failure  = True

                '''__________PROGRESS UPDATE____________________________________________________________________________________________________________________'''
                end_time = datetime.datetime.now().replace(microsecond=0) - start_time
                print(".........................................................................")
                print(str(round(100*float(idx_newCSV)/float(data_lenght2), 2))+" % completed (" + str(end_time) + ") | "
                      + str(round(100*float(dropped)/float(data_lenght2),2)) +" % failed")

                if counter == 1:
                    break

            last_int_particle_rows = row.copy()

            if i == data_lenght - 1 and alpha[a] == 1.:
                if (None in sample) == False:
                    for s in sample:
                        systems.append(s)


with open(r"info_Errors/Errors_Inversion_" + csv_index + ".txt", "w") as error_file:
    error_file.write("The following Systems have been dropped:\n\n")
    error_file.write(str(error_table))
    error_file.write("\n\nTotal failures: " + str(round(100*float(dropped)/float(data_lenght2), 2)) + " %")
field_names = list(p.keys())

with open(csv_ouput, 'w') as csvfile:  # creates TraingsSet.csv
    writer = csv.DictWriter(csvfile, fieldnames=field_names)
    writer.writeheader()
    writer.writerows(systems)

print("_______________________________\n")
print("Inversion finished!\n")
print(str(round(100*float(dropped)/float(data_lenght), 2)) + " % failures (see info_Errors). \n")
print("_______________________________\n")
print(test_x)
print(test_y)







