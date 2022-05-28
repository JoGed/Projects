import numpy as np
import OctCalc
import sys
import os
import csv
from prettytable import PrettyTable
import copy
from tqdm import tqdm

print("Hello from ", os.getcwd())
if len(sys.argv) < 3:
    print("_________HOPPLA!__________________________________________________")
    print("Arguments are missing! \nType 'python DensCreator.py <index> <integer> <integer1-integer2> <string>', \n"
          "where 'index' will be assigned to the DensData_<index>.csv output and 'integer' is the number of external potentials.\n"
          "'integer1-integer2' specifies the range of particles for which the densities will be calculated.\n"
          + "'string' specifies the calculation mode, e.g.: g - groundstate | e6 - excited mode with excited sates < 6")
    print("__________________________________________________________________")
    sys.exit()

mode = sys.argv[4].replace('e', '')
if mode != 'g':
    extraStates = int(sys.argv[4].replace('e', ''))

dim_range   = sys.argv[3].replace('-',' ').split()
low_dim     = int(dim_range[0])
end_dim     = int(dim_range[-1])
'''___________________________________TRAININGSSET OPTIONS______________________________________________________________________________'''
random_pots  = True # if true, random potentials will be created based on params['Z_range'], params['Z_range'] and params['nuc_pos_radius']
dims_to_calc = np.arange(low_dim, end_dim +1) # e.g. [1,2]: Dens of 1 and 2 particle system will be calculated for same potentials
samples      = int(sys.argv[2]) # samples * dims_to_calc == number of (Dens, v_xc)-pairs in the trainings set
csv_index    = int(sys.argv[1]) # index fpr .csv output
'''_______________________PARAMETERS_____________________________________________________________________________________________________'''
params = {
    #______SYSTEM (1D)______________________________________________________________________________________________________________________
    'radius':         11.5,                     # 2 * radius = 'box_length' in 1D
    'dx':             0.1,                       # grid spacing
    'total_charge':   3,                       # total charge of the nuclei
    'nuc_pos_radius': 4.,                      # max_range of positions of nuclei
    'a':              1,                      # Coulomb parameter
    'nucs_array':     np.array([[3] , [0]], dtype=str), # array containg [charge_arr , position_arr] (dtype = str) of each nuclei.
                                                         # Explicit definition is only necessary, if random_pot == False
    #______TO CALCULATE______________________________________________________________________________________________________________________
    'nucs':           None,                     # number of nuclei
    'dim':            None,                     # number of particles (defined by 'dim_to_calc' variable)
    'v_ext':          None,                     # external potential (defined by 'random_pots'. If False: 'nucs_array' will be needed
    'Dens_data_total':      None,               # exact Density (Up + Down) (calculated by octopus + SpinDens3Particle script)
    'Dens_data_up':   None,                     # exact Density_Up   (calculated by octopus + SpinDens3Particle script)
    'Dens_data_down': None,                     # exact Density_Down (calculated by octopus + SpinDens3Particle script)
    'Energy':         None,                     # energy of the groundstate
    'state':          None,                     # state index
    'nspindown':      None,                     # number of spin downs

    'dens3_total_wf_YP1': None,                 # total dens (down + up), projected wf on 1st young diagram as primitive spatial part used
    'dens3_up_wf_YP1':    None,                 # dens_up, projected wf on 1st young diagram as primitive spatial part used
    'dens3_down_wf_YP1':  None,                 # dens_down, projected wf on 1st young diagram as primitive spatial part used

    'dens3_total_wf_YP2': None,                 # total dens (down + up), projected wf on 2nd young diagram as primitive spatial part used
    'dens3_up_wf_YP2':    None,                 # dens_up, projected wf on 2nd young diagram as primitive spatial part used
    'dens3_down_wf_YP2':  None,                 # dens_down, projected wf on 2nd young diagram as primitive spatial part used

    'Energy_2dExState':  None,                 # Energy of 1st excited state for 2 particles (needed for Koopman's Theorem for open Shells)
}
csv_name = 'DensData_' + str(csv_index) + '.csv'  # creates new csv file
params['points']    = np.arange(-params['radius'], params['radius'] + params['dx'], params['dx']) # spatial points
params['N']         = len(params['points'])    # number of spartial points
'''________________________________________________________________________________________________________________________________________'''
field_names    = list(params.keys())

if mode == 'g':
    output_str = r"DensSets/" + csv_name
else:
    output_str = r"ExcitedDensSets/" + csv_name

with open(output_str, 'w') as csvfile: #creates empty csv file
    writer = csv.DictWriter(csvfile, fieldnames=field_names)
    writer.writeheader()

error_table = PrettyTable()
error_table.field_names = ["System", "at dim", "nucs_array", "Error"]
dropped = 0

for s in tqdm(range(samples)):
    if mode == 'g':
        group = (1+dims_to_calc[-1]- dims_to_calc[0]) * [None]
    else:
        group = (1+ extraStates) * [None]

    if random_pots:
        if params['total_charge'] != 3:
            print("Creation of random potentials only exists for 3 particles yet!")
            sys.exit()

        nucs = np.random.randint(1, 4)
        nucs_Z = None
        nucs_Z_int = np.zeros(nucs)

        if nucs == 1:
            nucs_Z_int = np.array([3])
        if nucs == 2:
            nucs_Z_int = np.array([1, 2])
        if nucs == 3:
            nucs_Z_int = np.ones(3, dtype=int)

        nucs_Z   = nucs_Z_int.astype(str)
        nucs_pos = (2 * params['nuc_pos_radius'] * np.random.rand(nucs) - params['nuc_pos_radius']).astype(str)
        params['nucs'] = nucs

    else:
        params['nucs'] = len(params['nucs_array'][0])

    try:

        for dim in dims_to_calc:  # safe for 1 partcile number yet if excited states mode used
            print("dim: ", dim)
            p_sample = copy.deepcopy(params)

            if random_pots:
                p_sample['nucs_array'] = np.array([nucs_Z, nucs_pos])

            p_sample['dim'] = dim

            if mode == 'g':
                extraStates = 2
                if dim == 4:
                    extraStates = 10

                # create inp (& returns v_ext if v_ext_return = True)
                p_sample['v_ext'] = OctCalc.Create_InpFile(p_sample, v_ext_return=True, plot=False,
                                                           extra_states=extraStates)

                # return data calculated by octopus:
                p_sample['Dens_data_total'], p_sample['Dens_data_up'], p_sample['Dens_data_down'], \
                Energy, p_sample['state'], p_sample['nspindown'], \
                p_sample['dens3_total_wf_YP1'], p_sample['dens3_up_wf_YP1'], p_sample['dens3_down_wf_YP1'], \
                p_sample['dens3_total_wf_YP2'], p_sample['dens3_up_wf_YP2'], p_sample['dens3_down_wf_YP2'] \
                    = OctCalc.return_Dens(p_sample, extra_states=extraStates, exc_mode=False)

                p_sample['Energy'] = Energy[0]
                p_sample['Energy_2dExState'] = Energy[1] if p_sample['dim'] == 2 else None

                if len(group) == 1:
                    group[0] = copy.deepcopy(p_sample)
                else:
                    group[dim - 1] = copy.deepcopy(p_sample)

            else:
                raise Exception('SPINDENS CREATION ONLY AVAILABLE ONLY FOR GROUND STATE!')

            if (dim == dims_to_calc[-1]) and (None in group) == False:

                with open(output_str, 'a') as csvfile:  # fills collumns of csv file
                    writer = csv.DictWriter(csvfile, fieldnames=field_names)
                    writer.writerows(group)
            
    except Exception as e:
       dropped += 1
       print("E: " + str(e) + " " + str(round(100*float(dropped)/float(samples), 2)) + " % failures so far.")
       error_table.add_row([str(s), str(p_sample['dim']), str(p_sample['nucs_array']).replace('\n', ' ').replace("'",''), str(e)])
    


with open(r"info_Errors/Errors_DensData_" + str(csv_index) + ".txt", "w") as error_file:
    error_file.write("The following Systems have been dropped:\n\n")
    error_file.write(str(error_table))
    error_file.write("\n\nTotal failures: " + str(round(100*float(dropped)/float(samples), 2)) + " %")
'''______________________________________________________________________________________________________________________________________________'''
print("_______________________________\n")
print("Calculation finished!\n")
print(str(round(100*float(dropped)/float(samples), 2)) + " % failures (see info_Errors). \n")
#df = pd.read_csv(r"DensSets/" + csv_name)
#print(df)
print('-------------------------------------------------------------------')







