import pandas as pd
import numpy as np
import pickle
import gzip
from argparse import ArgumentParser
import torch
import random
import sys

'''
Code to transform the trainings set (as .cvs file) into pickeled
dictionary containing pytorch tensors. The .gz can eventually be used as a
data set file for the training & testing (via train.py and test.py)
'''

def str_arrToNumpy(self):
    return np.fromstring(self[1:-1].replace("\n", ""), sep=' ').astype(np.float)

'''______________ARGUMENTS___________________________________________________________________________________________'''
parser = ArgumentParser()
parser.add_argument('--name',        type=str, required=True, help="Dictionary will be created named 'DictTorch_$name$'")
parser.add_argument('--data_dir',        type=str, required=True, help="Path to the csv data file")
parser.add_argument('--shuffle',     type=int, default=0,     help="If True, order of external potentials will be shuffled")
parser.add_argument('--NucsNumCorr', type=int, default=0,     help="If True, number of nuclei will be equalized")
args = parser.parse_args()

'''______________IMPORT DATA_________________________________________________________________________________________'''
def dim_info(Set):
    min_dim = np.amin(np.array(Set['dim']))
    max_dim = np.amax(np.array(Set['dim']))
    max_dim_1stidx = np.where(np.array(Set['dim']) == max_dim)[0][0]
    dim_arr = np.array(Set['dim'])[0:max_dim_1stidx + 1]
    dims_length = len(dim_arr)
    return min_dim, max_dim, dim_arr, dims_length

with open(args.data_dir) as data_file:
    DataSet = pd.read_csv(args.data_dir)
    points = str_arrToNumpy(DataSet['points'][0])
    min_dim, max_dim, dim_arr, dims_length = dim_info(DataSet)

'''______________DATA INFO___________________________________________________________________________________________'''
def info(data, status):
    print("________________" + str(status) +"______________________")
    for d in dim_arr:
        print("dim == " + str(d) + ": " + str(int(list(data['dim']).count(d)) ))
    print("_______________________________________________")
    for i in range(np.amax(np.array(data['nucs']))):
        print("nucs == " + str(i+1) + ": " + str(int(list(data['nucs']).count(i+1) / dims_length)))
    print("_______________________________________________")

    '''______________REMOVING________________________________________________________________________________________'''
    print("Data set shape:", data.shape)
    print("Dimensions:    ", dim_arr)
    print("_______________________________________________")

info(DataSet, "DATA_INFO")

'''___________________________TRANSFORMING TO DICT_TORCH______________________________________________________________'''

def remove_dim(Set, dimension):
    return Set.drop(Set[Set["dim"] == float(dimension)].index).reset_index(drop=True)

def shuffle_Set(Set):
    groups = [Set for _, Set in Set.groupby('nucs_array')]
    random.shuffle(groups)
    return pd.concat(groups).reset_index(drop=True)

def NucsNumCorr_Set(Set):
    groups_NucsNumbers_pre = [Set for _, Set in Set.groupby('nucs')]
    nucs_len_pre = [len(groups_NucsNumbers_pre[i]) for i in range(3)]
    groups_NucsNumbers_post = [groups_NucsNumbers_pre[i][nucs_len_pre[i]-np.amin(nucs_len_pre)::] for i in range(3)]
    print("NucsNumCorr: ", nucs_len_pre, " ---> ", [len(groups_NucsNumbers_post[i]) for i in range(3)], " (Sum = " + str(3*np.amin(nucs_len_pre))+")")
    Set_new = pd.concat(groups_NucsNumbers_post).reset_index(drop=True)
    return shuffle_Set(Set_new) if args.shuffle else Set_new

if args.shuffle:
    DataSet = shuffle_Set(DataSet)
    min_dim, max_dim, dim_arr, dims_length = dim_info(DataSet)
    print("DataSet shuffled!")

if args.NucsNumCorr:
    DataSet = NucsNumCorr_Set(DataSet)
    min_dim, max_dim, dim_arr, dims_length = dim_info(DataSet)
    print("Nucs number equaled!")

info(DataSet, "NEW DATA_INFO")

#print(DataSet[['dim', 'evals']][np.where(np.abs(1.-dim_arr)<1e-9)[0][0]::len(dim_arr)])
#print(DataSet[['dim', 'evals']][np.where(np.abs(2.-dim_arr)<1e-9)[0][0]::len(dim_arr)])

#___CONVERTING TO PYTORCH TENSOR____________________________________________________________________________________
points          = torch.from_numpy(np.stack(DataSet["points"].apply(str_arrToNumpy).to_numpy())).double()
Dens_total      = torch.from_numpy(np.stack(DataSet["Dens_data_total"].apply(str_arrToNumpy).to_numpy())).double()
v_xc_total      = torch.from_numpy(np.stack(DataSet['v_xc'].apply(str_arrToNumpy).to_numpy())).double()
Dens_up         = torch.from_numpy(np.stack(DataSet["Dens_data_up"].apply(str_arrToNumpy).to_numpy())).double()
v_xc_up         = torch.from_numpy(np.stack(DataSet['v_xc_up'].apply(str_arrToNumpy).to_numpy())).double()
Dens_down       = torch.from_numpy(np.stack(DataSet["Dens_data_down"].apply(str_arrToNumpy).to_numpy())).double()
v_xc_down       = torch.from_numpy(np.stack(DataSet['v_xc_down'].apply(str_arrToNumpy).to_numpy())).double()
v_ext           = torch.from_numpy(np.stack(DataSet["v_ext"].apply(str_arrToNumpy).to_numpy())).double()
E_xc_total      = torch.from_numpy(np.stack(DataSet["E_xc_total(kin)"].to_numpy())).double()
E_xc_spin       = torch.from_numpy(np.stack(DataSet["E_xc_spin(kin)"].to_numpy())).double()
E_tot           = torch.from_numpy(np.stack(DataSet['Energy'].to_numpy())).double()
evals_Int1      = torch.from_numpy(np.stack(DataSet["evals"][np.where(np.abs(1.-dim_arr)<1e-9)[0][0]::len(dim_arr)]
                                        .apply(str_arrToNumpy).to_numpy())).repeat(1, len(dim_arr)).view(E_tot.shape[0],-1).double()
evals_Int2      = torch.from_numpy(np.stack(DataSet["evals"][np.where(np.abs(2.-dim_arr)<1e-9)[0][0]::len(dim_arr)]
                                        .apply(str_arrToNumpy).to_numpy())).repeat(1, len(dim_arr)).view(E_tot.shape[0],-1).double()
dim             = torch.from_numpy(np.stack(DataSet['dim'].to_numpy())).double()
E_tot_Triplet   = torch.stack([E_tot[0::len(dim_arr)], E_tot[int(len(dim_arr) / 2)::len(dim_arr)], E_tot[(len(dim_arr)-1)::len(dim_arr)]], dim=1)
E_tot_Triplet   = E_tot_Triplet.repeat(1, len(dim_arr)).view(E_tot.shape[0],-1)

DictTorch = {
    "points":        points,
    "Dens_total":    Dens_total,
    "Dens_up":       Dens_up,
    "Dens_down":     Dens_down,
    "v_ext":         v_ext,
    "E_xc_NoSpin":   E_xc_total,
    "E_xc_Spin":     E_xc_spin,
    "v_xc_NoSpin":   v_xc_total,
    "v_xc_up":       v_xc_up,
    "v_xc_down":     v_xc_down,
    "E_tot":         E_tot,
    "evals_Int1":    evals_Int1,
    "evals_Int2":    evals_Int2,
    "E_tot_Triplet": E_tot_Triplet,
    "dim":           dim,
}

f = gzip.open(r"DictTorch_" + str(args.name) + ".gz", 'wb')
pickle.dump(DictTorch, f)
f.close()
print("File 'DictTorch_" + str(args.name) + ".gz'" + " successfully created!")
