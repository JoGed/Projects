from matplotlib import pyplot as plt
import sys
import numpy as np
import csv
import pandas as pd
from tqdm import tqdm
import matplotlib
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
if len(sys.argv) < 5:
    print("Type: python SimpleRowPlotterLatexSpin.py <path/file> <rows> <d (densities) or t (incl. inversion)> <p (optional)> ")
    sys.exit()
rows       = int(sys.argv[2])
#print(rows)
csv_input  = sys.argv[1]
dim_collumn = np.array(pd.read_csv(csv_input)['dim'])

fig_height  = 4
fig_width   = 7


#import os
#DataSet = pd.read_csv(csv_input).sort_values(by=["nucs_array"])
#DataSet.to_csv(os.getcwd() + "/TestSetsSpecial/DensData_Modified.csv", index=False)
#sys.exit()


def floatArray(string):
    return np.array(string.replace('[', '').replace(']', '').split()).astype(float)

with open(csv_input) as csvfile:
    reader = csv.DictReader(csvfile)

    for i, row in enumerate(reader):
        f, ax = plt.subplots(1, sharex=True, sharey=True)
        f.set_figheight(fig_height)
        f.set_figwidth(fig_width)

        points = floatArray(row['points'])
        v_ext = floatArray(row['v_ext'])
        Dens_data = floatArray(row['Dens_data_total'])
        Dens_data_up = floatArray(row['Dens_data_up'])
        Dens_data_down = floatArray(row['Dens_data_down'])

        if sys.argv[3] != 'd':
            Dens_inv = floatArray(row['Dens_inv'])
            v_H = floatArray(row['v_H'])
            v_xc = floatArray(row['v_xc'])
            v_xc_up = floatArray(row['v_xc_up'])
            v_xc_down = floatArray(row['v_xc_down'])
            MSE = np.format_float_scientific(float(row['MSE']), precision=2)
            maxErr = np.format_float_scientific(float(row['maxDensDiff']), precision=2)
            pl_Vxc_tot, = ax.plot(points, v_xc, color="k", linestyle=(0, (5, 1)))
            pl_Vxc_up, = ax.plot(points, v_xc_up, color="r", linestyle="dashed")
            pl_Vxc_down, = ax.plot(points, v_xc_down, color="g", linestyle="dotted")

        if sys.argv[-1] == 'p':
            ax.plot(points, v_ext)
            plt.ylim(-3, 2)
        else:
            plt.ylim(-1.25, 2)


        D_atBounds = [Dens_data[0], Dens_data[-1]]
        D_atBounds_str = "(" + np.format_float_scientific(float(D_atBounds[0]), precision=1) + ", " \
                         + np.format_float_scientific(float(D_atBounds[-1]), precision=1) + ")"
        info = "fermions : " + row['dim']  # + "\n" + "DaB = " + D_atBounds_str

        if sys.argv[3] != 'd':
            ax.set_ylabel(r'$\rho, v^{xc} \;[a.u.]$', fontsize=20)
        else:
            ax.set_ylabel(r'$\rho \;[a.u.]$', fontsize=20)

        ax.plot(points, Dens_data, label=r'$\rho_{tot}$', color="k")
        ax.plot(points, Dens_data_up, label=r'$\rho_{up}$', color="r", linestyle=(0, (5, 1)))
        ax.plot(points, Dens_data_down, label=r'$\rho_{down}$', color="g", linestyle="dashed")
        
        leg1 = ax.legend(loc='upper left', fontsize=14)
        #leg2 = ax.legend([pl_Vxc_tot, pl_Vxc_up, pl_Vxc_down],
         #                [r'$v^{xc}_{tot}$', r'$v^{xc}_{up}$', r'$v^{xc}_{down}$'],
         #                loc='upper right', fontsize=15)
        ax.add_artist(leg1)


        plt.xlabel(r'$x\;[a.u.]$', fontsize=20)

        ax.tick_params(labelsize=13)
        ax.grid()

        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        ax.text(0.03, 0.14, info, fontsize=12, bbox=props, horizontalalignment='left',
                               verticalalignment='center', transform=ax.transAxes)
        ax.set_aspect(aspect=3.5)
        plt.xlim(-11.5, 11.5)

        #plt.savefig('Pics/SinglePlots_dim='+row['dim'][0]+"_"+row['dim'][-1]+'.pdf', dpi=600)
        #plt.savefig('Pics/SinglePlots_dim=' + row['dim'][0] + "_" + row['dim'][-1] + '.png', dpi=300)
        #plt.savefig('Pics/H2_dis=' + str(int(i)) + '.png', dpi=600)
        #print(row["nucs_array"])
        #print(Dens_data[115])
        #print(Dens_data[0], Dens_data[-1])
        plt.show()


