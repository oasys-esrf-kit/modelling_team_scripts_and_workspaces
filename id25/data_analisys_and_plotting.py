import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import os
from decimal import Decimal



def read_csv(file_name):

    print(file_name)

    df = pd.read_csv(file_name, sep=',|\s+', comment = '#', engine='python')
        
    x = np.array(df.iloc[:,0])
    y = np.array(df.iloc[:,1])

    emittance = float(format(re.findall(r'\d+\.\d+', file_name)[0]))

    return x, y, emittance

def plot_profiles(type='single', file_name='BM_fresne_diff_bm_test_0.1_pm.csv'):

    if type == 'single':
        x, y, emittance = read_csv(file_name)
        plt.plot(x, y, label=f"{emittance} pm")

    elif type == 'all':
        dir_path = r'N:/OASYS/ID25'
        res = []
        for file in os.listdir(dir_path):
            if file.endswith('.csv'):
                res.append(file)
        for csv_file in res:
            x, y, emittance = read_csv(csv_file)
            plt.plot(x, y, label=f"{emittance} pm")
    else:
        raise RuntimeError("ERROR: failure reading CSV files")

    plt.legend()

    plt.xlabel("Vertical position on screen ($\mu$m)",fontsize=12)
    plt.ylabel("Intensity (a.u.)",fontsize=12)

    plt.grid(which='both', axis='y')

    plt.show()


def find_max_min(file_name):

    x, y, emittance = read_csv(file_name)

    mid_indx = int(len(x)/2)

    peak = np.amax(y)
    peak_ind = np.where(y == peak)

    if len(peak_ind[0]) == 1:

        peak_ind = int(peak_ind[0]) 

        if x[peak_ind] < 0:
            delta_indx = mid_indx - peak_ind
            y_min_range = y[peak_ind: mid_indx + delta_indx]
        elif x[peak_ind] > 0:
            delta_indx = peak_ind - mid_indx
            y_min_range = y[mid_indx - delta_indx: peak_ind]
        else:
            raise RuntimeError("ERROR: Something failed while allocating minimum and maximum indexes")
        
    elif len(peak_ind[0]) == 2:

        indxs = peak_ind[0]

        peak_ind_left = int(indxs[0])
        peak_ind_rigth = int(indxs[1])

        y_min_range = y[peak_ind_left: peak_ind_rigth]

    else:
        raise RuntimeError("ERROR: More than 2 maximum where found")

    valley = np.amin(y_min_range)

    return peak, valley, emittance


def plot_visibility(save_file=False):
    
    dir_path = r'N:/OASYS/ID25'
    res = []
    emittances = []
    visibilities = []
    for file in os.listdir(dir_path):
        if file.endswith('.csv'):
            res.append(file)
    for csv_file in res:
        peak, valley, emittance = find_max_min(csv_file)
        visibility = (peak - valley) / (peak + valley)
        emittances.append(emittance)
        visibilities.append(visibility)

    plt.plot(emittances, visibilities, 'bo-')
    plt.xlabel("Vertical emittance (pm rad)",fontsize=12)
    plt.ylabel("Visibility of center dip (a.u.)",fontsize=12)

    plt.xscale('log')
    plt.xlim(0.095, 10)
    plt.grid(which='both', axis='both')

    plt.show() 

    if save_file == True:
        #TODO saving a file #
        pass   

if __name__=="__main__":
    pass


    
    #x1, y1, label1 = read_csv('BM_fresnel_diff_bm_test_0.1_pm.csv')
    #x2, y2, label2 = read_csv('BM_fresnel_diff_bm_test_7.0_pm.csv')
#
    #plt.plot(x1, y1, label=label1)
    #plt.plot(x2, y2, label=label2)

    