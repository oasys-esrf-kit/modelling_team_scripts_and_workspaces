
import matplotlib.pyplot as plt
import pandas as pd

def plot_comparison(csv_files):
    """ Plot comparison between the different insertion devices, files is a list """
    f_size = 12

    for file_csv in csv_files:
        if 'CPMU19_2.5m' in file_csv:            
            label = 'CPMU19_2.5m'            
            
        elif 'CPMU20.5_2.5m' in file_csv:            
            label = 'CPMU20.5_2.5m'

        elif 'CPMU20.5_4m' in file_csv:            
            label='CPMU20.5_4m'

        else:
           raise RuntimeError("ERROR: Unidentified insertion device")
        
        df = pd.read_csv(file_csv, sep=',', comment = '#', engine='python')
        plt.plot(df.iloc[:, 0], df.iloc[:, 1], label=label)
    
    plt.xlabel("Photon energy (eV)", fontsize= f_size)    
    plt.ylabel(f"Total power (W)", fontsize= f_size)

    plt.xlim(4e3, 20.1e3)
    #plt.ylim(1e12, 5e15)
    plt.xticks(fontsize= f_size)
    plt.yticks(fontsize= f_size)
    #plt.yscale("log")
    plt.grid(which='both', axis='y')
    #plt.title("EBS ID06 CMPU18, On Resonance, $\epsilon$=0, $\delta$=0", fontsize=f_size)
    #plt.title("EBS ID06 CMPU18, On Resonance, $\delta$=0", fontsize=f_size)
    #plt.title("EBS ID20 Power on DCM for some future sources", fontsize=f_size)

    plt.legend(fontsize = f_size)
       
    plt.show()    
        

if __name__=="__main__":
    plot_comparison(['id20_pow_on_dcm_CPMU19_2.5m.csv', 'id20_pow_on_dcm_CPMU20.5_2.5m.csv', 'id20_pow_on_dcm_CPMU20.5_4m.csv'])