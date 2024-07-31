# -*- coding: utf-8 -*-
"""
Short script to get the the absorbed and trasmited power of a multilayer
from a h5 file created by saving a 3D undulator radiation OASYS widget
"""

__author__ = "Juan Reyes-Herrera"
__contact__ = "juan.reyes-herrera@esrf.fr"
__licence__ = "MIT"
__copyright__ = "ESRF - The European Synchrotron, Grenoble, France"
__date__ = "13/04/2023"
__version__ = "0.0.1"

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.constants as codata

### Defining the path to save the figure ###        
#script_dir = os.path.dirname(__file__)
script_dir = os.path.join(os.getcwd())
#images_dir = script_dir
images_dir = os.path.join(script_dir,"")
# font and label size for the plots #
f_size = 12

#TODO implement interpolation for now it work only with equal energy arrays

class MultilayerMirror:
    """
    A Multilayer Mirror gets the incoming flux/power from a 3D file from XOPPY
    undulator radiation and get its reflectivity from an external file, from IDM
    for example. It can create files and images with the absorbed and trasmitted
    3D flux and power density.

    Attributes
    ----------
    source_file : str
        Name of the h5 with the 3D radiation data
    ref_file : str
        Name of the txt multilayer reflectivity file
    angle: float or int
        Mirror grazing angle in deg    

    Methods
    -------
    und_rad_axis(axis)
        Get the values of energies 'e', horizontal 'x' or vertical 'y' axes

    reflec_axis(axis):
        Get an array of the reflectivity axis coordinates

    und_rad_data(slice):
        Get 3D data from undulator radiation, can get sigle value, a slice or full

    plot_ref()
        Plots the multilayer reflectivity

    plot_spectrum(self, type='incoming')    
    
        Plots different flux spectra.
          -'incoming': Flux from the source
          -'absorbed': Flux absobed on the mirror
          -'transmitted': Flux tramitted after the mirror reflection
          -'combined': Comparing incoming flux and mirror reflection file

    mirror_power()
        Calculates the absorbed power and trasmitted by the mirror

        Can plot the 2D power density of:
            - 'inc' incoming
            - 'abs' absorbed
            - 'trans'trasmitted
        Saves the 2D power density in CSV ('csv'), HDF5 ('h5') or the total
        flux density 3D ('h5_3d')      

    Notes
    -----
    - Tested for h5 files created by saving the 3D und_rad from
      Undulator Radiation XOPPY widget

    Examples
    --------
    >>> test_ml_mirror = MultilayerMirror('und_pow.h5', 'reflectivity.txt', 1.5, reflection='horizontal')

    >>> test_ml_mirror.mirror_power(self, kind='trans', plot=True, save_data='h5_3d')
    Plots the absorbed power in the mirror and saves the data as HDF5        
    """

    def __init__(self, source_file, ref_file, angle, reflection = 'vertical'):

        """
        Parameters
        ----------
        source_file : str
            Name of the h5 with the 3D radiation data
        ref_file : str
            Name of the txt multilayer reflectivity file
        angle: float or int
            Mirror grazing angle in degrees
        reflection: str
            Specify if the mirror is deflecting in the plane: 'vertical' or
            'horizontal'
        """ 
        self.source_file = source_file
        self.ref_file = ref_file        

        if type(angle) is float or int:
            self.angle = float(angle)
        else:
            raise RuntimeError("ERROR: please provide a number for mirror angle")
        
        self.reflection = reflection      
   
    def __str__(self):
        return "MultilayerMirror obj"    
    
    def und_rad_axis(self, axis):
        
        """ Get an array of the undulator radiation axis coordinates

        Parameters
        ----------
        axis : str
            axis str enegies('e'), horizontal ('x') or vertical ('y')

        Returns
        -------		
        coord : ndarray
            points at the given axis	
        """
        h5_file = h5py.File(f'{script_dir}/{self.source_file}','r')
        
        if axis == 'e':
            tmp = np.array(h5_file['/XOPPY_RADIATION/Radiation/axis0'])
            coord = np.round(tmp, 1)
        elif axis == 'x':
            coord = np.array(h5_file['/XOPPY_RADIATION/Radiation/axis1'])
        elif axis == 'y':
            coord = np.array(h5_file['/XOPPY_RADIATION/Radiation/axis2'])
        else:
            raise RuntimeError('ERROR: provided axis is not identified')
        h5_file.close()
        return coord
    
    def reflec_axis(self, axis):

        """ Get an array of the reflectivity axis coordinates

        Parameters
        ----------
        axis : str
            axis str energies('e') or reflectivity ('r')

        Returns
        -------		
        coord : ndarray
            points at the given axis	
        """

        df_ref = pd.read_csv(self.ref_file, sep=',|\s+', header=None, comment = ';', engine='python')
        
        if axis == 'e':
            tmp = np.array(df_ref.iloc[:,0])
            coord = np.round(tmp, 1)

        elif axis == 'r':
            coord = np.array(df_ref.iloc[:,1])
        else:
            raise RuntimeError('ERROR: provide axis is not identified')

        return coord


    def und_rad_data(self, slice='full'):
        """ Get the 3D data from the source file

        Parameters
        ----------
        slice : str or list
            full energy range ('full') or index range ([initial, final])

        Returns
        -------
        A 3D matrix with the power density distribution
        data : ndarray
        """
        h5_file = h5py.File(f'{script_dir}/{self.source_file}','r')

        if slice == 'full':
            data = np.array(h5_file['/XOPPY_RADIATION/Radiation/stack_data'])
        elif type(slice) is list:
            tmp = np.array(h5_file['/XOPPY_RADIATION/Radiation/stack_data'])
            data = tmp[slice[0]: slice[-1]]
        else:
            raise RuntimeError('ERROR: slice not reconized')
        h5_file.close() 

        return data

    def plot_spectrum(self, type='incoming'):
        """ Plots a spectrum from a given 3D data
        Parameters
        ----------
        type: str
            initial spectrum ('incoming'), absorbed by the mirror ('absorbed'),
            transmitted after the mirror ('transmitted') or compare the incoming
            spectrum with the multilayer reflectivity ('combined')  

        Returns
        -------
        Spectrum plot
        """        
        e = self.und_rad_axis('e')
        h = self.und_rad_axis('x')
        v = self.und_rad_axis('y')                

        if type == 'incoming':
            data3d = self.und_rad_data()
            f = data3d.sum(axis=2).sum(axis=1)*(h[1]-h[0])*(v[1]-v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(e[1]-e[0]), 1)
            plt.plot(e, f)
            plt.title(f'Total incoming power {tot_pow} W', fontsize = f_size)
            plt.ylabel("Flux [photons/s/0.1%bw]", fontsize=f_size)
            plt.xticks(fontsize=f_size)
            plt.yticks(fontsize=f_size)
            plt.show() 
        elif type == 'absorbed':
            data3d = self.mirror_power(kind='abs')
            f = data3d.sum(axis=2).sum(axis=1)*(h[1]-h[0])*(v[1]-v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(e[1]-e[0]), 1)
            plt.plot(e, f)
            plt.title(f'Total absorbed power {tot_pow} W', fontsize = f_size)
            plt.ylabel("Flux [photons/s/0.1%bw]", fontsize=f_size)
            plt.xticks(fontsize=f_size)
            plt.yticks(fontsize=f_size)
            plt.show() 
        elif type == 'transmitted':
            data3d = self.mirror_power(kind='trans')
            f = data3d.sum(axis=2).sum(axis=1)*(h[1]-h[0])*(v[1]-v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(e[1]-e[0]), 1)
            plt.plot(e, f)
            plt.title(f'Total transmitted power {tot_pow} W', fontsize = f_size)
            plt.xlabel("Photon energy [eV]", fontsize=f_size)
            plt.ylabel("Flux [photons/s/0.1%bw]", fontsize=f_size)
            plt.xticks(fontsize=f_size)
            plt.yticks(fontsize=f_size)
            plt.show() 

        elif type == 'combined':
            data3d = self.und_rad_data()
            f = data3d.sum(axis=2).sum(axis=1)*(h[1]-h[0])*(v[1]-v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(e[1]-e[0]), 1)
            r= self.reflec_axis('r')
            fig, ax1 = plt.subplots()
            ax1.set_xlabel('Photon energy [eV]')
            ax1.set_ylabel("Flux [photons/s/0.1%bw]", color='r', fontsize=f_size)
            ax1.plot(e, f, color='r')
            ax1.tick_params(axis='y', labelcolor='r')
            
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis            
            
            ax2.set_ylabel('Reflectivity [a.u.]', color='b', fontsize=f_size)  # we already handled the x-label with ax1
            ax2.plot(e, r, color='b')
            ax2.tick_params(axis='y', labelcolor='b', labelsize=f_size)
            
            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()
        else:
            raise RuntimeError('type of spectrum not reconized')                         


    def plot_ref(self):
        """ Plots the multilayer reflectivity from the input file """

        plt.ion()
        plt.figure()
        plt.plot(self.reflec_axis('e'), self.reflec_axis('r'))
        plt.xlabel("Photon energy [eV]", fontsize=f_size)
        plt.ylabel("Reflectivity [a.u.]", fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        plt.show()
    
    def mirror_power(self, kind='trans', plot=False, save_data=False):
        """Gets the 3D data weighted using multilayer mirror reflectivity,
        have the option to plot and directly save the power dentisty data
        as CSV or H5File

        Parameters
        ----------
        kind : str
            kind str incoming('inc'), absorbed ('abs') or transmitted ('trans')
        plot: bool
            plot bool True or False
        save_data: str
            save_data str power density CSV ('csv'), power density h5 ('h5')
            or flux density 3D ('h5_3d')

        Returns
        -------
        A 3D matrix with the incoming, transmitted or absorbed data by
        the multilayer        
        """
        e = self.und_rad_axis('e')
        h = self.und_rad_axis('x')        
        v = self.und_rad_axis('y')
        r = self.reflec_axis('r')

        energy_step = e[1] - e[0]

        input_data = self.und_rad_data()

        data3d = np.zeros_like(input_data)

        #calcualtions part

        if kind == 'inc':
            data3d = self.und_rad_data()
            data2d = data3d.sum(axis=0) * (energy_step) * codata.e * 1e3
            f = data3d.sum(axis=2).sum(axis=1)*(h[1] - h[0]) * (v[1] - v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(energy_step), 1)
            title = 'incoming'

        elif kind == 'trans':
            for i in range(len(e)):
                data3d[i] = input_data[i] * r[i]
            data2d = data3d.sum(axis=0) * (energy_step) * codata.e * 1e3
            f = data3d.sum(axis=2).sum(axis=1)*(h[1] - h[0]) * (v[1] - v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(energy_step), 1)
            title = 'transmitted'    
                
        elif kind == 'abs':
            for i in range(len(e)):
                data3d[i] = input_data[i] - input_data[i] * r[i]
            data2d = data3d.sum(axis=0) * (e[1] - e[0]) * codata.e * 1e3
            f = data3d.sum(axis=2).sum(axis=1)*(h[1] - h[0]) * (v[1] - v[0])
            tot_pow = round(f.sum()*1e3*codata.e*(energy_step), 1)
            #projection to the mirror surface
            if self.reflection == 'horizontal':
                h /= np.cos((90 - self.angle) * np.pi / 180)
                along_mirror = h
                width_mirror = v
                data3d *= np.cos((90 - self.angle) * np.pi / 180)
                data2d *= np.cos((90 - self.angle) * np.pi / 180)
                #data3d = np.transpose(data3d, (0, 2, 1))
                data2d = data2d.T 
                       
            elif self.reflection == 'vertical':
                v /= np.cos((90 - self.angle) * np.pi / 180)
                along_mirror = v
                width_mirror = h
                data3d *= np.cos((90 - self.angle) * np.pi / 180)
                data3d = np.transpose(data3d, (0, 2, 1))                
                data2d *= np.cos((90 - self.angle) * np.pi / 180)
                                
            else:
                raise RuntimeError('mirror reflection plane not reconized') 

            
            title = 'absorbed'
                
        else:
            raise RuntimeError('type of calcualtion not reconized') 

        
        #redefining the variables to plot or to save        
        if kind == 'inc' or kind == 'trans':

            x = h
            y = v
            axis_labels = ["Horizontal [mm]", "Vertical [mm]"]

        elif plot and kind == 'abs':

            x = along_mirror
            y = width_mirror
            axis_labels = ["Along the mirror [mm]", "Footprint width [mm]"]
        #plotting block    
        if plot:

            peak_power = round(data2d.max(), 2)
            plt.ion()
            plt.figure()
            plt.pcolormesh(x, y, data2d, cmap=plt.cm.viridis, shading='auto')	
            plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)	
            plt.title("Total {} power = {} W, Power Density Peak = {} W/mm$^2$".format(\
                      title, tot_pow, peak_power), fontsize=f_size)
            
            plt.xlabel(axis_labels[0], fontsize=f_size)
            plt.ylabel(axis_labels[1], fontsize=f_size)
            
            plt.xticks(fontsize=f_size)
            plt.yticks(fontsize=f_size)
            plt.xlim(x[0], x[-1])
            plt.ylim(y[0], y[-1])
            #plt.axis('equal')
            plt.show()
         
        #saving block
        #power density

        if save_data == 'csv':
            #Save data in milimiters
            save_matrix = np.zeros((data2d.shape[0] + 1, data2d.shape[1] + 1))         
            save_matrix[0][1:] += x
            save_matrix[: , 0][1:] += v
            save_matrix[1: , 1: ] += data2d

            np.savetxt(f'{title}_pow_dens.csv', save_matrix, delimiter=',')
            print(f"File written to disk: {title}_pow_dens.csv")
        
        elif save_data == 'h5':           

            h5_name2save = f'ml_{self.reflection}_power_density_{kind}.hdf5'
            h5_f = h5py.File(h5_name2save, 'a')

            nxentry = h5_f.create_group('pow_dens')
            nxentry.attrs['NX_class'] = 'NXentry'
            
            nxdata = nxentry.create_group('pow_dens_{}'.format(kind))
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['signal'] = "pow_dens"
            nxdata.attrs['axes'] = ['x', 'y']
                
            # Image data
            fs = nxdata.create_dataset("pow_dens", data=data2d if (kind=='inc' or kind=='trans') else data2d.T)
            fs.attrs['long_name'] = 'Power density'    
        
            # X axis data
            xs = nxdata.create_dataset('x', data=x)
            xs.attrs['units'] = 'mm'
            xs.attrs['long_name'] = axis_labels[0]    
        
            # Y axis data
            ys = nxdata.create_dataset('y', data=y)
            ys.attrs['units'] = 'mm'
            ys.attrs['long_name'] = axis_labels[1]    
        
            h5_f.close()
            
            print ('File {} save to the disk'.format(h5_name2save))

        # saving flux density
        elif save_data == 'h5_3d':

            h53d_name2save = f'ml_{self.reflection}_flux_density_{kind}.hdf5'
            h53d_f = h5py.File(h53d_name2save, 'a')

            nxentry = h53d_f.create_group('flux_dens')
            nxentry.attrs['NX_class'] = 'NXentry'
            
            nxdata = nxentry.create_group('flux_dens_{}'.format(kind))
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['signal'] = "flux_dens"
            nxdata.attrs['axes'] = ['energy', 'x', 'y']            
                
            # Image data
            fs = nxdata.create_dataset("flux_dens", data=data3d)
            fs.attrs['long_name'] = 'Flux density'

            # photon energies
            es = nxdata.create_dataset('energy', data=e)
            es.attrs['units'] = 'eV'
            es.attrs['long_name'] = 'Photon Energy (eV)'     
        
            # X axis data
            xs = nxdata.create_dataset('x', data=x)
            xs.attrs['units'] = 'mm'
            xs.attrs['long_name'] = axis_labels[0]   
        
            # Y axis data
            ys = nxdata.create_dataset('y', data=y)
            ys.attrs['units'] = 'mm'
            ys.attrs['long_name'] = axis_labels[1]    
        
            h53d_f.close()
            
            print ('File {} save to the disk'.format(h53d_name2save))

        else:
            pass
            #raise RuntimeError('To save the file please provide the save data file type (csv or h5)') 

        return data3d

def add_2d_plot(data2d_list, h, v, save_file=False):

    # just in case it has different energy steps    
    sum_data2d = np.zeros_like(data2d_list[0])    
    for data2d in data2d_list:
        sum_data2d += data2d
    peak_power = round(sum_data2d.max(), 2)
    total_power = round(sum_data2d.sum()*(h[0] - h[1])*(v[0] - v[1]), 1)
    plt.ion()
    plt.figure()
    plt.pcolormesh(h, v, sum_data2d.T, cmap=plt.cm.viridis, shading='auto')	
    plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)	
    plt.title("Total sum power = {} W, Power Density Peak = {} W/mm$^2$".format(\
              total_power, peak_power), fontsize=f_size)
    plt.ylabel("Vertical [mm]", fontsize=f_size)
    plt.xlabel("Horizontal [mm]", fontsize=f_size)
    plt.xticks(fontsize=f_size)
    plt.yticks(fontsize=f_size)
    plt.xlim(h[0], h[-1])
    plt.ylim(v[0], v[-1])
    #plt.axis('equal')
    plt.show()

    #Todo implement save file

if __name__=="__main__":
    pass
    #id29_ml_mirror = MultilayerMirror('11.56keV_und_rad_ESRF_ID29_EBS_CPMU16_4.h5', 'ml_ref_Pb-B4C-C_9_100keV.txt', 1.06, reflection='horizontal')