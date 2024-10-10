# -*- coding: utf-8 -*-
"""
Short script to get the the absorbed and trasmited power of a mask from a h5
file created by saving a 2D power density OASYS widget
"""

__author__ = "Juan Reyes-Herrera"
__contact__ = "juan.reyes-herrera@esrf.fr"
__licence__ = "MIT"
__copyright__ = "ESRF - The European Synchrotron, Grenoble, France"
__date__ = "05/04/2022"
__version__ = "0.0.1"

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os

### Defining the path to save the figure ###        
#script_dir = os.path.dirname(__file__)
script_dir = os.path.join(os.getcwd())
#images_dir = script_dir
images_dir = os.path.join(script_dir,"")
# font and label size for the plots #
f_size = 12

class SlitMask:
    """
    A slit mask that can be rectangular or circular that will be applied over the
    2D map h5file created by saving the image from an OASYS XOPPY power density
    widget. It can create images for the transmitted and masked power density.

    Attributes
    ----------
    source_file : str
        Name of the h5 with the 2D power density data
    kind : str
        Kind of the slit, can be 'rectangular' or 'circular'
    hor : float or int
        Horzontal aperture in mm (if circular is the diameter)
    ver : float, int or None
        Vertical aperture in mm (if circular, it is not necessary)

    Methods
    -------
    map_axis(axis)
        Get the values of horizontal 'x' or vertical 'y' axes
    
    map_data()
        Gets the 2D power density data

    abs_plot(save_fig=False)
        Plots the absorbed power density by the mask, option to save the figure

    trans_plot(safe_fig=False)
        Plots the masked power density, have the option to directly save the
        figure

    Notes
    -----
    - Tested for h5 files created by saving the 2D plot of any power density
      plot of XOPPY widgets, please note that the power density has to be in
      W/mm2 units

    Examples
    --------
    >>> mask_test = SlitMask('pow_dens_map_test.h5', 'rectangular', 2.5, 1.5)

    >>> mast_test.abs_plot(save_fig=True)
    Plots the absorbed power in the mask and saves the plot
    Figure abs_pow_dens_map_test.png has been saved in the folder
        
    """

    def __init__(self, source_file, kind, hor, ver=None):

        """
        Parameters
        ----------
        source_file : str
            Name of the h5 source file (e. g. 'pow_dens_map_test.h5')
        kind : str
            Kind of the slit, can be 'rectangular' or 'circular'
        hor : float or int
            Horzontal aperture in mm (if is circular is the diameter)
        ver : float, int or None
            Vertical aperture in mm (if circular, it is not necessary)
        """ 
        self.source_file = source_file
                
        if kind in ['rectangular', 'circular']:
            self.kind = kind
        else:
            raise RuntimeError("ERROR: please indicate 'rectangular' or 'circular'")

        if type(hor) is float or int:
            self.hor = float(hor)
        else:
            raise RuntimeError("ERROR: please provide a number for horizontal aperture")

        if ver is None:
            self.ver = 0.0            
        elif ver is not None:
                self.ver = float(ver)
        else:
            raise RuntimeError("ERROR: please provide a number for vertical aperture")
   
    def __str__(self):
        return "SlitMask obj"

    def __repr__(self):
        
        if self.kind == 'rectangular':
            lbl = "{} mask of {} mm widht and {} mm height".format(\
                                                self.kind, self.hor, self.ver)
        else:
            lbl = "{} mask of {} mm diameter".format(self.kind, self.hor)
            
        return lbl
    
    def map_axis(self, axis):
        
        """ Get a list of the axis coordinates

        Parameters
        ----------
        axis : str
            axis str horizontal ('x') or vertical ('y')

        Returns
        -------		
        coord : ndarray
            points at the given axis	
        """
        h5_file = h5py.File(f'{script_dir}/{self.source_file}','r')
        
        if axis == 'x':
            coord = np.array(h5_file['/entry/data0/x'])
        elif axis == 'y':
            coord = np.array(h5_file['/entry/data0/y'])
        else:
            raise RuntimeError('ERROR: provided axis not identified')
        return coord

    def map_data(self):
        """ Get the 2D data from the source file
        Returns
        -------
        A 2D matrix with the power density distribution
        data : ndarray
        """
        h5_file = h5py.File(f'{script_dir}/{self.source_file}','r')

        data = np.array(h5_file['/entry/data0/image'])

        return data
    
    def abs_plot(self, save_fig=False):
        """Plots the power density removed by the mask, have the option to
        directly save the figure as PNG

        Returns
        -------
        A 2D matrix with the absorbed power density distribution absorbed by
        the mask (optional)

        a_data : ndarray
        """
        x = self.map_axis('x')
        y = self.map_axis('y')        
        a_data = np.copy(self.map_data())

        if self.kind == 'rectangular':
            
            if self.hor/2 <= x[-1]:
                h_min = np.argmin(x <= -self.hor/2)
                h_max = np.argmax(x >= self.hor/2)
            else:
                h_min = 0
                h_max = len(x)

            if self.ver/2 <= y[-1]:
                v_min = np.argmin(y <= -self.ver/2)
                v_max = np.argmax(y >= self.ver/2)
            else:
                v_min = 0
                v_max = len(y)

            a_data[v_min:v_max, h_min:h_max] = 0

        elif self.kind == 'circular':
            
            xr, yr = np.meshgrid(x, y)
            cmask = np.sqrt((xr)**2 + (yr)**2)
            for x1 in range(0, len(x)):
                for y1 in range(0, len(y)):
                    if cmask[x1, y1] <= self.hor/2:
                        cmask[x1, y1] = 0
                    elif cmask[x1, y1] > self.hor/2:
                        cmask[x1, y1] = 1
                    else:
                        pass
                        
            a_data = np.multiply(a_data, cmask)
        
        else:
            raise RuntimeError("ERROR: please indicate 'rectangular' or 'circular'")
        
        plt.ion()
        total_power = round((a_data.sum() * (abs(x[1] - x[0])) * \
                             abs((y[1] - y[0]))) / 1000, 3)
        peak_power = round(a_data.max(), 2)
        plt.ion()
        plt.figure()
        plt.pcolormesh(x, y, a_data, cmap=plt.cm.viridis, shading='auto')	
        plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)	
        plt.title("Total Power = {} kW, Power Density Peak = {} W/mm$^2$".format(\
                                         total_power,peak_power), fontsize=f_size)
        plt.ylabel("Vertical [mm]", fontsize=f_size)
        plt.xlabel("Horizontal [mm]", fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        plt.xlim(x[0], x[-1])
        plt.ylim(y[0], y[-1])
        #plt.axis('equal')   
        
        if save_fig:
            
            save_name = "abs_"+self.source_file[:-3]+".png"
            plt.savefig(images_dir + save_name, dpi=300, bbox_inches='tight',
                        format='png')
            print("Figure {} has been saved in the folder {}".format(\
                  save_name, images_dir))
                  
        plt.show()
        #return a_data

    def trans_plot(self, save_fig=False):
        """Plots the masked power density, have the option to directly save the
        figure as PNG

        Returns
        -------
        A 2D matrix with the transmitted power density distribution after
        applying the mask (optional)

        t_data : ndarray
        """
        x = self.map_axis('x')
        y = self.map_axis('y')
        t_data=np.copy(self.map_data())

        if self.kind == 'rectangular':
            
            if self.hor/2 <= x[-1]:
                t_data[:, :np.argmin(x <= -self.hor/2)] = 0
                t_data[:, np.argmax(x >= self.hor/2):] = 0
            else:
                pass
            if self.ver/2 <= y[-1]:
                t_data[:np.argmin(y <= -self.ver/2) , :] = 0
                t_data[np.argmax(y >= self.ver/2): , :] = 0
            else:
                pass

        elif self.kind == 'circular':

            xr, yr = np.meshgrid(x, y)
            cmask = np.sqrt((xr)**2 + (yr)**2)
            for x1 in range(0, len(x)):
                for y1 in range(0, len(y)):
                    if cmask[x1, y1] > self.hor/2:
                        cmask[x1, y1] = 0
                    elif cmask[x1, y1] <= self.hor/2:
                        cmask[x1, y1] = 1
            t_data = np.multiply(t_data, cmask)

        
        total_power = round((t_data.sum() * (abs(x[1] - x[0])) * \
                             abs((y[1] - y[0]))) / 1000, 3)
        peak_power = round(t_data.max(), 2)
        plt.ion()
        plt.figure()
        plt.pcolormesh(x, y, t_data, cmap=plt.cm.viridis, shading='auto')	
        plt.colorbar().ax.tick_params(axis='y', labelsize=f_size)	
        plt.title("Total Power = {} kW, Power Density Peak = {} W/mm$^2$".format(\
                                        total_power,peak_power), fontsize=f_size)
        plt.ylabel("Vertical [mm]", fontsize=f_size)
        plt.xlabel("Horizontal [mm]", fontsize=f_size)
        plt.xticks(fontsize=f_size)
        plt.yticks(fontsize=f_size)
        plt.xlim(x[0], x[-1])
        plt.ylim(y[0], y[-1])
        #plt.axis('equal')       

        if save_fig:
            
            save_name = "trans_" + self.source_file[:-3]+".png"
            plt.savefig(images_dir + save_name, dpi=300, bbox_inches='tight',
                       format='png')
            print("Figure {} has been saved in the folder {}".format(\
                  save_name, images_dir))        
        plt.show()
        
        #return t_data
