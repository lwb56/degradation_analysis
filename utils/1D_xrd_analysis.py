import sys, os
import glob
import matplotlib.pylab as plt
import numpy as np
import scipy as sp
import pandas as pd



class xrd_analysis():
    # this class is to analyze the 1D xrd pattern match it with cif files

    # ALL analysis and plotting will be in q values

    def __init__(self, data_file, cif_dir, giwaxs_data = False):
        import matplotlib as mpl
        mpl.rcParams["axes.linewidth"] = 5
        mpl.rcParams["axes.labelsize"] = 40
        mpl.rcParams["axes.labelweight"] = 'bold'
        mpl.rcParams['xtick.major.size'] = 14
        mpl.rcParams['xtick.minor.size'] = 10
        mpl.rcParams['ytick.major.size'] = 14
        mpl.rcParams['ytick.minor.size'] = 10
        mpl.rcParams['xtick.major.width'] = 5
        mpl.rcParams['ytick.major.width'] = 5
        mpl.rcParams['xtick.minor.width'] = 4
        mpl.rcParams['ytick.minor.width'] = 4
        mpl.rcParams['font.size'] = 38
        mpl.rcParams["font.family"] = 'Arial'
        mpl.rcParams["font.weight"] = 'normal'
        mpl.rcParams['text.usetex'] = False


        # get cif file directory with all the simulated 1D xrd (xy files)
        self.cif_dir = cif_dir
        self.cif_files = glob.glob(os.path.join(cif_dir, '*.xy'))
        self.data_name = data_file
        self.giwaxs_data = giwaxs_data


    def convert_ras_data(self):
        # this function returns 2 values, the q and the respective intensity

        file1 = open(self.data_name, 'r')
        Lines = file1.readlines()
        ct = 0
        for line in Lines:
            ct+=1
            if '*RAS_INT_START' in line:
                index_begin = ct+1
            if '*RAS_INT_END' in line:
                index_end = ct-1

        theta2 = []
        intensity = []

        count = 0
        for line in Lines:
            count+=1
            if (count >= index_begin) & (count <= index_end):
                val = line.split(' ')
                theta2.append(float(val[0]))
                intensity.append(float(val[1]))

        theta2 = np.array(theta2)
        intensity = np.array(intensity)

        # convert to q values

        theta2q = lambda theta2: 4*np.pi/(1.5604)*np.sin(np.radians(theta2/2))

        return theta2q(theta2), intensity


    def plot_cif(self, which_cif = 'BA_n2'):

        if type(which_cif) != str:
            print('wrong plotting function')
            return None
        
        # finding the cif files in the set
        for i in self.cif_files:
            if which_cif in i:
                cif = i

        # loading the cif files
        cif_data = np.loadtxt(cif)

        cif_theta = cif_data[:,0]
        cif_intensity = cif_data[:,1]

        theta2q = lambda theta2: 4*np.pi/(.9184)*np.sin(np.radians(theta2/2))
        
        q, intensity = self.convert_ras_data()

        cif_q = theta2q(cif_theta)

        indx = np.argwhere(cif_q > q[-1])
        indx = indx[0][0]

        plt.figure(figsize = (20,10))
        plt.plot(cif_q[0:indx],cif_intensity[0:indx]/100, linewidth  = 4, label = which_cif)
        plt.plot(q,intensity/np.max(intensity), linewidth  = 4)
        plt.yticks([])
        plt.xlabel("$q (\AA^{-1})$")
        plt.ylabel('Intensity')
        plt.xticks(np.arange(.2,1.5,.2))
        plt.legend(fontsize = 20)

    

    def plot_multiple_cif(self, which_cif = []):

        theta2q = lambda theta2: 4*np.pi/(.9184)*np.sin(np.radians(theta2/2))
        
        q, intensity = self.convert_ras_data()

        plt.figure(figsize = (20,10))
        plt.plot(q,intensity/np.max(intensity), linewidth  = 4)

        
        for cif_file in which_cif:
            # finding the cif files in the set
            cif_found = False
            for i in self.cif_files:
                if cif_file in i:
                    cif = i
                    cif_found = True

            if cif_found == True:
                # loading the cif files
                cif_data = np.loadtxt(cif)

                cif_theta = cif_data[:,0]
                cif_intensity = cif_data[:,1]
                cif_q = theta2q(cif_theta)
                indx = np.argwhere(cif_q > q[-1])
                indx = indx[0][0]

                plt.plot(cif_q[0:indx],cif_intensity[0:indx]/100, linewidth  = 4, label = cif_file)
            else:
                print('cif file not found/matched')

        plt.yticks([])
        plt.xlabel("$q (\AA^{-1})$")
        plt.ylabel('Intensity')
        plt.xticks(np.arange(.2,1.5,.2))
        plt.legend(fontsize = 20)

    
    def plot_multiple_data(self,filenames_set, dat_labels, log_y = False):

        # plot the first dataset
        q, intensity = self.convert_ras_data()
        intensity = intensity - np.min(intensity) #normalization

        plt.figure(figsize = (20,10))
        
        if log_y == True:
            plt.plot(q,np.log10(intensity/np.max(intensity)), linewidth  = 4, label = dat_labels[0])
        else:
            plt.plot(q,(intensity/np.max(intensity)), linewidth  = 4, label = dat_labels[0])



        for indx,file in enumerate(filenames_set):
            self.data_name = file
            q, intensity = self.convert_ras_data()
            intensity = intensity - np.min(intensity) #normalization
            
            if log_y == True:
                plt.plot(q,np.log10(intensity/np.max(intensity)), linewidth  = 4, label = dat_labels[indx+1])
            else:
                plt.plot(q,(intensity/np.max(intensity)), linewidth  = 4, label = dat_labels[indx+1])

        plt.yticks([])
        plt.xlabel("$q (\AA^{-1})$")
        plt.ylabel('Intensity')
        plt.xticks(np.arange(.2,1.5,.2))
        plt.legend(fontsize = 20)
        
        return None

            

    
