import netCDF4 as nc
import pylab as plt
import numpy as np
import matplotlib.cm as mplcm
import matplotlib.colors as colors

def main(files):

	nc_dataset_list(files)
	
    close_datasets(files)
    return
    
def nc_dataset_list(files):
    files['datasets'] = []
    for p in files['paths']:
        files['datasets'].append(nc.Dataset(p,'r'))

    return

def close_datasets(files):
    for f in files['datasets']:
        f.close()

def test_files(files): 
    files['paths'] = []
    files['names'] = []

    files['paths'].append('./Output.Bomex.41a52/fields/5400/0.nc') # for one core
    files['names'].append('Rico ')
    

if __name__ == '__main__':

    files = {}
    
    test_files(files)
    
    main(files)
