#import modules
from __future__ import print_function
import copy
import time
#import gc
from aggregation import riming_runs
from aggregation.riming import gen_polydisperse_monomer
from aggregation import mcs
from aggregation import fallvelocity
import __fallspeed_relations as fallspeedMarkus
from matplotlib import pyplot as plt
import numpy as np
import os
import glob
from scipy import stats

for i_size,sizem in enumerate(np.logspace(np.log10(50.0),np.log10(2000.0),20)):
    size = sizem*1e-6
    # create folder to save the shapefiles properly
    dir_shape = "/work/lvonterz/SSRGA/multimonomer/shape_files/size_%E" % float(size)
    if not os.path.exists(dir_shape):
        os.mkdir(dir_shape)
        print("Directory " , dir_shape ,  " Created ")

        # Defining the properties of the dendrites    
    min_size = 100.0e-6; max_size = 3000.0e-6; separation_size = 1000.0e-6; grid_res = 10.0e-6; rimed = True; align = True; riming_lwp = 0.1; riming_mode="subsequent"; lwp_div=500
    columns = {"psd":"exponential",
            "size":size,
            "min_size":min_size,
            "max_size":separation_size,
            "mono_type":"column",
            "grid_res":grid_res, 
            "rimed":rimed,
            "debug":False}
        # Making a copy of dendrites and turn them to needles (keeping all the parameters unchanged)
    dendrites = {"psd":"exponential",
            "size":size,
            "min_size":separation_size,
            "max_size":max_size,
            "mono_type":"dendrite",
            "grid_res":grid_res, 
            "rimed":rimed,
            "debug":False}
    frac1 = stats.expon.cdf(separation_size,scale=size)-stats.expon.cdf(min_size,scale=size)
    frac2 = stats.expon.cdf(max_size,scale=size)-stats.expon.cdf(separation_size,scale=size) 
    ratio1 = (frac1)/(frac1+frac2)
    ratio2 = (frac2)/(frac1+frac2)
    #set up the distribution
    # Create a new polydisperse generator of dendrites and needles in unequal ratios
    polygen = gen_polydisperse_monomer(monomers=[columns,dendrites], ratios=[ratio1, ratio2])
    print('Polygen done')
    # create properties files, one for each size parameter
    with open(dir_shape+"/size_%E_properties.txt" %float(size),"wb") as f: #
        np.savetxt(f, np.vstack(("#", "size", "min_size", "max_size",  "rimed", "grid_res", "align", "riming_lwp", "riming_mode", "mixingratio columns:dendrites","lwp_div")).T, fmt="%s")
        np.savetxt(f, np.vstack(("#", size, min_size, max_size, rimed, grid_res, align, riming_lwp, riming_mode, ratio1*100, ratio2*100,lwp_div)).T, fmt="%s")      
        np.savetxt(f, np.vstack(("#","N_mono", "mass", "area", "D_max", "vel_HW", "vel_KC")).T, fmt="%s")
        print('with open')
        # Use my polygen as a generator of monomers
        Nmono = 10
        while Nmono<=800:
            t0 = time.time()
            agg_iter = riming_runs.generate_rimed_aggregate(polygen,N=Nmono,align=align,riming_lwp=riming_lwp,riming_mode=riming_mode,lwp_div=lwp_div,iter=True)
            t1 = time.time()
            print('agg_iter done in ', t1-t0)
            rho_i = 917.6
            evol_file = "/rimelwp_%5.1f_size_%E_Nmono_%d_evol.txt" %(float(riming_lwp),float(size),(Nmono))
            evol_file = evol_file.replace(" ", "")
            t0 = time.time()
            with open(dir_shape+evol_file,"wb") as f_iter:
                np.savetxt(f_iter,np.vstack(('mass','area','D_max','vel_HW','vel_KC')).T, fmt='%s')
                mass = []; area = []; D_max = []; vel_HW = []; vel_KC = [] 
                for agg in agg_iter:
                    agg = agg[0]
                    mass.append(rho_i*agg.X.shape[0]*agg.grid_res**3)
                    area.append(agg.vertical_projected_area())
                    D_max.append(mcs.minimum_covering_sphere(agg.X)[1]*2)
                    vel_HW.append(fallvelocity.fall_velocity(agg, method="HW"))
                    vel_KC.append(fallvelocity.fall_velocity(agg, method="KC"))
                np.savetxt(f_iter, np.vstack((mass, area, D_max, vel_HW, vel_KC)).T, fmt="%.6e")
                f_iter.close()
                    
            print('Nmono',Nmono)
            t1 = time.time()  
            print('save evol took', t1-t0)
            quit()
            # save shapefile
            savefile = dir_shape+"/multimonomer_rimelwp_%5.1f_size_%E_nmono_%d.txt" %(float(riming_lwp),float(size),(Nmono))
            savefile = savefile.replace(" ", "")
            np.savetxt(savefile, agg.X, fmt="%.6e")            
            # save properties
            np.savetxt(f, np.vstack((Nmono, mass, area, D_max, vel_HW, vel_KC)).T, fmt="%.6e")
      #      del agg
            #gc.collect
            if Nmono < 100:
                Nmono += 1
            elif Nmono >= 100:
                Nmono += 10
            print(Nmono)
            quit()
