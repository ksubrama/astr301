import numpy as np
import matplotlib.pyplot as plt

#read in the data columns: vel contains the radial velocites of all the stars, 
#and prob contains the probability that each star is a member of the dSph.

target_id, ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec, vmag, vicol, vel, vel_errs, prob = np.loadtxt("INSERT_PATH_HERE/dSph_Carina.dat",unpack=True)

#make a histogram of the velocity values with 50 bins
plt.figure(1)
plt.clf()
plt.hist(vel,50)
plt.xlabel("Radial Velocity [km/s]")
plt.ylabel("Number of stars")
plt.title("Carina: All Stars")
