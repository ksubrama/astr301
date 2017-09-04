"""Create a copy of PS1problem3starter.py as PS1problem3_yourname.py 
so that you can write python code and answers to the questions as comments in your 
own file.  Feel free to work together, but turn in your own work 
(2 plots + problem1_yourname.py with the edits you make)."""


"""This python file has all of the commands necessary for problem 1.
 Questions are written into the work sheet and are highlighted in green rather
 than the usual grey so you will be sure to notice and answer them as comments  
 in the program you submit. """

import numpy as np
import matplotlib.pyplot as plt

# Part 1 - Read in the file, breaking the data into individual arrays. #

filename = r'/INSERT_PATH_HERE/massdecompdata.txt'
data = np.loadtxt(filename) 
# The columns are radius (kpc), stellar mass (m_sun), 
# neutral gas mass (m_sun) and velocity (km/s)

# Below we make the arrays for the data set

# data
rad = data[:,0] # Distance in kiloparsecs (kpc)
mstars = data[:,1] # Stellar mass in solar masses (m_sun)
mgas = data[:,2] # Neutral gas mass, in solar masses (m_sun)
vel = data[:,3] # Velocity in kilometers/second (km/s)

#//////////////////////////////////////////////////////////////////////////////////////

# Part 2 - Plot stellar mass and HI mass on the same plot using stars and dots, respectively. 
# You will need to plot them in log units!!!

plt.figure(1)

plt.plot(rad,np.log10(mgas),'b.',markersize=3)
plt.plot(rad,np.log10(mstars),'r*',markersize=7)
plt.title('Plot of Neutral Gas Mass and Stellar Mass versus Radius')
plt.xlabel('Radius (kpc)')
plt.ylabel('Log Mass (m_sun)')
ymin,ymax = plt.ylim() # get y axis limits
plt.ylim(0,ymax) # fix lower end of y axis to start at zero

# This plot shows the incremental stellar and gas mass. To understand what that means,
# refer back to the picture of the galaxy in the handout. On the picture of the galaxy
# you drew a ring (annulus) around the galaxy at 10kpc. This ring contains one increment
# of mass.

"""!!QUESTION!! How would you describe the radial distributions of stellar and 
neutral gas mass?"""

#///////////////////////////////////////////////////////////////////////////////////////////

# Part 3 - Plot velocity (km/s) versus distance (kpc), in other words, the 
# rotation curve of the galaxy.

plt.figure(2)

plt.plot(rad,vel,'g.',markersize = 4)
plt.title('Galaxy Rotation Curve')
plt.xlabel('Radius (kpc)')
plt.ylabel('Velocity (km/s)')
ymin,ymax = plt.ylim() # get y axis limits
plt.ylim(0,ymax) # fix lower end of y axis to start at zero

#/////////////////////////////////////////////////////////////////////////////////////////

# Part 4 - Compute the velocity that is expected based on the Motions Find Mass 
# equation at the edge of the visible galaxy (about rad=10kpc). Call this vel10.
# Recall that G=4.28X10-6 kpc(km/s)^2/Msun. Note that strictly speaking the Motions 
# Find Mass equation applies to a spherical mass distribution, but the error in 
# approximating this galaxy as spherical is small. Also remember that (as you figured
# out in Part 2) the radial gas and stellar mass distributions are incremental,
# so you will need to use the numpy sum function to obtain the enclosed mass.

Grav = 4.28e-6 # define G in units of kpc*km^2/(M_sun*s^2)

# Compute the velocity at rad=10 kpc

selrad10=np.where(rad < 10)
masssum10=np.sum(mstars[selrad10])+np.sum(mgas[selrad10])

# Find the expression for the velocity expected at 10kpc using the Motions Find Mass 
# equation and write the Python expression for vel10 below

#vel10 = ????

# Repeat the calculation of the velocity at rad=20kpc filling in your own python
# expressions below

#selrad20 = ???
#masssum20= ???
#vel20 = ???

# Now compare with the galaxy rotation curve by plotting vel10 and vel20 on 
# your rotation curve plot.

plt.figure(2)

#plt.plot(10,vel10,'xm',markersize = 7)
#plt.plot(20,vel20,'xm',markersize = 7)
"""!!Question!! How do vel10 and vel20 compare to the rotation curve at their 
respective radii?"""

#///////////////////////////////////////////////////////////////////////////////////////////

# Part 5 - Calculate the total mass of the galaxy enclosed within each radius 
# using the Motions Find Mass equation.

# Here please note the words "total" and "enclosed" refer to two different things.
# The TOTAL mass of the galaxy comes from the many components of the galaxy: stars,
# HI gas, and dark matter. The word TOTAL used here refers to the addition of all
# of these components. 
# The word "enclosed" refers to summing all the annuli interior to a given radius.
# We don't use the word "total" here, because we could add up all the enclosed 
# stellar mass without including the gas mass (for example).  
# So far we have plotted only incremental radial distributions of mass as in part 1
# (i.e., mass within individual annuli). The radial distribution of enclosed mass 
# is referred to as the CUMULATIVE mass distribution. 
# Again, it is key to understand that the words TOTAL and CUMULATIVE used in 
# this situation have different meanings. TOTAL refers to the addition of all the 
# galaxy mass components, while cumulative refers to sum of the masses in all the
# annuli up to a certain radius. So TOTAL mass may be plotted as either an incremental 
# or a cumulative distribution, and CUMULATIVE mass distributions can be plotted for 
# stellar, gas, or total mass.

totmassEnc = (vel**2)*rad/Grav

"""!!QUESTIONS!! Is totmassEnc a cumulative or incremental radial distribution of the
total mass of the galaxy? (Hint: think about what the mass in the motions find mass
equation is.) Why does it NOT make sense to include a plot of totmassEnc in plot window
1 along with our stellar and gas mass distributions? What kind of distribution 
(cumulative or incremental) of total mass should we plot along with the stellar and 
gas mass distributions?"""


# We can find the largest radius we have data for using the "len" function.  The len fuction
# returns the number of elements in the array.  The last element in the array will have an 
#index of len(array)-1.

maxR = rad[len(rad)-1] #maxR = largest radius in data set = 21.978 kpc 

# Find the total mass, the stellar mass, and the gas mass enclosed in this radius. 

#totalmass= ???
#starmass= ???
#gasmass = ???

"""!!QUESTIONS!! What fraction of the total mass enclosed within 22 kpc is dark matter?  
Make your python script display this result"""

#//////////////////////////////////////////////////////////////////////////////////////////////////

# Part 6 - Calculate the incremental radial distribution of total mass so that we can
# include it with the other incremental mass distributions (stellar and gas) in plot
# window 1. Use the commands below to calculate the incremental mass. 

#incrementaltotalmass = np.zeros(len(totmassEnc)) # step 1
#incrementaltotalmass[0] = totmassEnc[0] # step 2
#offsetindex = np.arange(1,len(totmassEnc)) #step 3
#incrementaltotalmass[offsetindex] = totmassEnc[offsetindex] - totmassEnc[offsetindex - 1] #step 4

"""!!QUESTIONS!! Why is the first element of incrementaltotalmass the enclosed mass?
Explain in words how we compute the incrementaltotalmass in the last step."""

plt.figure(1)
#plt.plot(rad,np.log10(incrementaltotalmass),'gx',markersize=5)
#plt.ylim(??,??) 


#////////////////////////////////////////////////////////////////////////////////////

# Part 7 - Calculate the radial dark matter distribution from the incrementaltotalmass,
# mstars, and mgas arrays that we have made in this Python script. The dark matter
# is the leftover mass when you subtract stellar and gas mass from total mass.

# Write out how to compute the incremental dark matter distribution in the Python code below.
# Be sure to check that the dark matter mass is always postive.  If the sum of stellar and gas
# mass exceeds the total mass inferred from the rotation speed, then we should conclude that
# there is no evidence for dark matter in that annulus.  

#darkmatter = ??

# Plot the dark matter distribution in plot window 1 with the other incremental
# mass distributions.

plt.figure(1)
#plt.plot(rad,np.log10(darkmatter),'m+',markersize=10)

"""!!QUESTIONS!! How does the radial distribution of dark matter compare to the radial 
distribution of stars and gas? What questions does this raise to your mind?"""

