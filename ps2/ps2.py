"""This python file has all of the commands necessary for problem 2.
 Questions are written into the work sheet and are highlighted in green rather
 than the usual grey so you will be sure to notice and answer them as comments
 in the program you submit. """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import sys


GRAV = 4.28e-6  # define G in units of kpc*km^2/(M_sun*s^2)


def load_file(filename):
  """ Part 1 - Read in the file, breaking the data into individual arrays.

  returns (np.array, np.array, np.array, np.array)
  """

  data = np.loadtxt(filename)

  # Data describes breakdown of galactic mass distribution in terms of annular rings
  # that are ~0.8kpc thick (as seen by us?).
  rad = data[:, 0]     # Radius of body in kiloparsecs (kpc) from... inferred center?
  mstars = data[:, 1]  # Stellar mass in annulus, in solar masses (m_sun)
  mgas = data[:, 2]    # Neutral gas mass in annulus, in solar masses (m_sun)
  vel = data[:, 3]     # Velocity (tangential? how?) in kilometers/second (km/s)

  return (rad, mstars, mgas, vel)


def plot_incremental_mass(ax, rad, mstars, mgas):
  """ Part 2 - Plot stellar mass and HI mass on the same plot using stars and dots, respectively.

  You will need to plot them in log units!!!
  """
  ax.plot(rad, np.log10(mgas), 'b.', markersize=4, label='Gas mass')
  ax.plot(rad, np.log10(mstars), 'r*', markersize=4, label='Star mass')
  ax.set_title('Incremental Mass versus Radius')
  ax.set_xlabel('Radius (kpc)')
  ax.set_ylabel('Log Mass (m_sun)')

  # This plot shows the incremental stellar and gas mass. To understand what that means,
  # refer back to the picture of the galaxy in the handout. On the picture of the galaxy
  # you drew a ring (annulus) around the galaxy at 10kpc. This ring contains one increment
  # of mass.

  """!!QUESTION!! How would you describe the radial distributions of stellar and
  neutral gas mass?

  It looks bizarre...  This is a log-linear plot.  The gas incremental mass is rising
  but not as r^2 which is what I'd expect from uniform distribution.  Instead, the gas
  distribution flattens out after 5Mpc.  The stellar mass at that radius also drops - but
  it drops exponentially (because I see a straight line of stars).  What in the world is
  going on?

  I don't know what Newtonian mechanics predicts.  It just seems weird that galaxies seem
  to have these flat curve for incremental masses for every section of annulus.  That would
  mean that the density is going as 1/r^2...  for a sphere...  which this isn't... oh boy.
  """


def plot_rotation_curve(ax, rad, vel):
  """ Part 3 - Plot velocity (km/s) versus distance (kpc), in other words, the
  rotation curve of the galaxy.
  """
  ax.plot(rad, vel, 'g.', markersize=4, label='Observed')
  ax.set_title('Galaxy Rotation Curve')
  ax.set_xlabel('Radius (kpc)')
  ax.set_ylabel('Velocity (km/s)')


def plot_expected_velocity(ax, rad, mstars, mgas, vel):
  """ Part 4 - Compute the velocity that is expected based on the Motions Find Mass
  equation at the edge of the visible galaxy (about rad=10kpc). Call this vel10.
  Recall that G=4.28X10-6 kpc(km/s)^2/Msun. Note that strictly speaking the Motions
  Find Mass equation applies to a spherical mass distribution, but the error in
  approximating this galaxy as spherical is small. Also remember that (as you figured
  out in Part 2) the radial gas and stellar mass distributions are incremental,
  so you will need to use the numpy sum function to obtain the enclosed mass.
  """

  # I'm grabbing elements from this list left and right.  Just assume that we have atleast 2
  # rows in our dataset so that we can compute averages.
  assert len(rad) >= 2

  bounds = (rad[:-1] + rad[1:]) / 2
  # To be nice, let the last bound stretch past the last mean annular radius by the last difference
  # in annular radial difference.  It's crude but whatever...
  bounds = np.append(bounds, rad[-1] + (rad[-1] - rad[-2]) / 2)
  cum_mass = np.cumsum(mstars + mgas)
  pred_v = (cum_mass * GRAV / rad) ** 0.5
  ax.plot(bounds, pred_v, 'xm', markersize=4, label='Expected (without DM)')
  ax.legend()

  """!!Question!! How do vel10 and vel20 compare to the rotation curve at their
  respective radii?"""
  """ I see what you did thar!  10 and 20 kpc conveniently lined up just in between the 9.7 and 10.4
  bins and the 19.7 and the 20.4 bins.  Those bins are ~0.8 kpc apart.  So it makes sense to
  ask for the total mass at those radii without having to worry about how to distribute the mass
  within the annunular extent of a particular bin.  I instead have plotted for you the entire
  predicted rotation curve but I picked points inbetween each bins radial center.  I have a
  suspicion that the arithmetic mean is probably not the correct label here - it should be a
  better function that involves the fact that the volume in a shell of thickness dr goes as r^2.
  But this feels good enough for now.  You should find two entries close to 10 kpc and 20 kpc.
  """


def plot_cumulative_mass(ax, rad, mstars, mgas, vel):
  """ I wanted to do this - this shows the cumulative mass distribution from each of these three
  sources so that I can see what the heck is going on."""
  mstars_tot = np.cumsum(mstars)
  mgas_tot = np.cumsum(mgas)
  mass_tot = (vel ** 2) * rad / GRAV
  # Truncate the actual dark matter density to a small but positive value in order to not
  # confuse our log-lin plotter.
  darkmatter_tot = np.clip(mass_tot - mgas_tot - mstars_tot, a_min=1e-20, a_max=None)

  ax.plot(rad, np.log10(mass_tot), 'gx', markersize=4, label='Total inferrred mass')
  ax.plot(rad, np.log10(mstars_tot), 'r*', markersize=4, label='Star mass')
  ax.plot(rad, np.log10(mgas_tot), 'b.', markersize=4, label='Gas mass')
  ax.plot(rad, np.log10(darkmatter_tot), 'm+', markersize=4, label='DM mass')
  ax.set_title('Cumulative Mass versus Radius')
  ax.set_xlabel('Radius (kpc)')
  ax.set_ylabel('Log Mass (m_sun)')


def log_total_enclosed_mass(rad, mstars, mgas, vel):
  """ Part 5 - Calculate the total mass of the galaxy enclosed within each radius
  using the Motions Find Mass equation.
  """

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

  """!!QUESTIONS!! Is totmassEnc a cumulative or incremental radial distribution of the
  total mass of the galaxy? (Hint: think about what the mass in the motions find mass
  equation is.) Why does it NOT make sense to include a plot of totmassEnc in plot window
  1 along with our stellar and gas mass distributions? What kind of distribution
  (cumulative or incremental) of total mass should we plot along with the stellar and
  gas mass distributions?"""
  """ It's cumulative because the motions-find-mass uses the entire mass enclosed within the
  given radius r (assuming it is spherically symmetric) to calculate the tangential velocity
  or given mass of an object.  You would want to compare it against other cumulative masses."""

  # We can find the largest radius we have data for using the "len" function.  The len fuction
  # returns the number of elements in the array.  The last element in the array will have an
  # index of len(array)-1.

  maxR = rad[-1]  # maxR = largest radius in data set = 21.978 kpc

  # Find the total mass, the stellar mass, and the gas mass enclosed in this radius.

  totalmass = (vel[-1] ** 2) * rad[-1] / GRAV
  print('Total inferred mass within %g kpc: %.3e m_sun' % (maxR, totalmass))
  starmass = sum(mstars)
  print('Star mass within %g kpc: %.3e m_sun' % (maxR, starmass))
  gasmass = sum(mgas)
  print('Gas mass within %g kpc: %.3e m_sun' % (maxR, gasmass))

  """!!QUESTIONS!! What fraction of the total mass enclosed within 22 kpc is dark matter?
  Make your python script display this result"""
  dmmass = totalmass - starmass - gasmass
  print('DM mass within %g kpc: %.3e m_sun' % (maxR, dmmass))


def plot_total_incremental_mass(ax, rad, vel):
  """ Part 6 - Calculate the incremental radial distribution of total mass so that we can
  include it with the other incremental mass distributions (stellar and gas) in plot
  window 1. Use the commands below to calculate the incremental mass.
  """
  totalmass = (vel ** 2) * rad / GRAV
  incrementaltotalmass = np.array(totalmass)
  incrementaltotalmass[1:] -= totalmass[:-1]

  """!!QUESTIONS!! Why is the first element of incrementaltotalmass the enclosed mass?
  Explain in words how we compute the incrementaltotalmass in the last step."""
  """ I changed things around.  We're basically shifting the array one index to the right
  in order to compute m_r - m_(r-1) to get the incremental mass.  Still seems a bit fishy
  in terms fo the actual radial bands whose mass this difference calculates but *shrug* -
  I presume it's close enough to not matter.
  """
  ax.plot(rad, np.log10(incrementaltotalmass), 'gx', markersize=4, label='Total mass')


def plot_dm_incremental_mass(ax, rad, mstars, mgas, vel):
  """ Part 7 - Calculate the radial dark matter distribution from the incrementaltotalmass,
  mstars, and mgas arrays that we have made in this Python script. The dark matter
  is the leftover mass when you subtract stellar and gas mass from total mass.
  """
  # Write out how to compute the incremental dark matter distribution in the Python code below.
  # Be sure to check that the dark matter mass is always postive.  If the sum of stellar and gas
  # mass exceeds the total mass inferred from the rotation speed, then we should conclude that
  # there is no evidence for dark matter in that annulus.
  totalmass = (vel ** 2) * rad / GRAV
  incrementaltotalmass = np.array(totalmass)
  incrementaltotalmass[1:] -= totalmass[:-1]
  darkmatter = np.clip(incrementaltotalmass - mstars - mgas, a_min=1e-20, a_max=None)

  ax.plot(rad, np.log10(darkmatter), 'm+', markersize=4, label="DM mass")

  """!!QUESTIONS!! How does the radial distribution of dark matter compare to the radial
  distribution of stars and gas? What questions does this raise to your mind?"""
  """ Gas mass is almost irrelevant.  It is multiple orders of magnitude smaller than the
  mass present in stars within the radial window we are considering.  That said, gas mass
  seems to be more prominent outside the galactic core.  Stars dominate the inernal bulk of
  the galaxy.  Dark Matter (or missing matter) seems to just linger along the outskirts
  of the galaxy and seems to increase in abundance as we look away from the galactic core.

  Does it ever stop getting more abundant?  It's still interacting gravitationally...  why
  such a peculiar distribution where it doesn't coalesce.  Is it because coalescing in the
  middle of a galaxy requires collisions (i.e. EM interactions)?
  """


def main():
  rad, mstars, mgas, vel = load_file('massdecompdata.txt')
  figures = []

  fig1, ax1 = plt.subplots(1, 1)
  figures.append(fig1)
  plot_incremental_mass(ax1, rad, mstars, mgas)

  fig2, ax2 = plt.subplots(1, 1)
  figures.append(fig2)
  plot_rotation_curve(ax2, rad, vel)

  plot_expected_velocity(ax2, rad, mstars, mgas, vel)

  fig3, ax3 = plt.subplots(1, 1)
  figures.append(fig3)
  plot_cumulative_mass(ax3, rad, mstars, mgas, vel)

  log_total_enclosed_mass(rad, mstars, mgas, vel)

  plot_total_incremental_mass(ax1, rad, vel)

  plot_dm_incremental_mass(ax1, rad, mstars, mgas, vel)

  with pdf.PdfPages("ps2plots.pdf") as out:
    for fig in figures:
      for ax in fig.get_axes():
        _, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax * 1.4)
        ax.legend()

      out.savefig(fig)


if __name__ == '__main__':
  sys.exit(main())
