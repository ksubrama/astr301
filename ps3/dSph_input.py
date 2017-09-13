import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import os
import re
import sys


GRAV = 4.28e-6  # define G in units of kpc*km^2/(M_sun*s^2)


class DwarfSpheroidals:

  def __init__(self, name, file_path):
    self.name = name
    # read in the data columns: vel contains the radial velocites of all the stars,
    # and prob contains the probability that each star is a member of the dSph.
    (self.target_id,
     self.ra_hr, self.ra_min, self.ra_sec,
     self.dec_deg, self.dec_min, self.dec_sec,
     self.vmag, self.vicol, self.vel, self.vel_errs,
     self.prob) = np.loadtxt(file_path, unpack=True)

  def plot_hist_90(self, ax, bins=50):
    vel = self.vel[self.vel > 0.9]
    ax.hist(vel, bins=bins)
    ax.set_title(self.name)
    # Don't set axis labels here - they'll be set as a super label.
    return vel


def super_label(fig, xlabel, ylabel):
  """Make some commmon axis labels for a grid of sub-plots.

  This uses an artificial super-plot that covers all subplots but turns off its
  axis rendering.
  """
  ax = fig.add_subplot(111, frameon=False)
  ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  return ax


def main():
  figures = []
  dSphs = [DwarfSpheroidals(
      re.match("dSph_(.*)\.dat", os.path.basename(path))[1],
      path) for path in glob.glob(os.path.join("dSph", "*.dat"))]

  fig1, axes1 = plt.subplots((len(dSphs) + 1) // 2, 2, sharex=True, sharey=True)
  figures.append(fig1)
  super_label(fig1, "Radial Velocity [km/s]", "Number of stars P > 0.9")

  # make a histogram of the velocity values with 50 bins
  print("Stars with P > 0.9")
  for (i, dSph) in enumerate(dSphs):
    nstar = len(dSph.plot_hist_90(axes1[i // 2, i % 2]))
    print("%12s : %s" % (dSph.name, nstar))

  with pdf.PdfPages("ps3plots.pdf") as out:
    for fig in figures:
      # for ax in fig.get_axes():
        # _, ymax = ax.get_ylim()
        # ax.set_ylim(0, ymax * 1.4)
        # ax.legend()

      out.savefig(fig)


if __name__ == '__main__':
  sys.exit(main())
