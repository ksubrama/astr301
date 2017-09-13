import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import os
import sys


GRAV = 4.28e-6  # define G in units of kpc*km^2/(M_sun*s^2)


class DwarfSpheroidals:

  def __init__(self, name, radius, luminosity, file_path):
    self.name = name
    self.radius = radius
    self.luminosity = luminosity
    # read in the data columns: vel contains the radial velocites of all the stars,
    # and prob contains the probability that each star is a member of the dSph.
    (self.target_id,
     self.ra_hr, self.ra_min, self.ra_sec,
     self.dec_deg, self.dec_min, self.dec_sec,
     self.vmag, self.vicol, self.all_vel, self.vel_errs,
     self.prob) = np.loadtxt(file_path, unpack=True)
    self._vel = None
    self._mean_vel = None
    self._vel_disp = None
    self._virial_mass = None
    self._mass_light_ratio = None

  def _is_member(self):
    return self.prob > 0.9

  def vel(self):
    if self._vel is None:
      self._vel = self.all_vel[self._is_member()]
    return self._vel

  def plot_hist_90(self, ax, bins=50):
    ax.hist(self.vel(), 50)
    ax.set_title(self.name)
    # Don't set axis labels here - they'll be set as a super label.
    return self.vel()

  def mean_vel(self):
    if self._mean_vel is None:
      self._mean_vel = np.mean(self.vel())
    return self._mean_vel

  def vel_disp(self):
    if self._vel_disp is None:
      self._vel_disp = np.sum((self.vel() - self.mean_vel()) ** 2) / len(self.vel())
    return self._vel_disp

  def virial_mass(self):
    if self._virial_mass is None:
      self._virial_mass = 5 * self.vel_disp() * self.radius / GRAV
    return self._virial_mass

  def mass_light_ratio(self):
    if self._mass_light_ratio is None:
      self._mass_light_ratio = self.virial_mass() / self.luminosity
    return self._mass_light_ratio


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


def table_row(name, members, mean_vel, vel_disp, virial_mass, mass_light_ratio):
  return "%15s %15s %15s %15s %15s %15s" % (name, members, mean_vel, vel_disp, virial_mass, mass_light_ratio)


def table_data_row(name, members, mean_vel, vel_disp, virial_mass, mass_light_ratio):
  return table_row(name, members, "%.3e" % mean_vel, "%.3e" % vel_disp, "%.3e" % virial_mass, "%.3e" % mass_light_ratio)


def load_data():
  # name, radius (kpc), luminosity (LSun)
  data = [
      ("Carina", 0.85, 4.3e5),
      ("Fornax", 2.70, 15.5e6),
      ("Sculptor", 1.63, 2.15e6),
      ("Sextans", 4.00, 5e5)]
  return [DwarfSpheroidals(name, rad, lum, os.path.join("dSph", "dSph_%s.dat" % name))
          for (name, rad, lum) in data]


def main():
  figures = []
  dSphs = load_data()

  fig1, axes1 = plt.subplots((len(dSphs) + 1) // 2, 2, sharex=True, sharey=True)
  figures.append(fig1)
  super_label(fig1, "Radial Velocity [km/s]", "Number of stars P > 0.9")

  print(table_row("name", "# of members", "mean velocity", "vel dispersion", "virial mass", "mass/light"))
  print(table_row("", "", "km/s", "(km/s)^2", "Msun", "Msun/Lsun"))
  for (i, dSph) in enumerate(dSphs):
    nstar = len(dSph.plot_hist_90(axes1[i // 2, i % 2]))
    mean_vel = dSph.mean_vel()
    vel_disp = dSph.vel_disp()
    virial_mass = dSph.virial_mass()
    mass_light_ratio = dSph.mass_light_ratio()
    print(table_data_row(dSph.name, nstar, mean_vel, vel_disp, virial_mass, mass_light_ratio))

  with pdf.PdfPages("ps3plots.pdf") as out:
    for fig in figures:
      # for ax in fig.get_axes():
        # _, ymax = ax.get_ylim()
        # ax.set_ylim(0, ymax * 1.4)
        # ax.legend()

      out.savefig(fig)


if __name__ == '__main__':
  sys.exit(main())
