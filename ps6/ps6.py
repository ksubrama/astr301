import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import numpy as np
import scipy.constants as sc
import scipy.optimize as opt
import sys


def planck_lam(lam, T):
  return ((2 * sc.h * (sc.c ** 2) / (lam ** 5)) /
          (np.exp(sc.h * sc.c / (lam * sc.k * T)) - 1))


def planck_nu(nu, T):
  return ((2 * sc.h * (nu ** 3) / (sc.c ** 2)) /
          (np.exp(sc.h * nu / (sc.k * T)) - 1))


def plot_planck_wavelength(ax, temps):
  lams = np.linspace(100, 2000, 1000)  # In nm
  for idx, (T, style) in enumerate(temps):
    ax.plot(lams, planck_lam(lams / 1e9, T), style, label='T = %s K' % T)
    max_res = opt.minimize_scalar(
        lambda l: - planck_lam(l / 1e9, T),
        method='brent',
        bracket=(lams[0], lams[-1]))
    lam_max = max_res.x
    B_max = planck_lam(lam_max / 1e9, T)
    ax.annotate(
        "lam = %s nm" % int(round(lam_max)),
        xy=(lam_max, B_max), xycoords='data',
        xytext=(0.75, 0.2 + idx * 0.2), textcoords='axes fraction',
        arrowprops=dict(arrowstyle='->', connectionstyle="angle,angleA=0,angleB=90,rad=5", edgecolor=style[0]),
        bbox=dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor=style[0]))
    print("Lambda max for %s K = %s nm (%s THz)" % (T, lam_max, (sc.c / 1e3) / lam_max))

  ax.set_title('Black body spectral radiance by wavelength')
  ax.set_xlabel('Wavelength (nm)')
  ax.set_ylabel('Spectral radiance (W m-2 sr-1 m-1)')


def plot_planck_frequency(ax, temps):
  nus = np.linspace(1, 1500, 1000)  # In nm
  for (T, style) in temps:
    ax.plot(nus, planck_nu(nus * 1e12, T), style, label='T = %s K' % T)
  ax.set_title('Black body spectral radiance by frequency')
  ax.set_xlabel('Frequency (THz)')
  ax.set_ylabel('Spectral radiance (W m-2 sr-1 Hz-1)')


def print_root():
  # Our initial guess using 580nm at 5000K
  x0 = sc.h * sc.c / (580e-9 * sc.k * 5000)
  x = opt.fsolve(lambda x: x * np.exp(x) / (np.exp(x) - 1) - 5, x0)
  print("x to get lambda max: %s" % x[0])
  wein_factor = sc.h * sc.c * 1000 / (x * sc.k)
  print(" => lambda = %s (mm/K) / T" % wein_factor)


def main():
  figures = []
  temps = [(4000, 'r-'), (5000, 'g--'), (6000, 'b-.')]

  fig1, axes1 = plt.subplots(1, 1)
  figures.append(fig1)
  plot_planck_wavelength(axes1, temps)

  fig2, axes2 = plt.subplots(1, 1)
  figures.append(fig2)
  plot_planck_frequency(axes2, temps)

  with pdf.PdfPages("ps6plots.pdf") as out:
    for fig in figures:
      for ax in fig.get_axes():
        _, xmax = ax.get_xlim()
        ax.set_xlim(0, xmax)
        _, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax * 1.3)
        ax.legend()
      out.savefig(fig)

  print_root()


if __name__ == '__main__':
  sys.exit(main())
