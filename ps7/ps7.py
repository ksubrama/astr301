import astropy.constants as ac
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import numpy as np
import scipy.constants as sc
import scipy.optimize as opt
import sys


def make_plot(ax, rmax):
  Mbh = (10 ** 7) * ac.M_sun.value   # Mass of central BH.
  R = 2 * 1000 * sc.parsec           # Radius of galaxy.
  v_max = 200000                     # Max rotational velocity from star matter (m/s).

  rs = np.linspace(0.001, rmax, 10000)  # In kpc
  rs_m = rs * 1000 * sc.parsec

  ax.set_xlabel("Radius from center (Kpc)")
  ax.set_ylabel("Rotational velocity (Km/s)")

  v_bh = np.sqrt(Mbh * sc.G / rs_m)
  ax.plot(rs, v_bh / 1000, 'b-', label="Black Hole only")

  v_curve = np.clip(v_max * rs_m / R, 0, v_max)
  ax.plot(rs, v_curve / 1000, 'g--', label="Galactic matter only")

  M_enc = (v_curve ** 2 * rs_m / sc.G) + Mbh
  v_enc = np.sqrt(M_enc * sc.G / rs_m)
  ax.plot(rs, v_enc / 1000, 'r:', label="Total matter")

  # Radius of influence
  r_inf = opt.fsolve(lambda r: np.sqrt(Mbh * sc.G / r) - (v_max * r / R), sc.parsec)[0]
  v_inf = v_max * r_inf / R
  r_inf_kpc = r_inf / (1000 * sc.parsec)
  v_inf_kms = v_inf / 1000
  ax.annotate(
      "radius of influence = %s kpc" % round(r_inf_kpc, 3),
      xy=(r_inf_kpc, v_inf_kms + 10), xycoords='data',
      xytext=(0.54, 0.3), textcoords='axes fraction',
      arrowprops=dict(arrowstyle='->', connectionstyle="angle,angleA=0,angleB=90,rad=5"),
      bbox=dict(boxstyle="square,pad=0.3", facecolor="white"))

  return r_inf


def main():
  figures = []

  fig1, ax1 = plt.subplots(1, 1)
  figures.append(fig1)
  make_plot(ax1, 0.5)

  fig2, ax2 = plt.subplots(1, 1)
  figures.append(fig2)
  make_plot(ax2, 3)

  fig3, ax3 = plt.subplots(1, 1)
  make_plot(ax3, 20)
  figures.append(fig3)

  with pdf.PdfPages("ps7plots.pdf") as out:
    for fig in figures:
      for ax in fig.get_axes():
        _, xmax = ax.get_xlim()
        ax.set_xlim(0, xmax)
        _, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax * 1.3)
        ax.legend()
      out.savefig(fig)


if __name__ == '__main__':
  sys.exit(main())
