import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
  # Oh boy, this where I self-troll.
  # Does "number letters" include hyphens or not?
  myarray = np.arange(1, len('Cating-Subramanian') + 1)
  print(myarray)
  rootarray = myarray ** 0.5
  print(rootarray)
  ratio = myarray / rootarray
  print(ratio)
  print(rootarray * ratio)
  # sqrt(a) * (a / sqrt(a)) == a

  data = np.loadtxt('ps2tutorial.data')
  temperature = data[:, 0]
  humidity = data[:, 1]
  print(humidity[humidity < 20])
  print(temperature[humidity < 20])

  plt.plot(humidity, temperature, 'b.', markersize=12)
  plt.title('Jane Doe Python Tutorial')
  plt.xlabel('humidity (%)')
  plt.ylabel('temperature (F)')
  plt.xlim(10, 60)
  plt.ylim(75, 100)
  plt.plot(
      humidity[(temperature > 80) & (temperature < 100)],
      temperature[(temperature > 80) & (temperature < 100)],
      'g*',
      markersize=15)
  plt.plot(
      humidity[(humidity >= 10) & (humidity <= 40)],
      temperature[(humidity >= 10) & (humidity <= 40)],
      'r+',
      markersize=15)
  plt.savefig('ps2tutorial.png')


if __name__ == '__main__':
  sys.exit(main())
