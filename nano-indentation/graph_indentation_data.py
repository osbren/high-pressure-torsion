import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm

names = ['W20Cu5T', 'W20Cu20T']  # List of converted nano-indentation profiles
fig, axes = plt.subplots(1, 2)
colourmap = cm.get_cmap('brg')  # See https://matplotlib.org/stable/tutorials/colors/colormaps.html

for name in names:
    df = pd.read_pickle(name + '.pkl')
    frac = (names.index(name) + 1) / len(names)
    df.plot.scatter(x='Radius', xlabel=r'Radius ($\mu$m)', y='Hardness', ylabel='Hardness (GPa)', ax=axes[0],
                    legend=True, label=name, color=colourmap(frac))
    df.plot.scatter(x='Radius', xlabel=r'Radius ($\mu$m)', y='Modulus', ylabel='Modulus (GPa)', ax=axes[1],
                    legend=True, label=name, color=colourmap(frac))

plt.show()
