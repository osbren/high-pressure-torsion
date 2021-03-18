#!/usr/bin/env python
'''from ovito.io import import_file
from ovito.modifiers import DislocationAnalysisModifier
from ovito.data import DislocationNetwork'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import csv
from sklearn import linear_model

SUFFIX ="_Cu_Mishin2001EAM_40"
LINEAR_DISLOC_LIMIT = 3E16

plot_linear_fit = False
plot_ransac = True

with open('elastic' + SUFFIX + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file)
    l = -1  # Ignore first row which is the column headings
    for row in csv_reader:
        l += 1
disloc_density = np.zeros(l)
C_cubic = np.zeros(shape=(3,l))

'''for i in range(l):

    filename = "Cu_" + str(i) + ".cfg"

    pipeline = import_file(filename)
    # Extract dislocation lines from a crystal with diamond structure:
    modifier = DislocationAnalysisModifier()
    modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
    pipeline.modifiers.append(modifier)
    data = pipeline.compute()

    total_line_length = data.attributes['DislocationAnalysis.total_line_length'] * 1E-10
    cell_volume = data.attributes['DislocationAnalysis.cell_volume'] * 1E-30
    disloc_density[i] = total_line_length / cell_volume
    print("Dislocation density: {:.3e}".format(disloc_density[i]))
    data = [i, disloc_density[i]]
    with open('disloc_ovito.csv', 'a') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(data)'''

with open('elastic' + SUFFIX + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file)
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            i = int(row[0])
            C_cubic[:,i] = np.asarray([float(x) for x in row[1:]])
        line_count += 1
    print(f'Processed {line_count} lines.')

with open('disloc' + SUFFIX + '.csv') as csv_file:
    csv_reader = csv.reader(csv_file)
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            i = int(row[0])
            disloc_density[i] = float(row[1])
        line_count += 1
    print(f'Processed {line_count} lines.')

ransac = linear_model.RANSACRegressor()
ransacX = np.vstack(disloc_density[disloc_density < LINEAR_DISLOC_LIMIT])
ransacY = np.vstack(C_cubic[:, disloc_density < LINEAR_DISLOC_LIMIT])
slope, intercept, r, p, se = stats.linregress(disloc_density, C_cubic[0])
ax = plt.subplot(1, 3, 1)
#ax.set_title('Elastic modulus change with dislocations')
ax.set_xlabel('Dislocation Density / m$^{-2}$')
ax.set_ylabel("Young's Modulus / GPa")
plt.scatter(disloc_density, C_cubic[0])
if plot_linear_fit:
    plt.plot(disloc_density, disloc_density*slope + intercept, 'r')
if plot_ransac:
    ransac.fit(ransacX, ransacY[0])
    line_y_ransac = ransac.predict(np.vstack(disloc_density))
    plt.plot(disloc_density, line_y_ransac, 'g', label='Slope = {:.2e} GPam$^2$'.format(ransac.estimator_.coef_[0]))
    plt.legend()
print(p)


slope, intercept, r, p, se = stats.linregress(disloc_density, C_cubic[1])
ax = plt.subplot(1, 3, 2)
ax.set_xlabel('Total dislocations / m$^{-2}$')
ax.set_ylabel("C12 / GPa")
plt.scatter(disloc_density, C_cubic[1])
if plot_linear_fit:
    plt.plot(disloc_density, disloc_density*slope + intercept, 'r')
if plot_ransac:
    ransac.fit(ransacX, ransacY[1])
    line_y_ransac = ransac.predict(np.vstack(disloc_density))
    plt.plot(disloc_density, line_y_ransac, 'g', label='Slope = {:.2e} GPam$^2$'.format(ransac.estimator_.coef_[0]))
    plt.legend()
print(p)

slope, intercept, r, p, se = stats.linregress(disloc_density, C_cubic[2])
ax = plt.subplot(1, 3, 3)
ax.set_xlabel('Total dislocations / m$^{-2}$')
ax.set_ylabel("Shear Modulus / GPa")
plt.scatter(disloc_density, C_cubic[2])
if plot_linear_fit:
    plt.plot(disloc_density, disloc_density*slope + intercept, 'r')
if plot_ransac:
    ransac.fit(ransacX, ransacY[2])
    line_y_ransac = ransac.predict(np.vstack(disloc_density))
    plt.plot(disloc_density, line_y_ransac, 'g', label='Slope = {:.2e} GPam$^2$'.format(ransac.estimator_.coef_[0]))
    plt.legend()
print(p)
plt.show()