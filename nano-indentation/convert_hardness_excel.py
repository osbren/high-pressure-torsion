import pandas as pd
from pathlib import Path
import numpy as np

name = 'W20Cu20T'
path = Path(name + '.xlsx')
# Assume indexing goes from centre to far right, then left of centre to far left
centre = 2
left_of_centre = 82
far_left = 160
far_right = left_of_centre - 1
step = 50

radial_profile = pd.DataFrame({'Radius': [], 'Hardness': [], 'Modulus': []})
print('Starting at centre point')
for i in range(centre, far_right + 1):  # Right points

    sheet = pd.read_excel(path, sheet_name=f'Test {i:03d}', usecols='B:G', dtype=np.float32, skiprows=[1])

    mean_hardness = sheet[(sheet['Displacement Into Surface'].between(500, 1500)) &
                          sheet['Hardness'].between(0, 20)]['Hardness'].mean()
    mean_modulus = sheet[(sheet['Displacement Into Surface'].between(500, 1500)) &
                         sheet['Modulus'].between(50, 500)]['Modulus'].mean()

    radial_profile.loc[i] = [(i - centre) * step, mean_hardness, mean_modulus]
print('Reached rightmost point')
for i in range(left_of_centre, far_left + 1):  # Left points
    sheet = pd.read_excel(path, sheet_name=f'Test {i:03d}', usecols='B:G', dtype=np.float32, skiprows=[1])

    mean_hardness = sheet[(sheet['Displacement Into Surface'].between(500, 1500)) &
                          sheet['Hardness'].between(0, 20)]['Hardness'].mean()
    mean_modulus = sheet[(sheet['Displacement Into Surface'].between(500, 1500)) &
                         sheet['Modulus'].between(50, 500)]['Modulus'].mean()
    radial_profile.loc[i] = [(left_of_centre - i - 1) * step, mean_hardness, mean_modulus]
print('Reached leftmost point')

# Store pandas DataFrame (getting Excel data is slow)
radial_profile.to_pickle(name + '.pkl')
print(radial_profile)
