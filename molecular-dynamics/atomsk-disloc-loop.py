import os
import numpy as np
import glob
import csv
from pathlib import Path

class Element:
    def __init__(self, symbol, structure):
        self.symbol = symbol
        if structure == 'FCC':
            self.slip_planes = np.asarray(((1, 1, 1), (1, -1, 1), (1, 1, -1), (1, -1, -1)))
            self.slip_dirs = np.asarray(((1, 1, 0), (1, 0, 1), (1, -1, 0), (1, 0, -1), (0, 1, 1), (0, 1, -1)))
        elif structure == 'BCC':
            self.slip_planes = np.asarray(((1, 1, 0), (1, 0, 1), (1, -1, 0), (1, 0, -1), (0, 1, 1), (0, 1, -1)))
            self.slip_dirs = np.asarray(((1, 1, 1), (1, -1, 1), (1, 1, -1), (1, -1, -1)))
        else:
            raise ValueError

    def setLatticeConstant(self, latticeConst):
        self.latticeConst = latticeConst
        self.burgersVector = np.sqrt(3)/2 * latticeConst

    def setPoisson(self, poisson):
        self.poisson = poisson


class DislocLoops:
    def __init__(self, element):
        self.element = element

    def set_loop_sizes(self, rmin, rmax):
        self.rmin = rmin
        self.rmax = rmax

    def set_max_total_radii(self, rtotal):
        self.rtotal = rtotal

    def set_radius_step(self, rstep):
        self.rstep = rstep

    def set_path(self, path):
        self.path = path
        self.empty_filename = self.path + self.element.symbol + '_0'

    def set_cell_size(self, size):
        self.size = size
        self.volume = (size * self.element.latticeConst * 1E-10) ** 3

    def s(self, x):
        y = " "
        for i in x:
            y += str(i)
        return y

    def rand_shift(self):  # Random shift
        y = " "
        for i in np.round(np.random.random(3) / 2, 3):
            y += str(i) + "\\*BOX "
        return y

    def create_empty(self, filename):
        os.system("atomsk --create bcc " + str(self.element.latticeConst) + " " + self.element.symbol + " -duplicate " + (
                str(self.size) + " ") * 3 + filename + ".cfg")  # Create supercell

    def add_loop_to_existing(self, filename, R):
        slip_plane = self.element.slip_planes[np.random.randint(0, len(self.element.slip_planes))]
        slip_dir = np.asarray([7,8,9])
        while np.dot(slip_plane, slip_dir) != 0:
            slip_dir = self.element.slip_dirs[np.random.randint(0, len(self.element.slip_dirs))]  # Generate random valid slip system

        perp_dir = np.cross(slip_dir, slip_plane)
        # Atomsk does not like first direction to be negative
        perp_dir = perp_dir * np.sign(perp_dir[0])

        # Add dislocation loop orientated in random slip system
        os.system("atomsk " + filename + ".cfg -orient 100 010 001" + self.s(slip_dir) + self.s(perp_dir) + self.s(slip_plane) +
                  " -disloc loop 0.501\\*box 0.501\\*box 0.501\\*box Z " + str(R) +
                  " " + str(self.element.burgersVector) + " 0 0 " + str(self.element.poisson) + " \\" +
                  "-alignx  -shift" + self.rand_shift() + " -wrap cfg -ow")

    def generate_lammps_output(self, filename):
        # Generate LAMMPS data file
        os.system("atomsk " + filename + ".cfg lammps")
        # Add support for triclinic co-ordinates - otherwise LAMMPS gives error
        # This works for Mac, needs to be changed for windows?
        os.system("gsed -i '/zhi/a \\      0.0000000    0.00000000    0.00000000  xy xz yz' " + filename + ".lmp")

    def copy_rename(self, prev_filename, filename):
        os.system("atomsk " + prev_filename + ".cfg " + filename + ".cfg ")

    def remove_files(self):
        # Remove previous files
        for file in glob.glob(self.path + '*.cfg') + glob.glob(self.path + '*.lmp'):
            os.remove(file)

    def generate_loops(self, lammps=True):
        # Undeformed lattice
        count = 0
        # Create folder in given path if not already there
        Path(self.path).mkdir(parents=True, exist_ok=True)
        self.create_empty(self.empty_filename)
        disloc_density = 0
        if lammps:
            self.generate_lammps_output(self.empty_filename)
        data = [count, disloc_density]
        with open('disloc.csv', 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(data)

        # Adding loops
        for R in range(self.rmin, self.rtotal, self.rstep):
            count += 1
            disloc_density = 2E-10 * np.pi * R / self.volume
            filename = self.path + self.element.symbol + '_' + str(count)
            r = R
            self.copy_rename(self.empty_filename, filename)

            while r > self.rmax:
                if r < self.rmax + self.rmin:
                    r -= self.rmin
                    self.add_loop_to_existing(filename, self.rmin)
                else:
                    r -= self.rmax
                    self.add_loop_to_existing(filename, self.rmax)

            self.add_loop_to_existing(filename, r)

            self.generate_lammps_output(filename)
            data = [count, disloc_density]

            with open('disloc.csv', 'a') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(data)
