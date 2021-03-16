from atomsk_disloc_loop import *

molybdenum = Element('Mo', 'BCC')
molybdenum.setLatticeConstant(3.14)
molybdenum.setPoisson(0.28)

mo_loops = DislocLoops(molybdenum)
mo_loops.set_cell_size(40)
mo_loops.set_path('/private/tmp/LAMMPS/')
mo_loops.set_loop_sizes(rmin=20, rmax=40)
mo_loops.set_max_total_radii(300)
mo_loops.set_radius_step(20)

mo_loops.remove_files()
mo_loops.generate_loops(lammps=True)
