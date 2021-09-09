"""
Downloads test trajectories and resaves a shortened version in xtc and the mdtraj hdf5 format
"""

import mdshare
import pathlib  
import mdtraj

THIS_DIR=pathlib.Path(__file__).parent

files = mdshare.fetch('alanine-dipeptide-*-250ns-nowater.xtc', working_directory=THIS_DIR, force=True)
pdb = mdshare.fetch('alanine-dipeptide-nowater.pdb', working_directory=THIS_DIR, force=True)

for file in files:
    file_name = str(pathlib.Path(file).stem)
    traj = mdtraj.load(file, top=pdb)
    pathlib.Path(file).unlink()
    traj = traj[:int(traj.n_frames*0.2)]
    traj.save_hdf5(str(THIS_DIR.joinpath(file_name+'.h5')))
    traj.save_xtc(str(THIS_DIR.joinpath(file_name+'.xtc')))


