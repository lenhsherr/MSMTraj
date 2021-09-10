import mdtraj as md
import numpy as np
from pathlib import Path
from typing import Union, Dict
import pandas as pd


SELF_CONTAINED_FMTS = ['h5', 'hdf5', 'npy', 'npz', 'parquet']


class Trajectory(): 
    def __init__(self, configuration: Dict[str, Union[str, Path, None]]) -> None:
        self._coordinates_path: Path = None
        self._topology_path: Path = None
        self._coordinates: np.ndarray = None
        self._timestep: float = None
        self._file_type: str = None
        self._parse_configuration(configuration)

    @property
    def file_type(self) -> str:
        return self._file_type

    @file_type.setter
    def file_type(self, value) -> None:
        self._file_type = value

    @property
    def coordinates(self) -> np.ndarray:
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value) -> None: 
        self._coordinates = value

    @property
    def coordinates_path(self) -> Path:
        return self._coordinates_path

    @coordinates_path.setter
    def coordinates_path(self, value) -> None:
        try: 
            value = Path(value)
        except TypeError: 
            raise TypeError('Coordinates path must be str or Path type')
        self._coordinates_path = value

    @property
    def topology_path(self) -> Path:
        return self._topology_path
    
    @topology_path.setter
    def topology_path(self, value) -> None: 
        try: 
            value = Path(value)
        except TypeError: 
            raise TypeError('Topology path must be str or Path type')
        self._topology_path = value
        
    def _parse_configuration(self, configuration: Dict[str, Union[str, Path, None]]) -> None: 
        self.coordinates_path = configuration['coordinates_path']
        self.file_type = self.coordinates_path.suffix[1:]
        
        if not (self.file_type in SELF_CONTAINED_FMTS): 
            try: 
                self.topology_path = configuration['topology_path']
            except KeyError:
                raise KeyError(f'topology_path needed for coordinate files of type {self.file_type}')

    def read(self) -> None:
        if self.file_type in ['h5', 'hdf5', 'xtc', 'pdb']:
            kwargs = {}
            if not (self.file_type in SELF_CONTAINED_FMTS):
                kwargs['top'] = str(self.topology_path)

            traj = md.load(str(self.coordinates_path), **kwargs)
            n_frames = traj.n_frames
            self.coordinates = traj.xyz.reshape(n_frames, -1)

        elif self.file_type in ['npy']: 
            traj = np.load(self.coordinates_path)
            self.coordinates = traj

        else: 
            raise NotImplementedError(f"file type {self.file_type} is not supported yet.")
