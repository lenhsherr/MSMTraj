import mdtraj as md
import numpy as np
from pathlib import Path
from typing import Union, Dict

SELF_CONTAINED_FMTS = ['h5', 'hdf5', 'npy', 'npz']


class Trajectory(): 
    def __init__(self, configuration: Dict[str, Union[str, Path, None]]) -> None:
        self._coordinates_path: Path = None
        self._topology_path: Path = None
        self._coordinates: np.ndarray = None
        self._timestep: float = None

        self._parse_configuration(configuration)


    @property
    def coordinates_path(self) -> Path:
        return self._coordinates_path

    @property.setter
    def coordinates_path(self, value) -> None:
        try: 
            value = Path(value)
        except TypeError: 
            raise TypeError('Coordinates path must be str or Path type')
        self._coordinates_path = value

    @property
    def topology_path(self) -> Path:
        return self._topology_path
    
    @property.setter
    def topology_path(self, value) -> None: 
        try: 
            value = Path(value)
        except TypeError: 
            raise TypeError('Topology path must be str or Path type')
        self._topology_path = value
        
    def _parse_configuation(self, configuration: Dict[str, Union[str, Path, None]]) -> None: 
        self.coordinates_path = configuration['configuation_path']
        file_type = self.coordinates_path.suffix
        if not file_type in SELF_CONTAINED_FMTS: 
            try: 
                self.topology_path = configuration['topology_path']
            except KeyError:
                raise KeyError(f'topology_path needed for coordinate files of type {file_type}')
        


    # def _parse_configuration(self, config: Mapping) -> None: 
    #     self.coordinates_path = Path(config['coordinates_path'])
    #     if 'topology_path' in config.keys(): 
    #         if config['topology_path'] is not None: 
    #             self.topology_path = Path(config['topology_path'])
    #         else: 
    #             raise ValueError('topology_path can not be None')

    # def _requires_topology(self):
    #     return not self.coordinates_path.suffix in SELF_CONTAINED_FMTS
    

    # def read(self) -> None:

