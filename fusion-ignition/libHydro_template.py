################################################################################
################################################################################
#  File:           libHydro.py
#  Description:    Set of functions for Lagrangian hydro evolution in 1D
#  Notes:
#  Author:         Britton Olson
#  Date:           11/3/2021
################################################################################
################################################################################
import numpy as np
import sys

sys.path.append("../../libCommon/python")
from libCommon import makeMeshVariables


def setupHydroState(Nzones):
    """
    inputs:
      Nzones - total number of zones in domain
    returns:
      hydro_state - a dictionary with the hydro state nodals/zonals
    """
    # List of variables; Nodal and Zonal
    zonals = ['p','den','e','T','vol','dx','xz','Zmass','q','cs2','mat','work','cv']
    nodals = ['x','xdot','Nmass','Farea','dxn','acc']
    hydro_state = makeMeshVariables(Nzones,zonals,nodals)
    return hydro_state 
    

def updateVolume(hd):
    """
    inputs:
      hd   - hydro state created in setup_hydro
    returns:
      None - hydro state volume-related variables updated
    """
    geom = hd.get('geom', 0)                                 # 0-planar, 1-cylindrical, 2-spherical
    assert geom == 0, "Error: updateVolume only support geom=0 (planar)"

    hd['dx']  = np.diff( hd['x'] )                           # X (radial) spacing
    hd['xz']  = hd['dx']/2.0 + hd['x'][:-1]                  # Zone centered coordinate
    hd['vol'] = hd['dx']                                     # Volume of the zone
    hd['Farea'] = np.ones_like(hd['x'])                      # Face area
    hd['dxn'] = hd['x']*0.0 
    hd['dxn'][1:-1] = np.diff( hd['xz'] )                    # Node centered spacing
    hd['dxn'][0] = hd['dx'][0]
    hd['dxn'][-1] = hd['dx'][-1]
    
    return None


def updateMass(hd):
    """
    inputs:
      hd   - hydro state
    returns:
      None - hydro state 'Zmass' and 'Nmass' updated
    """
    hd['Zmass'] = hd['vol'] * hd['den']
    hd['Nmass'][1:-1] = ( hd['Zmass'][:-1] +
                          hd['Zmass'][1:]  ) / 2.0
    hd['Nmass'][0]  = hd['Zmass'][0]  / 2.0
    hd['Nmass'][-1] = hd['Zmass'][-1] / 2.0

    return None


def updateQ(hd):
    """
    inputs:
      hd   - hydro state
    return:
      None   - hydro state 'q' (artificial viscosity) updated
    """
    
    print("Warning: UpdateQ is not implemented")

    return None
    
def updateAcceleration(hd):
    """
    inputs:
      hd  - hydro state
    returns:
      None  - hydro state 'acc' (accleration) updated
    """
    
    print("Warning: UpdateAcceleration is not implemented")

    return None

def updateWork(hd):
    """
    inputs:
      hd   - hydro state
    returns:
      None   - hydro state 'work' (zonal based pressure work term... F*dV) updated
    """

    print("Warning: UpdateWork is not implemented")
    
    return None

def hydroStep(dt,hydro_state):

    
    print("Warning: hydroStep is not implemented")

    return None



