################################################################################
################################################################################
#  File:           libCommmon.py
#  Description:    Set of common functions for multi-physics library
#  Notes:
#  Author:         Britton Olson
#  Date:           04/28/2021
################################################################################
################################################################################
import numpy as np


def makeMeshVariables(Nzones,zonals,nodals):
    """
    inputs:
      Nzones - total number of zones in domain
      zonals - list of strings, the names of the zone centered variables
      nodals - list of strings, the names of the node centered variables
    returns:
      state - a dictionary with numpy arrays set to zero
    """
    # Make the data structures
    state = {}
    for nv in nodals:
        state[nv] = np.zeros( Nzones + 1 )
        
    for nz in zonals:
        state[nz] = np.zeros( Nzones )

    return state


def gradZ(zonal,dx):
    """
    inputs:
      zonal - Zone centered field of which to take the derivative
      dx    - Zone centered grid spacing
    returns:
      gradN - Nodal value of gradient of the zone centered field, 2nd order
    """

    gradN = np.zeros(np.size(zonal)+1)
    
    gradN[1:-1] = ( zonal[1:] - zonal[0:-1] ) / (.5*(dx[1:] + dx[0:-1] ) )
    gradN[0]    = gradN[1]
    gradN[-1]   = gradN[-2]

    return gradN

def zoneToNode(zonal,dx):
    """ 
    zoneToNode interpolates the zonal field "zonal" to a node centered field using the grid spacing "dx".
      :param zonal: zonal centered field 
      :param dx: array of grid size
      :return nodal: Interpolated nodal, centered field.
    """
    
    nodal = np.zeros(np.size(zonal)+1)
    
    for i in range(1,np.size(zonal)):
        nodal[i]=(zonal[i]-zonal[i-1])/(dx[i]/2.+dx[i-1]/2.)*(dx[i-1]/2.)+zonal[i-1]
    nodal[0] = (zonal[1]-zonal[0])/(dx[1]/2.+dx[0]/2.)*(-dx[0]/2.)+zonal[0]
    nodal[-1] = (zonal[-1]-zonal[-2])/(dx[-1]/2.+dx[-2]/2.)*(dx[-2]/2.+dx[-1])+zonal[-1]
    return nodal 
