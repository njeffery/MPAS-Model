#!/usr/bin/env python
# Generate initial conditions for hydro-margin land ice test case

import sys
from netCDF4 import Dataset as NetCDFFile
import numpy as np

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup dome", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
gridfile = NetCDFFile(options.filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
xCell = gridfile.variables['xCell']
yCell = gridfile.variables['yCell']
xEdge = gridfile.variables['xEdge']
yEdge = gridfile.variables['yEdge']
xVertex = gridfile.variables['xVertex']
yVertex = gridfile.variables['yVertex']
thickness = gridfile.variables['thickness']
bedTopography = gridfile.variables['bedTopography']
layerThicknessFractions = gridfile.variables['layerThicknessFractions']
SMB = gridfile.variables['sfcMassBal']


# Find center of domain
x0 = xCell[:].min() + 0.5 * (xCell[:].max() - xCell[:].min() )
y0 = yCell[:].min() + 0.5 * (yCell[:].max() - yCell[:].min() )
# Calculate distance of each cell center from dome center
r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5

# Center the dome in the center of the cell that is closest to the center of the domain.
#   NOTE: for some meshes, maybe we don't want to do this - could add command-line argument controlling this later.
putOriginOnACell = True
if putOriginOnACell:
   centerCellIndex = np.abs(r[:]).argmin()
   #print x0, y0, centerCellIndex, xCell[centerCellIndex], yCell[centerCellIndex]
   xShift = -1.0 * xCell[centerCellIndex]
   yShift = -1.0 * yCell[centerCellIndex]
   xCell[:] = xCell[:] + xShift
   yCell[:] = yCell[:] + yShift
   xEdge[:] = xEdge[:] + xShift
   yEdge[:] = yEdge[:] + yShift
   xVertex[:] = xVertex[:] + xShift
   yVertex[:] = yVertex[:] + yShift
   # Now update origin location and distance array
   x0 = 0.0
   y0 = 0.0
   r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5


r0=25000.0
thickness[0, r<r0] = 500.0 * (1.0 - (r[r<r0] / r0)**2)
thickness[0, r>22500.0] = 0.0


# flat bed
bedTopography[:] = 0.0

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# melt
gridfile.variables['meltInput'][:] = 0.0
gridfile.variables['meltInput'][:] = 0.0002
# Use this line to only add a source term to the center cell - useful for debugging divergence
#gridfile.variables['meltInput'][0,r==0.0] = 2.0e-9  + 0.06/3.0e5 + (1.0e-7 * 10.0**5)   # From Ian Hewitt email

# velocity
velo = r*0.0
velo = (r-5000.0)*100.0/3.14e7/22500.0
velo[r<5000.0]=0.0
gridfile.variables['uReconstructX'][0,:,-1] = velo
gridfile.variables['uReconstructX'][0,thickness[0,:]==0.0,:] = 0.0

# IC on thickness
gridfile.variables['waterThickness'][0,:] = (r-5000.0)/r0 * 0.5 *0.0
##gridfile.variables['waterThickness'][0,r<5000.0] = 0.0
#gridfile.variables['waterThickness'][0,:] = 0.05

gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

