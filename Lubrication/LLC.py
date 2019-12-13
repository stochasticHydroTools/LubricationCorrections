'''
Small class to handle linked list cell (LLC [but not limited liability coporation]) data structure
'''
import numpy as np
import copy
import sys

class LLC(object):
  '''
  Small class to linked list cell (LLC)
  '''  
  def __init__(self, locations):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 3
    self.r_locs = locations
    self.numberOfParticlesInCell = None
    self.particleIndexInCell = None
    self.cellParticle = None
    self.numberOfCellsX = None
    self.numberOfCellsY = None
    self.numberOfCellsZ = None
    self.numberOfCells = None
    self.maxNumberOfParticlesInCell = 100
    print "NO PERIODIC IMPLEMENTED BUT IS EASY TO DO"

  
  def make_particle_linked_lists(self, cutoff, maxDx):
    '''
    make linked list data structures
    '''
    r_locs = self.r_locs
    numberOfParticles = len(r_locs)
    
    lx = abs(np.amax(r_locs[:,0]))
    ly = abs(np.amax(r_locs[:,1]))
    lz = 10*abs(np.amax(r_locs[:,2]))
    
    cutoff_max = 2*cutoff + 2*maxDx
    self.numberOfCellsX = int(lx / cutoff_max)   
    self.numberOfCellsY = int(ly / cutoff_max)
    self.numberOfCellsZ = int(lz / cutoff_max) 
    
    self.numberOfCellsX = max(self.numberOfCellsX,3)
    self.numberOfCellsY = max(self.numberOfCellsY,3)
    self.numberOfCellsZ = max(self.numberOfCellsZ,3)
    
    numberOfCells = self.numberOfCellsX * self.numberOfCellsY * self.numberOfCellsZ
    dxCell = lx / self.numberOfCellsX
    dyCell = ly / self.numberOfCellsY
    dzCell = lz / self.numberOfCellsZ
    
    self.numberOfParticlesInCell = np.zeros(numberOfCells)
    self.particleIndexInCell = np.empty(numberOfCells * self.maxNumberOfParticlesInCell)
    self.cellParticle = np.empty(numberOfParticles)
    
    for i in range(numberOfParticles):
      x = r_locs[i,0]
      y = r_locs[i,1]
      z = r_locs[i,2]
      kx = int(x / dxCell) 
      ky = int(y / dyCell) 
      kz = int(z / dzCell) 
      cell = self.numberOfCellsX*self.numberOfCellsY*kz + self.numberOfCellsX*ky + kx
      if(cell >= numberOfCells):
	print "problem particle: " + str(i) + " cell " + str(cell) + " coors " +  str(x) + " " + str(y) + " " + str(z)
      self.numberOfParticlesInCell[cell] += 1 
      if(self.numberOfParticlesInCell[cell] >= self.maxNumberOfParticlesInCell):
	print "num particles in cell " + str(cell) + "exceeds max" 
      perm = int((cell)*self.maxNumberOfParticlesInCell + self.numberOfParticlesInCell[cell])
      self.particleIndexInCell[perm] = i
      self.cellParticle[i] = cell
      
  def query_particle(self, i):
    icell = self.cellParticle[i]

    kz = int((icell) / (self.numberOfCellsX*self.numberOfCellsY))
    ky = int((icell - (self.numberOfCellsX*self.numberOfCellsY*kz)) / self.numberOfCellsX)
    kx = int(icell - self.numberOfCellsX*self.numberOfCellsY*kz - self.numberOfCellsX*ky)
    
    index = []
    
    for ikx in [kx-1, kx, kx+1]:
      iikx = ikx
      if(iikx == -1 or iikx == self.numberOfCellsX):
	#iikx = (self.numberOfCellsX-1)*(iikx == -1)
	continue
      for iky in [ky-1, ky, ky+1]:
	iiky = iky
	if(iiky == -1 or iiky == self.numberOfCellsY):
	  #iiky = (self.numberOfCellsY-1)*(iiky == -1)
	  continue
	for ikz in [kz-1, kz, kz+1]:
	  iikz = ikz
	  if(iikz == -1 or iikz == self.numberOfCellsZ):
	    #iikz = (self.numberOfCellsZ-1)*(iikz == -1)
	    continue 
	  cell = int(self.numberOfCellsX*self.numberOfCellsY*iikz + self.numberOfCellsX*iiky + iikx)
	  if(cell >= self.numberOfCellsX*self.numberOfCellsY*self.numberOfCellsZ):
	    print "problem particle: " + str(i) + " cell " + str(cell) + " coors " +  str(iikx) + " " + str(iiky) + " " + str(iikz)
          # Find neighbors of particle i in cell "cell"
          num_part = int(self.numberOfParticlesInCell[cell])
          
          for j in range(num_part):
	    perm = int((cell)*self.maxNumberOfParticlesInCell + j + 1)
            index.append(self.particleIndexInCell[perm])
    return np.array(index)