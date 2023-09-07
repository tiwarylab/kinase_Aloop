import numpy as np
from openmmplumed import PlumedForce
from openmm.app import *
from openmm import *
from openmm.unit import *
import pdbfixer
import openmm_utils as op            # openmm functions
import sys 
import os
import mdtraj as md
import kinaseCVs as kcv
from kinaseCVs import resids_DDR1

oss=os.system


if __name__=='__main__':
    '''
    Python file to perform Molecular dynamics using openmm - basic example to use openmm_utils package
    
        (A better script will be updated)
    
    Args:
    
    Input:
     -file            : PDB file to start 
     -nsteps          : Number simulation steps (default=25000000)
     -temp            : Temperature of simulation (default=300K)
     -pressure        : Pressure of simulation (default=1bar)

    '''
    
    if '-name' in sys.argv:
        name = sys.argv[sys.argv.index('-name')+1]
        flag=1
    else:
        flag=0
        
    if '-nsteps' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-nsteps')+1])

    elif 'tns' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-tns')+1])*500000
        flag=1
    else:
        nsteps=5000000
    if '-temp' in sys.argv:
        temp = int(sys.argv[sys.argv.index('-temp')+1])   #Integer value
    else:
        temp = 300              # 300K temperature

        
        
    if flag:
            
             
            if '-path' in sys.argv:
                os.chdir(sys.argv[sys.argv.index('-path')+1])
            pdb=PDBFile("%s.pdb"%name)
          
            topology=pdb.topology
            positions=pdb.positions
            
            ff = ForceField('amber99sbildn.xml', 'tip3p.xml')
            simulation=op.get_LangevinM_system(topology,ff,temp,0.002)
            
            top=md.Topology.from_openmm(simulation.topology)
            traj=md.load("%s.pdb"%name)
            
            from openmm import XmlSerializer

            with open('%s.state'%name, 'r') as f:
                state_xml = f.read()

            state = XmlSerializer.deserialize(state_xml)
            velocities = state.getVelocities()
            classk=kcv.kinase_cvs("WT",traj.top,resids_DDR1,traj)
            atoms=classk.getatoms()
            # NPT - no restrain
            positions,velocities,simulation=op.run_restart(simulation,positions,velocities,nsteps,name,ens='NPT',printdists=atoms)                            





