from openmm import *
from openmm.app import *
from openmm.unit import *
import mdtraj as md
from mdtraj.reporters import XTCReporter 
import re
import pdbfixer
import numpy as np
import pickle as pckl

#from openmmplumed import PlumedForce

################################# File Management ###################################

class metadp():
    def __init__(self,w,atoms,width,gridedges):
        self.weights=w
        self.atompairs=atoms
        self.gridedges=gridedges
        self.widths=width

        
class biasparams():
    def __init__(self,biasFactor,height,frequency,widthfactor,biasDir):
        self.biasFactor=biasFactor
        self.height=height
        self.frequency=frequency
        self.biasDir=biasDir
        self.widthfactor=widthfactor

def check_file(fname):
    '''
    Checks the existence of a file in the given pathway and gives a newfile name
    filename format  - <string_indentifier>_<int_identifier>.<extension>

    '''
    if os.path.isfile(fname):
        ident=fname.split('.')[0].split('_')
        fname=f'{ident[0]}_{int(ident[1])+1}.{fname.split(".")[1]}'
        fname=check_file(fname)
    return fname

def check_dir(fname):
    '''
    Checks the existence of a file in the given pathway and gives a newfile name
    filename format  - <string_indentifier>_<int_identifier>.<extension>

    '''
    if os.path.isdir(fname):
        ident=fname.split('_')
        ident[-1]=str(int(int(ident[-1])+1))
        fname="_".join(ident)
        fname=check_dir(fname)
    return fname

############################ Structure Editing Functions ############################

# Will be added into a class soon!

#structuture
def fix_pdb(file_n):
    """
    fixes the raw pdb from colabfold using pdbfixer.
    This needs to be performed to cleanup the pdb and to start simulation 
    Fixes performed: missing residues, missing atoms and missing Terminals
    """
    raw_pdb=file_n;

    # fixer instance
    fixer = pdbfixer.PDBFixer(raw_pdb)

    #finding and adding missing residues including terminals
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    outfile=check_file('fixed_0.pdb')
    PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile,'w'), keepIds=True)
    return fixer.positions, fixer.topology


def add_hydrogen(positions,topology,forcefield,write_file=True,):
    """
    Adds missing hydrogen to the pdb for a particular forcefield
    """
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield);
    if write_file:
        hydfile=check_file('fixedH_0.pdb')
        PDBFile.writeFile(modeller.topology, modeller.positions, open(hydfile, 'w'))
    return modeller.positions, modeller.topology

def solvate_me(positions,topology,forcefield,write_file=True,padding=1,\
               water_model='tip3p',positiveIon='Na+',negativeIion='Cl-'):
    '''
    Creates a box of solvent with 1 nm padding and neutral charge
    '''
    modeller = Modeller(topology, positions)
    modeller.addSolvent(forcefield, padding=padding*nanometers, model=water_model, neutralize=True, positiveIon=positiveIon, negativeIon=negativeIion)
    if write_file:
        solvfile=check_file('solvated_0.pdb')
        PDBFile.writeFile(modeller.topology, modeller.positions, open(solvfile, 'w'))
    return modeller.positions, modeller.topology


############################ Perform Dynamics using openMM ############################

# Will be added into a class soon!



def get_LangevinM_system(topology,forcefield,temp=300,dt=0.002,state=False):
    system = forcefield.createSystem(topology,nonbondedMethod=PME,nonbondedCutoff=1.2*nanometer,switchDistance=1.0*nanometer,\
                                     constraints=HBonds);
    integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picoseconds,dt*picoseconds);
    #platform = Platform.getPlatformByName('CUDA');
    #print(platform)
    properties = {'Precision': 'double'} #change if required after setting up openmm
    simulation = Simulation(topology, system, integrator)
    if state:
        simulation.loadState(state)
    return simulation


def get_NoseHoover_system(positions,topology,forcefield,temp=300,dt=0.002):
    system = forcefield.createSystem(topology,nonbondedMethod=PME,nonbondedCutoff=1.2*nanometer,switchDistance=1.0*nanometer,\
                                     constraints=HBonds);
    integrator = NoseHooverIntegrator(temp*kelvin, 1/picoseconds,dt*picoseconds);
    platform = Platform.getPlatformByName('CUDA');
    properties = {'Precision': 'double'} #change if required after setting up openmm
    simulation = Simulation(topology, system, integrator, platform)
    
    return simulation


def add_pos_res(positions,topology,simulation,k=1000):
    '''
    
    Adds an harmonic potential to the heavy atoms of the system(proteins) with an user defined force constant 'k'
    
    '''
    
    AA=['ALA','ASP','CYS','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ARG','PRO','GLN','ASN','SER','THR','VAL','TRP','TYR']
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2") # Harmonic potential for position restrain
    force.addGlobalParameter("k",k*kilocalories_per_mole/angstroms**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    index=0;
    for i, res in enumerate(topology.residues()):
        if res.name in AA:                              # Required to select only the protein atoms
            for at in res.atoms():
                if not re.search(r'H',at.name):         # All heavy Atoms -(exculdes Hydrogens)
                    force.addParticle(index,positions[index].value_in_unit(nanometers))
                index+=1;
                
    posres_sys=simulation.context.getSystem()                  # A gets System for a simulation instance
    posres_sys.addForce(force)                          # Modifies system with custom Force
    simulation.context.reinitialize()                          # initializes the simulation instance with the modified system
    
    return simulation

def add_hel_res(diheds,simulation,k=1):
    phase=[-1.1053785,-0.7255615]
    force = PeriodicTorsionForce()
    for i in [0,1]:
        for dihed in diheds[i]:
            force.addTorsion(dihed[0],dihed[1],dihed[2],dihed[3], 1, phase[i], k)
                
    posres_sys=simulation.context.getSystem()                  # A gets System for a simulation instance
    posres_sys.addForce(force)                          # Modifies system with custom Force
    simulation.context.reinitialize()                          # initializes the simulation instance with the modified system
    
    return simulation


def ener_Minimize(positions,simulation,tolerance=10,n_iter=1500,write_file=True):
    '''
    Energy minimization step for simulation object
    '''
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(tolerance=tolerance*kilojoule/mole,maxIterations=n_iter)
    minim_positions = simulation.context.getState(getPositions=True).getPositions()
    if write_file:
        minimfile=check_file('minim_0.pdb')
        pdb_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
        PDBFile.writeFile(simulation.topology, pdb_positions, open(minimfile, 'w'))
    return minim_positions,simulation



class DistanceReporter():
    def __init__(self, file, reportInterval, atomIndex1, atomIndex2):
        self._file = open(file, 'a')
        self._reportInterval = reportInterval
        self._atomIndex1 = atomIndex1
        self._atomIndex2 = atomIndex2

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions()
        distance = (positions[self._atomIndex1] - positions[self._atomIndex2]).value_in_unit(nanometer)
        magnitude = (distance**2).sum()**0.5  # Compute the magnitude of the distance vector
        self._file.write('%10.5f\n' % magnitude)

    def __del__(self):
        self._file.close()
        
class MultipleDistanceReporter():
    def __init__(self, file, reportInterval, atomPairs):
        
        print("yayyay")
        self._file = open(file, 'a')
        self._reportInterval = reportInterval
        self._atomPairs = atomPairs

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        positions = state.getPositions(asNumpy=True)  # Ensure you're getting positions as a numpy array
        distances = []
        for atomIndex1, atomIndex2 in self._atomPairs:
            distance_vector = np.array(positions[atomIndex1] - positions[atomIndex2])
            magnitude = np.linalg.norm(distance_vector)
            distances.append(magnitude)
        self._file.write('\t'.join(['%10.5f' % d for d in distances]) + '\n')

    def __del__(self):
        self._file.close()

#class md_wrapper():
#    def __init__(positions,simulation,nsteps,ens='NVT',run_type='equil',temp=300,pressure=1,metadparams=False,\
           #diheds=False,helixk=False,save_chkpt_file=True,outfreq=500,writeXTC=True,velocities=None,cont=False,restart=False,printdists=False,biasparams=False):

def run_MD(positions,simulation,nsteps,ens='NVT',run_type='equil',temp=300,pressure=1,metadparams=False,\
           diheds=False,helixk=False,save_chkpt_file=True,outfreq=500,writeXTC=True,velocities=None,cont=False,restart=False,printdists=False,biasparams=False,checksteps=False,checkstate=False):
    '''
    Function to run simulation for a particular ensemble. The function can also incorporate plumed file. 
    There is a restart option to start a simulation from checkpoint file
    There is also a continue option if multiple simulation is required to be run in a single script
    
    positions, simulation and nsteps are required parameters
    
    '''
    
    ## Add forces according to the ensemble or plumed file
    
    if ens=='NPT':
        simulation.context.getSystem().addForce(MonteCarloBarostat(pressure*bar, temp*kelvin)) # Pressure coupling
    
    if helixk:
        simulation=add_hel_res(diheds,simulation,helixk)
    
    if metadparams:                                       # Adding custom force via plumed or getting CVs via plumed
        
        w=metadparams.weights
        atoms=metadparams.atompairs
        gridedges=metadparams.gridedges
        width=metadparams.widths
        

        cvs=[CustomBondForce("c*r") for i in range(2)]
        [cv.addPerBondParameter("c") for cv in cvs]

        for i,atompair in enumerate(atoms):
            
            [cvs[j].addBond(atompair[0],atompair[1] , [w[j][i]]) for j in range(2)]
        
        
        CVs = [BiasVariable(cv, gridedges[i][0], gridedges[i][1], width[i]*biasparams.widthfactor/100, False) for i,cv in enumerate(cvs)]
        system=simulation.context.getSystem()
        meta = Metadynamics(system, CVs, 300.0*kelvin, biasparams.biasFactor, biasparams.height*kilojoules_per_mole, biasparams.frequency, biasDir=biasparams.biasDir, saveFrequency=biasparams.frequency)
#        meta = Metadynamics(system, CVs, 300.0*kelvin, 10.0, 1.5*kilojoules_per_mole, 1000, biasDir=".", saveFrequency=1000)
        #print(simulation.context.getState(getForces=True).getForces())
        simulation.context.reinitialize()
    ## Initialize the simulation with positions and velocity (if you are continuing to run a simulation )
    
    simulation.context.setPositions(positions)
    if cont:
        simulation.context.setVelocities(velocities)
    
    
    ## If restarting from checkpoint file and saving the same
    
    if restart:
        simulation.loadCheckpoint('chkptfile.chk')
    if save_chkpt_file:
        chkpt_freq=0.05*nsteps
        chkfile=check_file(f'{run_type}_0.chk')
        simulation.reporters.append(CheckpointReporter(chkfile, chkpt_freq))
    
    
    ## Append reporters for the simulation output and output files

    outfile=check_file(f'{run_type}_0.pdb')
    logfile=check_file(f'{run_type}_0.txt')          #output files
    xmlfile=check_file(f'{run_type}_0.state')
    

    simulation.reporters=[]
    outlog=open(logfile,'w')
    simulation.reporters.append(StateDataReporter(outlog, outfreq*2, step=True,potentialEnergy=True,kineticEnergy=True,separator='\t|\t',progress=True,speed=True,totalSteps=nsteps))
    
    if writeXTC:                                     #default is True
        outfname=check_file(f'{run_type}_0.xtc')
        topology=md.Topology.from_openmm(simulation.topology)
        XTC_information='protein'
        python_expression=topology.select_expression(XTC_information)
        req_indices=np.array(eval(python_expression))
        simulation.reporters.append(XTCReporter(outfname, outfreq, atomSubset=req_indices))
    
    ## RUN the simulation
    
    if printdists:
        distfile=check_file(f'{run_type}_0_dists.txt')
        simulation.reporters.append(MultipleDistanceReporter(distfile, 100, printdists))

     
    
    if metadparams:
        meta.step(simulation,nsteps)
    else:
        simulation.step(nsteps)
    
    ## Output PDB and simulation state
    
    pdb_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(simulation.topology, pdb_positions, open(outfile, 'w'))
    simulation.saveState(xmlfile)
    
    
    ## Get the positions and velocities that could be used to continue the simulation
    
    positions = simulation.context.getState(getPositions=True).getPositions()        #Note PBC condition note enforced -depends on long simulations
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    
    
    
    ## Remove the forces added to the simulation by plumed and monte carlo -- better solution will be updated (required for multiple runs in a single script)
    
    if ens=='NPT':
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)
    if helixk:
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)
    if metadparams:
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)  #from plumedforce    
    return positions,velocities,simulation


def run_restart(simulation,positions,velocities,nsteps,name,ens='NVT',temp=300,pressure=1,metadparams=False,\
           diheds=False,helixk=False,save_chkpt_file=True,outfreq=500,writeXTC=True,printdists=False,biasparams=False,checkstate=False):
    '''
    Function to run simulation for a particular ensemble. The function can also incorporate plumed file. 
    There is a restart option to start a simulation from checkpoint file
    There is also a continue option if multiple simulation is required to be run in a single script
    
    positions, simulation and nsteps are required parameters
    
    '''
    
    ## Add forces according to the ensemble or plumed file
    
    if ens=='NPT':
        simulation.context.getSystem().addForce(MonteCarloBarostat(pressure*bar, temp*kelvin)) # Pressure coupling
    
    if helixk:
        simulation=add_hel_res(diheds,simulation,helixk)
    
    if metadparams:                                       # Adding custom force via plumed or getting CVs via plumed
        
        w=metadparams.weights
        atoms=metadparams.atompairs
        gridedges=metadparams.gridedges
        width=metadparams.widths*biasparams.widthfactor/100
        

        cvs=[CustomBondForce("c*r") for i in range(2)]
        [cv.addPerBondParameter("c") for cv in cvs]

        for i,atompair in enumerate(atoms):
            
            [cvs[j].addBond(atompair[0],atompair[1] , [w[j][i]]) for j in range(2)]
        
        
        CVs = [BiasVariable(cv, gridedges[i][0], gridedges[i][1], width[i], False) for i,cv in enumerate(cvs)]
        system=simulation.context.getSystem()
        if biasparams:
            meta = Metadynamics(system, CVs, 300.0*kelvin, biasparams.biasFactor, biasparams.height*kilojoules_per_mole, biasparams.frequency, biasDir=biasparams.biasDir, saveFrequency=1000)

        else:
            meta = Metadynamics(system, CVs, 300.0*kelvin, 10.0, 1.5*kilojoules_per_mole, 1000, biasDir=".", saveFrequency=1000)
        #print(simulation.context.getState(getForces=True).getForces())
        
        simulation.context.reinitialize()
    ## Initialize the simulation with positions and velocity (if you are continuing to run a simulation )
    simulation.context.setPositions(positions)

    simulation.context.setVelocities(velocities)
    if save_chkpt_file:
        chkpt_freq=0.05*nsteps
        simulation.reporters.append(CheckpointReporter('chkptfile.chk', chkpt_freq))
    
    ## Append reporters for the simulation output and output files

    outfile='%s.pdb'%name
    logfile='%s.txt'%name        #output files
    xmlfile='%s.state'%name
    

    simulation.reporters=[]
    outlog=open(logfile,'a')
    simulation.reporters.append(StateDataReporter(outlog, outfreq*2, step=True,potentialEnergy=True,kineticEnergy=True,separator='\t|\t',progress=True,speed=True,totalSteps=nsteps))
    
    if writeXTC:                                     #default is True
        outfname='%s.xtc'%name
        topology=md.Topology.from_openmm(simulation.topology)
        XTC_information='protein'
        python_expression=topology.select_expression(XTC_information)
        req_indices=np.array(eval(python_expression))
        simulation.reporters.append(XTCReporter(outfname, outfreq, atomSubset=req_indices))
    
    ## RUN the simulation
    
    if printdists:
        distfile='%s_dists.txt'%name
        simulation.reporters.append(MultipleDistanceReporter(distfile, 100, printdists))

    
    if metadparams:
        meta.step(simulation,nsteps)
    else:
        simulation.step(nsteps)
    
    ## Output PDB and simulation state
    
    pdb_positions = simulation.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
    PDBFile.writeFile(simulation.topology, pdb_positions, open(outfile, 'w'))
    simulation.saveState(xmlfile)
    
    
    ## Get the positions and velocities that could be used to continue the simulation
    
    positions = simulation.context.getState(getPositions=True).getPositions()        #Note PBC condition note enforced -depends on long simulations
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    
    
    
    ## Remove the forces added to the simulation by plumed and monte carlo -- better solution will be updated (required for multiple runs in a single script)
    
    if ens=='NPT':
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)
    if helixk:
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)
    if metadparams:
        simulation.context.getSystem().removeForce(simulation.context.getSystem().getNumForces()-1)  #from plumedforce    
    return positions,velocities,simulation