#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 16:07:54 2022

@author: bvani
"""
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md
import random

resids_DDR1={'HRDxN': 172, 'DFGAsp': 185, 'ChelE': 73, 'X2': 80, 'X3': 164, 'sbridgeR': 171, 'sbridgeK': 56, 'PloopN1': 17, 'PloopC1': 22, 'DFGPhe': 186, 'DFGGly': 187, 'ChelX': 77, 'X4': 190, 'PloopN2': 18, 'PloopC2': 21}

#topology=md.load("/home/bvani/static_analysis/nmr/6xr6/6xr6_01.gro").topology


class distance():
    def __init__(self, name, description, resids, atomids, top, traj=0):
        self.name=name
        self.description=description
        self.resids=resids
        self.atomids=atomids
        
        self.select0=top.select("resid %i and name %s"%(resids[0],atomids[0]))[0]
        self.select1=top.select("resid %i and name %s"%(resids[1],atomids[1]))[0]
        if traj:
            self.traj=self.compute(traj)
        else:
            self.traj=0
            
    def compute(self,traj):
        self.traj=md.compute_distances(traj, np.asarray(((self.select0,self.select1))).reshape(1,2),periodic=False).reshape(len(traj))
        return self.traj

    def printplumed(self,f):
        f.write("\n%s: DISTANCE ATOMS=%i,%i \n"%(self.__class__.__name__,self.select0+1,self.select1+1))


#%%Levy distances
class levy1(distance):
    def __init__(self,resdicts,top,traj=0):
        
        des="Levy D1, the Cαatom distance between the Asn of the HRDxxxxN motif (the first Asn residue that follows the HRD motif) and Phe of the DFG motif"
        name="DFG-HRDxN" #NF
        
        resids=[resdicts["HRDxN"],resdicts["DFGPhe"]]
        
        atomids=["CA","CA"]
        super().__init__(name, des, resids, atomids, top, traj)

class levy2(distance):
    def __init__(self,resdicts,top,traj=0):
        des="Levy D2, the Cαatom distance between the conserved Glu belonging to the αC-helix, and Phe of the DFG motif"
        name="DFG-Chelix"#FE
        resids=[resdicts["ChelE"],resdicts["DFGPhe"]]
        atomids=["CA","CA"]
        super().__init__(name, des, resids, atomids, top, traj)

#%% Gervasio distances

class gervasio1(distance):
    def __init__(self,resdicts,top,traj=0):
        des="Salt bridge with lysinewhen active DFG-in"
        name="DFG-N lobe"
        resids=[resdicts["sbridgeK"],resdicts["DFGAsp"]]
        atomids=["CB","CB"]
        super().__init__(name, des, resids, atomids, top, traj)

class gervasio2(distance):
    def __init__(self,resdicts,top,traj=0):
        des="Aloop-Clobe (Gervasio)"
        name="Aloop-Clobe"
        resids=[resdicts["X3"],resdicts["X4"]]
        atomids=["CB","CB"]
        super().__init__(name, des, resids, atomids, top, traj)

#%% Benoit

class benoit1(distance):
    def __init__(self,resdicts,top,traj=0):
        des="First (outer) P-loop distance"
        name="P-loop1"
        resids=[resdicts["PloopN1"],resdicts["PloopC1"]]
        atomids=["CA","CA"]
        super().__init__(name, des, resids, atomids, top, traj)
        
class benoit2(distance):
    def __init__(self,resdicts,top,traj=0):
        des="Second (inner) P-loop distance"
        name="P-loop2"
        resids=[resdicts["PloopN2"],resdicts["PloopC2"]]
        atomids=["CA","CA"]
        super().__init__(name, des, resids, atomids, top, traj)
        
class benoit3(distance):
    def __init__(self,resdicts,top,traj=0):
        des="C gamma D(FG) and C alpha-N, DFG config"
        name="DFG1"
        resids=[resdicts["DFGAsp"],resdicts["HRDxN"]]
        atomids=["CG","CA"]
        super().__init__(name, des, resids, atomids, top, traj)
        
class benoit4(distance):
    def __init__(self,resdicts,top,traj=0):
        des="C gamma-D(F)G and C alpha-I, DFG config"
        name="DFG2"
        resids=[resdicts["DFGPhe"],resdicts["X2"]]
        atomids=["CG","CA"]
        super().__init__(name, des, resids, atomids, top, traj)    

class benoit5(distance):
    def __init__(self,resdicts,top,traj=0):
        des="C eta-R and C delta-E, alphaC helix, salt bridge"
        name="Chelix salt bridge1"
        resids=[resdicts["sbridgeR"],resdicts["ChelE"]]
        atomids=["CZ","CD"]
        super().__init__(name, des, resids, atomids, top, traj)    

class benoit6(distance):
    def __init__(self,resdicts,top,traj=0):
        des="N eta-K and C delta-E, alphaC helix, salt bridge"
        name="Chelix salt bridge2"
        resids=[resdicts["sbridgeK"],resdicts["ChelE"]]
        atomids=["NZ","CD"]
        super().__init__(name, des, resids, atomids, top, traj)    

class benoit7(distance):
    def __init__(self,resdicts,top,traj=0):
        des="O DF(G) and N R, A-loop N-terminal helix"
        name="Aloop"
        resids=[resdicts["DFGGly"],resdicts["sbridgeR"]]
        atomids=["O","N"]
        super().__init__(name, des, resids, atomids, top, traj)    

#%%
        
class dunbrack1(distance):
    def __init__(self,resdicts,top,traj=0):
        des="dist(alphaChelix-Glu(+4)-Calpha, DFG-Phe-Ceta)"
        name="DFG-Clobe"
        resids=[resdicts["ChelX"],resdicts["DFGPhe"]]
        atomids=["CA","CZ"]
        super().__init__(name, des, resids, atomids, top, traj)    

class dunbrack2(distance):
    def __init__(self,resdicts,top,traj=0):
        des="dist(beta3-Lys-Calpha, DFG-Phe-Ceta)"
        name="DFG-Nlobe"
        resids=[resdicts["sbridgeK"],resdicts["DFGPhe"]]
        atomids=["CA","CZ"]
        super().__init__(name, des, resids, atomids, top, traj)    

class dunbrack3(distance):
    def __init__(self,resdicts,top,traj=0):
        des="C-helix–Glu–Cbeta and beta3–Lys–Cbeta"
        name="Chelix salt-bridge"
        resids=[resdicts["ChelE"],resdicts["sbridgeK"]]
        atomids=["CB","CB"]
        super().__init__(name, des, resids, atomids, top, traj)    
        
#%%
        
class kinase_cvs():
    def __init__(self,name,topology,resids,traj=0):
        self.name=name
        self.rdict=resids
        
        self.rdict.update({"DFGPhe":resids["DFGAsp"]+1,"DFGGly":resids["DFGAsp"]+2, "ChelX":resids["ChelE"]+4, "X4":resids["DFGAsp"]+5, "PloopN2":resids["PloopN1"]+1, "PloopC2":resids["PloopC1"]-1 })
        
        self.top=topology
        
        self.levy1=levy1(self.rdict,self.top,traj)
        self.levy2=levy2(self.rdict,self.top,traj)
        self.gervasio1=gervasio1(self.rdict,self.top,traj)
        self.gervasio2=gervasio2(self.rdict,self.top,traj)
        self.benoit1=benoit1(self.rdict,self.top,traj)
        self.benoit2=benoit2(self.rdict,self.top,traj)
        self.benoit3=benoit3(self.rdict,self.top,traj)
        self.benoit4=benoit4(self.rdict,self.top,traj)
        self.benoit5=benoit5(self.rdict,self.top,traj)
        self.benoit6=benoit6(self.rdict,self.top,traj)
        self.benoit7=benoit7(self.rdict,self.top,traj)
        self.dunbrack1=dunbrack1(self.rdict,self.top,traj)
        self.dunbrack2=dunbrack2(self.rdict,self.top,traj)
        self.dunbrack3=dunbrack3(self.rdict,self.top,traj)
        self.allcvs=[self.levy1,self.levy2,self.gervasio1,self.gervasio2,self.benoit1,self.benoit2,self.benoit3,self.benoit4,self.benoit5,self.benoit6,self.benoit7,self.dunbrack1,self.dunbrack2,self.dunbrack3]
      
    def getnames(self):
        return [a.name for a in self.allcvs]
    
    def gettraj(self):
        return [a.traj for a in self.allcvs]
    
    def traj(self,i):
        return self.allcvs[i].traj
    
    def name(self,i):
        return self.allcvs[i].name
    
    def des(self,i):
        return self.allcvs[i].description


    def printplumed(self):
        f=open("%s_CVs.dat"%self.name,"w+")
        
        f.write("\nWHOLEMOLECULES ENTITY0=1-%i \n"%(self.top.n_atoms+1))
        [x.printplumed(f) for x in self.allcvs]
        listcvs=",".join([x.__class__.__name__ for x in self.allcvs])
        f.write('\nPRINT ARG=%s STRIDE=500 FILE=COLVAR_%s'%(listcvs,self.name))
        f.write('\nDUMPMASSCHARGE FILE=mc_%s'%self.name)

    def make_biased_plumed(self,spibdt,weights,height,biasfactor,width1,width2,gridmin1,gridmin2,gridmax1,gridmax2,temperature):
        f=open("%s_metad_%i.dat"%(self.name,spibdt),"w+")
        
        f.write("\nWHOLEMOLECULES ENTITY0=1-%i \n"%(self.top.n_atoms+1))
        [x.printplumed(f) for x in self.allcvs]
        listcvs=",".join([x.__class__.__name__ for x in self.allcvs])
        
        w0=",".join([str(weights[0][i]) for i in range (len(weights[0]))])
        w1=",".join([str(weights[1][i]) for i in range (len(weights[1]))])
        f.write("\nsigma1: COMBINE ARG=%s COEFFICIENTS=%s PERIODIC=NO"%(listcvs,w0))
        f.write("\nsigma2: COMBINE ARG=%s COEFFICIENTS=%s PERIODIC=NO"%(listcvs,w1))

        f.write("\nMETAD ...\n \
          LABEL=metad\n \
          ARG=sigma1,sigma2\n \
          PACE=1000 HEIGHT=%f TEMP=%i\n \
          BIASFACTOR=%i\n \
          SIGMA=%f,%f\n \
          FILE=HILLS GRID_MIN=%f,%f GRID_MAX=%f,%f GRID_BIN=200,200\n \
          CALC_RCT RCT_USTRIDE=500\n \
          ... METAD\n"%(height,temperature,biasfactor,width1,width2,gridmin1,gridmin2,gridmax1,gridmax2))
      
        f.write('\nDUMPMASSCHARGE FILE=mc_%s'%self.name)

        f.write("\nPRINT ARG=%s,sigma1,sigma2,metad.rbias STRIDE=500 FILE=COLVAR_biased.dat"%listcvs)

        f.close()

    def getatoms(self):
      atoms=[[x.select0,x.select1] for x in self.allcvs]
      return atoms
        
def RegSpaceClustering(z, min_dist, max_centers=200, batch_size=100,randomseed=0,periodicity=0):
    '''Regular space clustering.
    Args:
        data: ndarray containing (n,d)-shaped float data
        max_centers: the maximum number of cluster centers to be determined, integer greater than 0 required
        min_dist: the minimal distances between cluster centers
    '''
    random.seed(5)
    num_observations, d = z.shape
    p = np.hstack((0,np.random.RandomState(seed=randomseed).permutation(num_observations-1)+1))
    data = z[p]
    center_list = data[0, :].copy().reshape(d,1)
    centerids=[p[0]+1]
    i = 1
    while i < num_observations:
        x_active = data[i:i+batch_size, :]
        differences=np.abs(np.expand_dims(center_list.T,0) - np.expand_dims(x_active,1))
        differences=np.max(np.stack((differences,periodicity-differences)),axis=0)
        distances = np.sqrt((np.square(differences)).sum(axis=-1))
        indice = tuple(np.nonzero(np.all(distances > min_dist, axis=-1))[0])
        if len(indice) > 0:
            # the first element will be used
            #print(center_list.shape,x_active.shape,x_active[indice[0]].shape)
            center_list = np.hstack((center_list, x_active[indice[0]].reshape(d,1)))
            centerids.append(p[i+indice[0]]+1)
            i += indice[0]
        else:
            i += batch_size
        if len(centerids) >= max_centers:
            print("%i centers: Exceeded the maximum number of cluster centers!\n"%len(centerids))
            print("Please increase dmin!\n")
            raise ValueError
    print("Found %i centers!"%len(centerids))
    return center_list,centerids