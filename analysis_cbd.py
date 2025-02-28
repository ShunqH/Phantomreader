import numpy as np
import index as ind
from PhantomDataReader import PhantomDataReader

class analysis_cbd(PhantomDataReader):
    def __init__(self, filename):
        super().__init__(filename)

    def move_to_center(self, masspts=None):
        if masspts is None:
            masspts = self.sink
        self.cal_cm(masspts)
        
        position_indices = {ind.x: 0, ind.y: 1, ind.z: 2}
        velocity_indices = {ind.vx: 0, ind.vy: 1, ind.vz: 2}

        for key in [ind.x, ind.y, ind.z, ind.vx, ind.vy, ind.vz]:
            if key in position_indices:  
                self.sink[:, key]    -= self.rcm[position_indices[key]]
                self.gas[:, key]     -= self.rcm[position_indices[key]]
                self.dead[:, key]    -= self.rcm[position_indices[key]]
                self.unknown[:, key] -= self.rcm[position_indices[key]]
            elif key in velocity_indices:  
                self.sink[:, key]    -= self.vcm[velocity_indices[key]]
                self.gas[:, key]     -= self.vcm[velocity_indices[key]]
                self.dead[:, key]    -= self.vcm[velocity_indices[key]]
                self.unknown[:, key] -= self.vcm[velocity_indices[key]]
        print("All particles are moved to center of mass coordinate. ")

    def cal_cm(self, masspts=None):
        if masspts is None:
            masspts = self.sink
        
        # Initialize sums for position and velocity
        totalmass = 0
        rcm = np.zeros(3)
        vcm = np.zeros(3)
        
        # Loop through mass points and accumulate the weighted sums
        for pt in masspts:
            mass = pt[ind.m]
            rcm += pt[ind.x:ind.z+1] * mass  # r[x, y, z] position weighted sum
            vcm += pt[ind.vx:ind.vz+1] * mass  # v[x, y, z] velocity weighted sum
            totalmass += mass

        # Compute the center of mass (position and velocity)
        self.rcm = rcm / totalmass
        self.vcm = vcm / totalmass

    def cross_section(self, term = None):
        if term is None:
            term = np.abs(self.gas[:, ind.z])<1
        self.disc = np.array([self.gas[i, :] for i in range(self.gas.shape[0]) if term[i]])