import numpy as np
import index as ind 

class PhantomDataReader:
    def __init__(self, filename):
        if self._checkfilename(filename):
            self.filename = filename
            self.time = None
            self.time_unit = None
            self.npart = {}
            self.units = []
            self.columns = []
            self._parse_file()
            print(self)
        else:
            print("Please input ascii file, try 'splash to ascii disc_00000'. ")
    
    def _parse_file(self):
        self._parse_header()
        self._parse_data()
    
    def _parse_header(self):
        with open(self.filename, 'r') as file_handle:
            
            for line in file_handle:
                if not line.startswith('#'):
                    break  
                    
                line = line.strip('# \n')
                
                if line.startswith('time:'):
                    time_values = next(file_handle).strip('#\n').split()
                    self.time = float(time_values[0])
                    self.time_unit = float(time_values[1])
                    
                elif line.startswith('npart:'):
                    particle_counts = next(file_handle).strip('#\n').split()
                    particle_types = [
                        'gas', 'dust (old)', 'sink', 'ghost', 'star',
                        'dark matter', 'bulge', 'type 8', 'type 9', 'type 10',
                        'type 11', 'type 12', 'type 13', 'type 14', 'type 15',
                        'type 16', 'type 17', 'type 18', 'type 19', 'type 20',
                        'type 21', 'type 22', 'type 23', 'unknown/dead'
                    ]
                    self.npart = dict(zip(particle_types, map(int, particle_counts)))
                    
                elif line.startswith('units:'):
                    next(file_handle)
                    unit_values = next(file_handle).strip('[]#\n').split()
                    self.units = list(map(str, unit_values))
                    
                elif line.startswith('x [au]'):
                    self.columns = ["x", "y", "z", "particle mass", "h", "density", "v_x", "v_y", "v_z", "div_v", "dt", "type"]
    
    def _parse_data(self):
        data = np.loadtxt(self.filename, skiprows=12)
        self.gas = data[data[:, -1] == 1]
        self.sink = data[data[:, -1] == 3]
        self.dead = data[data[:, -1] == 24]
        self.unknown = data[~np.isin(data[:, -1], [1, 3, 24])]

    def __repr__(self):
        return (f"PhantomDataReader(filename='{self.filename}')\n"
                f"Time: {self.time} yr\n"
                f"Particles: sink({self.npart["sink"]}), gas({self.npart["gas"]}), dead({self.npart["unknown/dead"]}), else({len(self.unknown)})\n"
                f"Data for each particle: {self.columns}")
    
    def spherical(self, xc=0, yc=0, zc=0):
        if ind.R == -1 or (self.gas.shape[1]<13): 
            R = self.dist(self.gas[:,ind.x],self.gas[:,ind.y],self.gas[:,ind.z], xc,yc,zc)
            phi = np.arctan2(self.gas[:,ind.y]-yc, self.gas[:,ind.x]-xc)
            theta = np.arcsin((self.gas[:,ind.z]-zc)/self.R)
            self.gas = np.column_stack((self.gas, R, phi, theta))
            ind.R = self.gas.shape[1] - 3
            ind.phi = self.gas.shape[1] - 2
            ind.theta = self.gas.shape[1] - 1
            # self.vR = 
            print("Now you can use '*.gas[:, ind.R]', '*.gas[:, ind.phi]', and '*.gas[:, ind.theta]'. ")
        else:
            print("Already done last time, try '*.gas[:, ind.R]', '*.gas[:, ind.phi]', and '*.gas[:, ind.theta]'. ")
    
    def cylinder(self, xc=0, yc=0, zc=0):
        if (ind.r == -1) or (self.gas.shape[1]<13): 
            r = self.dist(self.gas[:,ind.x],self.gas[:,ind.y],0, xc,yc,zc)
            phi = np.arctan2(self.gas[:,ind.y]-yc, self.gas[:,ind.x]-xc)
            vr = (self.gas[:,ind.x]*self.gas[:,ind.vx] + self.gas[:,ind.y] * self.gas[:,ind.vy])/r
            vphi = (-self.gas[:,ind.x]*self.gas[:,ind.vy] + self.gas[:,ind.y] * self.gas[:,ind.vx])/r
            # z_cyl = self.gas[:,ind.z] - zc 
            self.gas = np.column_stack((self.gas, r, phi, vr, vphi))
            ind.r = self.gas.shape[1] - 4
            ind.phi = self.gas.shape[1] - 3
            ind.vr = self.gas.shape[1] - 2
            ind.vphi = self.gas.shape[1] - 1
            # ind.z_cyl = self.gas.shape[1] - 1
            print("Now you can use '*.gas[:, ind.r]', '*.gas[:, ind.phi]', '*.gas[:, ind.vr]', and '*.gas[:, ind.vphi]'. ")
        else: 
            print("Already done last time, try '*.gas[:, ind.r]', '*.gas[:, ind.phi]', '*.gas[:, ind.vr]', and '*.gas[:, ind.vphi]'. ")
    
    def _checkfilename(self, filename):
        if filename[-5:]=="ascii":
            return True 

    @staticmethod
    def dist(x1, y1, z1, x2=0, y2=0, z2=0):
        return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
