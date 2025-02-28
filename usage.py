from PhantomDataReader import PhantomDataReader
import index as ind 

Path = "/Users/shunq/Documents/cluster/wind2cbd/run39/"     # path to your file
filename = Path + "disc_00100.ascii"                        # your filename 

frame = PhantomDataReader(filename)     # read data into frame, it will also print some iformation
print("\n"+"The mass of the first sink particle: ")
print(frame.sink[0, ind.ptm])           # print the mass (ind.ptm) of the first (0) sink particle     
print("\n"+"The x position for the first 5 gas particles: ")
print(frame.gas[:5, ind.x])             # print x position (ind.x) of the 5 (:5) gas particles   
   
