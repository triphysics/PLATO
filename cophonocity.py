import numpy as np

data=np.loadtxt("outfile.phonon_dos")
qp=data.shape[0]
freq=data[:,0]

gw=data[:,1]

gw_Mg=data[:,2]
gw_O=data[:,3]

order=1
moment = 0.0

fmin=0.01
fmax=25
norm = 0.0
norm_Mg= 0.0
norm_O = 0.0

moment_Mg = 0.0
moment_O= 0.0

for i in range(0,qp):
    if fmin < freq[i] and freq[i] < fmax :
        moment += freq[i]**order *gw[i]
        norm +=  gw[i]  

        moment_Mg += freq[i]**order *gw_Mg[i]
        norm_Mg +=  gw_Mg[i]  

        moment_O += freq[i]**order *gw_O[i]
        norm_O +=  gw_Mg[i]

print(moment/norm)
print(moment_Mg/norm_Mg)
print(moment_O/norm_O)
