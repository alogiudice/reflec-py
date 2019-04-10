import numpy as np
import matplotlib.pyplot as plt
outfile = "XRD_Real1396.dat"
fout = open(outfile, "w+")
with open("/home/agostina/Muestras Si/Muestras Real/RX/data_XRD_Real1396.dat") as f:
     data = f.readlines()
     for line in data:
         values = line.split()

for i in range(0, len(values)):
    values[i] = float(values[i])

theta = np.linspace(7.8370, 23.3320, len(values))
for i in range(0, len(values)):
    fout.write("%f %f\n" % (theta[i], values[i]))

fout.close()
fig = plt.figure()
plt.scatter(theta,values, s = 1)
plt.show()


