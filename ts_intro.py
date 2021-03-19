import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.image as mpimg
from mpl_toolkits import mplot3d
from Data_Point import Data_Point
# read in data
intensity  = np.loadtxt('intensity.dat')	# 2D array of CCD counts data
wavelength = np.loadtxt('lambda.dat')		# 2D array of wavelength data
radius     = np.loadtxt('radius.dat')		# 1D array of major radii
angle      = np.loadtxt('angle.dat')		# 1D array of scattering angles


#imgplot = plt.imshow(intensity)
#plt.colorbar()
#plt.show()

N_first_D = len(intensity)
N_sec_D = len(intensity[0])

first_array = np.zeros([N_first_D])
sec_array = np.zeros([N_sec_D])

for i in range(N_first_D):
    for j in range(N_sec_D):
        first_array[i]+=intensity[i][j]
for i in range(N_sec_D):
    for j in range(N_first_D):
        sec_array[i]+=intensity[j][i]
        
first_array_x = np.arange(0, N_first_D, 1)
sec_array_x = np.arange(0, N_sec_D, 1)

# plt.plot(first_array_x, first_array, label='Flattened on axis 1')
# plt.plot(sec_array_x, sec_array, label='Flattened on axis 2')
# plt.legend()
# plt.xlabel('Cell Number')
# plt.ylabel('Integrated Intensity')
# plt.show()
#print(wavelength)
#imgplot = plt.imshow(wavelength, origin='lower')
#plt.colorbar()
#plt.show()
#print(len(wavelength))
#print(len(wavelength[0]))

cm = plt.cm.get_cmap('RdYlBu')
xy = range(20)
#print(xy)
z = xy

data_array = [[None]*N_sec_D]*N_first_D
#figure = plt.figure()
#axis = mplot3d.Axes3D(figure)
xdata = []
ydata = []
zdata = []
for i in range(N_first_D):
    for j in range(N_sec_D):
        ij_data_point = Data_Point(intensity[i][j], j, i, wavelength[i][j])
        data_array[i][j] = ij_data_point
        xdata.append(ij_data_point.wavelength_angstrom)
        ydata.append(ij_data_point.cell_y_position)
        zdata.append(ij_data_point.intensity)
#axis.scatter(xdata, ydata, zdata)
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(xdata, ydata, c=zdata, vmin=min(zdata), vmax=max(zdata))#, s=35, cmap=cm)
plt.colorbar(sc)
plt.xlim([550, 950])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Position [arb.]')
plt.show()
#plt.show()
count = 0
for i in zdata:
    count+= i
print(count)