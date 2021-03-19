import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.image as mpimg
from mpl_toolkits import mplot3d
from Data_Point import Data_Point
from scipy import interpolate, optimize
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

#construct a data array
data_array = []
for i in range(N_first_D):
    data_array_line = []
    for j in range(N_sec_D):
        ij_data_point = Data_Point(intensity[i][j], j, i, wavelength[i][j], radius[i], angle[i])
        data_array_line.append(ij_data_point)
    data_array.append(data_array_line)
    

def get_slice(datas, index):
    xs = []
    ys = []
    for data in datas[index]:
        xs.append(data.wavelength_angstrom)
        ys.append(data.intensity)
    return xs, ys

def remove_filter_section(xs, ys, filter_range):
    new_xs = []
    new_ys = []
    for i in range(len(xs)):
        x = xs[i]
        y = ys[i]
        if x<filter_range[0] or x>filter_range[1]:
            new_xs.append(x)
            new_ys.append(y)
    return new_xs, new_ys


def remove_lines(xs, ys, wavelength_min, wavelength_max, is_plot_lines=True):
    Balmer_wavelength = 656.279 # in nm
    laser_wavelength = 694.3  # in nm
    BUFFER = 12 #in nm
    if is_plot_lines:
        plt.axvline(x=Balmer_wavelength, color= 'r', linestyle ='--', label='Balmer-alpha Line')
        plt.axvline(x=Balmer_wavelength+BUFFER, color= 'r', linestyle =(0,(1,10)), label='Balmer-alpha trim')
        plt.axvline(x=Balmer_wavelength-BUFFER, color= 'r', linestyle =(0,(1,10)), label='Balmer-alpha trim')
        plt.axvline(x=laser_wavelength, color= 'b', linestyle ='--', label='Ruby Laser Line')
        plt.axvline(x=laser_wavelength+BUFFER, color= 'b', linestyle =(0,(1,10)), label='Ruby Laser trim')
        plt.axvline(x=laser_wavelength-BUFFER, color= 'b', linestyle =(0,(1,10)), label='Ruby Laser trim')
    xs, ys = remove_filter_section(xs, ys,
                                   [Balmer_wavelength-BUFFER,
                                    Balmer_wavelength+BUFFER])
    xs, ys = remove_filter_section(xs, ys,
                                   [laser_wavelength-BUFFER,
                                    laser_wavelength+BUFFER])
    xs, ys = remove_filter_section(xs, ys,[0, wavelength_min])
    xs, ys = remove_filter_section(xs, ys,[wavelength_max, wavelength_min*2])
    return xs, ys

def gaussian(fit_vals,x):
    #note linear under log
    return fit_vals[0]*np.exp(-pow(x-fit_vals[1], 2)/(2*pow(fit_vals[2],2)))
    
def perform_fit(xs, ys):
    abc0=[2000, 694.3, 10]
    fit_vals, cov_vals = optimize.curve_fit(
                        lambda x,a,b,c: gaussian([a,b,c],x),
                        xs,  ys,  p0=(abc0[0], abc0[1],abc0[2]),
                        maxfev=100000)
    
    errors = [np.sqrt(abs(cov_vals[0][0])),
              np.sqrt(abs(cov_vals[1][1])),
              np.sqrt(abs(cov_vals[2][2])),
              np.sqrt(abs(cov_vals[1][2]))
              ]
    return fit_vals, errors

def plot_trimmed_data(xs,ys, trimxs, trimys, wavelength_min,wavelength_max):

    plt.plot(xs, ys, color='k', label='Underlying Data')
    plt.plot(trimxs, trimys, label='Trimmed Data')
    
    
    plt.xlim([wavelength_min,wavelength_max])
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Intensity [counts]')
    plt.legend()
    plt.show()

def plot_fitted_data(xs, ys, fit_vals, wavelength_min,wavelength_max):
    N = 100
    x_fits = np.linspace(wavelength_min, wavelength_max, N)
    y_fits = gaussian(fit_vals, x_fits)
    plt.plot(xs, ys, color='k', label='Underlying Data')
    plt.plot(x_fits, y_fits, label='Fitted Data')
    plt.xlim([wavelength_min,wavelength_max])
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Intensity [counts]')
    plt.legend()
    plt.show()

def find_Te(B, B_err, C, C_err, BC_cov, index,
            is_in_Kelvin=False, is_new_method=False):
    # convert B,C from nm to m
    B*=10**(-9)
    B_err*=10**(-9)
    C*=10**(-9)
    C_err*=10**(-9)
    BC_cov*=10**(-18)
    
    lambda_i = None
    if is_new_method:
        lambda_i = B
    else:
        lambda_i = 694.3*10**(-9) #m
    theta = angle[index]#dependent on slice
    m_e = 9.11*10**(-31) #kg
    c = 3.00*10**(8) #m/s
    kB = 1.38*10**(-23) #m^2 kg s^(-2) K^(-1)
    const = m_e*c*c/(lambda_i*lambda_i*4*kB)
    #print(const)
    Te = (const)*pow(C/np.sin(theta/2),2)
    if not is_in_Kelvin:
        K_to_eV_const = 11604.45
        Te/=K_to_eV_const
    # find error
    Te_err = None
    if is_new_method:
        Te_err = 2*Te*np.sqrt(pow(B_err/B,2) +
                              pow(C_err/C,2) - 
                              2*(BC_cov/(B*C)))
    else:
        Te_err = 2*Te*C_err/C
    return Te, Te_err

def perform_analysis(index, is_plot=True):
    # get data slice
    xs, ys = get_slice(data_array, index)
    
    # provide limits
    wavelength_min = 550
    wavelength_max = 950
    
    # trim data
    trimxs, trimys = remove_lines(xs, ys, wavelength_min, wavelength_max, is_plot_lines=is_plot)
    if is_plot:
        plot_trimmed_data(xs,ys, trimxs, trimys, wavelength_min,wavelength_max)
    
    # fit data
    fit_vals, fit_errors = perform_fit(trimxs, trimys)
    if is_plot:
        plot_fitted_data(xs, ys, fit_vals, wavelength_min,wavelength_max)
    B = [fit_vals[1], fit_errors[1]]
    C = [fit_vals[2], fit_errors[2]]
    BC_cov = fit_errors[3]
    
    return B, C, BC_cov

def plot_Te_vs_radius(xss, yss, y_errss, fmts, labels, num):
    for i in range(num):
        plt.errorbar(xss[i], yss[i], y_errss[i],
                     fmt=fmts[i],
                     label=r"{}".format(labels[i])) 
    #plt.xlim([0.3,1.4])
    plt.ylim([0,600])
    plt.xlabel('Radius [m]')
    plt.ylabel(r'$T_e$ [eV]')
    plt.legend()
    plt.show()

def plot_Bs_vs_radius(xs, ys, y_errs):
    plt.errorbar(xs, ys, y_errs, fmt='.k', label='With Propagated Error') 
    #plt.xlim([0.3,1.4])
    plt.ylim([680,700])
    plt.xlabel('Radius [m]')
    plt.ylabel(r'$B$-fit Parameter [nm]')
    laser_wavelength = 694.3  # in nm
    plt.axhline(y=laser_wavelength, label='Ruby Laser Wavelength')
    plt.legend()
    plt.show()

# perform_analysis(100)
Tes = []
Te_errs = []
Tes_new = []
Te_errs_new = []
Tes_new2 = []
Te_errs_new2 = []
Bs = []
B_errs = []
for i in range(len(data_array)):
    B, C, BC_cov = perform_analysis(i, is_plot=False)
    Te, Te_err = find_Te(B[0], B[1], C[0],C[1],BC_cov,i)
    Tes.append(Te)
    Te_errs.append(Te_err)
    
    Te_new, Te_err_new = find_Te(B[0], B[1], C[0],C[1],BC_cov,i, is_new_method=True)
    Tes_new.append(Te_new)
    Te_errs_new.append(Te_err_new)
    
    lambda_i = 694.3 #nm
    B_error = abs(B[0] - lambda_i)
    Te_new2, Te_err_new2 = find_Te(B[0], B_error, C[0],C[1],BC_cov,i, is_new_method=True)
    Tes_new2.append(Te_new2)
    Te_errs_new2.append(Te_err_new2)
    
    Bs.append(B[0])
    B_errs.append(B[1])
    
#plot_Bs_vs_radius(radius, Bs, B_errs)
plot_Te_vs_radius([radius,radius,radius],
                  [Tes, Tes_new, Tes_new2],
                  [Te_errs, Te_errs_new, Te_errs_new2],
                  ['.k', '.b', '.r'],
                  [r"Fixed $\lambda_i$", r"$\lambda_i = B$",
                   r"$\lambda_i = B$, $\sigma_B = |B-\lambda_i|$"],
                  3)
