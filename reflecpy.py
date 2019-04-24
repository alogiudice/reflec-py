import numpy as np
#from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy import stats
import statsmodels.api as sm
from math import pow

auto_optim = True

class Reflect():
    
    def __init__(self, x, counts):
        self.x = x
        self.counts = counts
        
    def thetacrit(self, maxindex_crop = 1.6):
        # Esta función devuelve el valor del ángulo crítico del sistema, teniendo 
        # en cuenta a éste calculado como el ángulo en el que I = Imax / 2.
        max_value = max(self.counts)
        max_index = self.counts.index(max_value)
        thetacrit_value = max_value / 2
        # Buscamos el valor más cercano al del I/2 y le asignamos el valor de 
        # theta crítico.
        maxindex_scan = int(max_index * maxindex_crop)
        # Croppeamos el array hasta este valor.
        scan = self.counts[max_index + 1:maxindex_scan]
        scan = [abs(value - thetacrit_value) for value in scan]
        mini = np.argmin(scan)
        thetacritic = x[max_index + 1 + mini]
        
        return thetacritic
    
    def smoothcounts(self, fraclowess, fig_smooth, ax, tester, x_cutoffindex = 1.705):
        # Esta función grafica los datos sin smoothear, los ya smootheados por 
        # lowess, y a los picos encontrados. Devuelve el array peaks_index, en el 
        # que están especificadas las posiciones de los picos en el rango 
        # [max_index:x_cutoff]
        
        # La idea es ver ahora cómo obtener los máx y mín de la curva. Es un proceso
        # que parece complicado, porque vamos a tener que smoothear la función.
        # Probé con savgol y no dio muy bien para el tipo de señal que tenemos.
        # También habría que poner hasta qué pico queremos smoothear.
        
        
        max_value = max(self.counts)
        max_index = self.counts.index(max_value)
        x_cutoff = self.x.index(x_cutoffindex)
    
        # Croppeamos los arrays para realizar los ajustes en este rango.
        
        xx = self.x[max_index:x_cutoff]
        yy = self.counts[max_index:x_cutoff]
        
        
    
        lowess = sm.nonparametric.lowess(yy, xx, frac=fraclowess)
        
        peaks_index = argrelextrema(lowess[:,1], np.greater)
        peaks_index_zero = [max_index + y for y in peaks_index[0]]
        print("Found %d local maxima." % len(peaks_index[0]) )
        
        peak_list = []
        xpeaks = []
        for peak in peaks_index[0]:
            peak_list.append(lowess[peak, 1])
            xpeaks.append(xx[peak])
            
        if tester is False:
            line1, = ax.plot(xx, yy, color='orange', label='data')
            line2, = ax.plot(xx, lowess[:,1], color='red', label="Smooth")
            #line3, = ax.plot(xpeaks, 2, color='grey', marker='s', label= 'Found maxima')
            ax.set_ylabel('counts')
            ax.set_xlabel('2theta[deg]')
            ax.set_yscale('log')
            ax.legend(loc='lower left')
        else:
            line2.set_xdata(xx)
            line2.set_ydata(lowess[:, 1])
            #line3.set_xdata()
            #line3.set_ydata(peak_list)
            fig_smooth.canvas.draw()
            
            
        tester = True
        return peaks_index_zero, tester


    def get_slope(self, peaks_index_zero, n):
        # Esta función obtiene los parámetros del ajuste lineal del sin**2 de los 
        # picos con los valores de n. Toma un array c con las posiciones en 2theta
        # (en deg) de los picos y el valor n inicial para generar la sucesión de 
        # de etiquetas.
        
        angpeaks = [0] * len(peaks_index_zero)
        angpeaks = [np.sin(pow(self.x[peak] * (np.pi / 360), 2)) for peak in peaks_index_zero]
        npeaks = np.arange(n, n + len(peaks_index_zero), 1)
        npeaks = [pow(nn * lambdax / 2, 2) for nn in npeaks]
        
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(npeaks, angpeaks)
        thetacrit_int = self.thetacrit()
        thetacrit_exp = pow(intercept, 1/2) * (360 / np.pi)
        thetacrit_diff = ( (thetacrit_exp - thetacrit_int) / thetacrit_int ) * 100
        
        fittedline = [m * slope + intercept for m in npeaks]
        
        line3, = ax2.plot(npeaks, angpeaks, 'o')
        line4, = ax2.plot(npeaks, fittedline, color='green')
    
        print('Fitting results for data:')
        print("slope: %e,\nIntercept: %e,\nR-squared: %f" % (slope, intercept, r_value))
        print("2theta_crit value obtained is %f . Value obtained from counts was %f" % (thetacrit_exp, thetacrit_int))
        print("Percentage difference between both critical angles is %f%%" % thetacrit_diff)
        return slope, intercept, r_value

    
    def prompt_n(n):
        # Esta función la definimos para ir cambiando el valor de n según lo que 
        # decida el usuario.
        while 1 != 0:
            prompt = input("Change n value? (Current n is %f):  " % n)
            if prompt == 'y':
                n_new = input("Specify new n value:  ")
                n = float(n_new)
                break
            elif prompt == 'n':
                break
            else:
                print('Specify y/n.')
        return n


##############################################################################
 ############################################################################
  ########################### PROGRAM START ###############################
#plt.use("TkAgg")
plt.ion()
tester = False
fileopen = "/home/agostina/git/reflec-py/XRR_Real1396.dat"
file = open(fileopen, "r")
x = []
counts = []
n = 1
lambdax = 1.5406 # La longitud de onda de los rayos X usados (Cu)

fig_smooth = plt.figure()
ax = fig_smooth.add_subplot(111)

#fig = plt.figure()


with open(fileopen) as f:
     data = f.readlines()
     for line in data:
         values = line.split()
         x.append(float(line.split(' ')[0]))
         counts.append(float(line.split(' ')[1]))
         
file.close()
