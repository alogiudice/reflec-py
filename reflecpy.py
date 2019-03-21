import numpy as np
#from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy import stats
import statsmodels.api as sm
from math import pow


def thetacrit(counts, x):
    # Esta función devuelve el valor del ángulo crítico del sistema, teniendo 
    # en cuenta a éste calculado como el ángulo en el que I = Imax / 2.
    max_value = max(counts)
    max_index = counts.index(max_value)
    print('Max counts value is: %f' % max_value)
    thetacrit_value = max_value / 2
    # Buscamos el valor más cercano al del I/2 y le asignamos el valor de theta
    # crítico.
    # scanemos hasta 2 veces el ángulo en el que I = max.
    maxindex_scan = int(max_index * 1.6)
    # Croppeamos el array hasta este valor.
    scan = counts[max_index + 1:maxindex_scan]
    
    for i in range(0, len(scan)):
        scan[i] = abs(scan[i] - thetacrit_value)    
        
    mini = np.argmin(scan)
    thetacritic = x[max_index + 1 + mini]
    return thetacritic

def get_slope(c, n, fig):
    # Esta función obtiene los parámetros del ajuste lineal del sin**2 de los 
    # picos con los valores de n. Toma un array c con las posiciones en 2theta
    # (en deg) de los picos y el valor n inicial para generar la sucesión de 
    # de etiquetas.
    angpeaks = np.zeros(len(c))
    npeaks = np.arange(n, n + len(c), 1)

    for i in range(0, len(c)):
        angpeaks[i] = np.sin(pow(c[i] * (np.pi / 360), 2))
        npeaks[i] = pow(npeaks[i] * lambdax / 2, 2)

    slope, intercept, r_value, p_value, std_err = stats.linregress(npeaks, angpeaks)
    thetacrit_exp = pow(intercept, 1/2) * (360 / np.pi)
    ax[1].clear()
    ax[1].grid()
    ax[1].plot(npeaks, angpeaks, 'o')
    ax[1].plot(npeaks, npeaks * slope + intercept, color='green')
    
    print('Fitting results for data:')
    print("slope: %e intercept: %e, R-squared: %f" % (slope, intercept, r_value))
    print("2theta_crit value obtained is %f . Value obtained from counts was %f" % (thetacrit_exp, thetacrit_int))
    prompt1 = input("Change max to min? (y/n): ")
    prompt2 = input("Increment n by one? (Current n is %f): " % n)
    return slope, intercept, r_value, prompt1, prompt2


def smoothcounts(counts, x, fraclowess, x_cutoff, fig):
    # Esta función grafica los datos sin smoothear, los ya smootheados por 
    # lowess, y a los picos encontrados. Devuelve el array peaks_index, en el 
    # que están especificadas las posiciones de los picos en el rango 
    # [max_index:x_cutoff]
    
    #La idea es ver ahora cómo obtener los máx y mín de la curva. Es un proceso
    # que parece complicado, porque vamos a tener que smoothear la función.
    # Probé con savgol y no dio muy bien para el tipo de señal que tenemos.
    #También habría que poner hasta qué pico queremos smoothear.
    
    max_value = max(counts)
    max_index = counts.index(max_value)

    xx = x[max_index:x_cutoff]
    yy = counts[max_index:x_cutoff]

    
    lowess = sm.nonparametric.lowess(yy, xx, frac=fraclowess)

    # Ojo: Estoy casi segura que la frac del lowess deberá cambiarse dependiendo
    # del caso. Agregar un input para esto.
    ax[0].clear()
    ax[0].grid()
    ax[0].scatter(xx, yy, color='orange')
    ax[0].plot(xx, lowess[:, 1], color='red')
    ax[0].set_yscale('log')

    peaks_index = argrelextrema(lowess[:,1], np.greater)
    print("Found %d local maxima." % len(peaks_index[0]) )
    c = np.zeros(len(peaks_index[0]))

    for i in range(0, len(peaks_index[0])):
        k = peaks_index[0][i]
        c[i] = xx[k]
        z = lowess[k,1]
        ax[0].plot(c[i], z, color='grey', marker='s')
        plt.draw()
    
    ax[0].legend(('Smooth', 'Data', 'Found maxima'), loc='upper right')
    ques = input("Counts smoothing is OK? (y/n) (If not, specify a new LOWESS fraction, current is %f)" % fraclowess)
    return peaks_index, ques, c


##############################################################################
 ############################################################################

fileopen = "XRR_Real1396.dat"
file = open(fileopen, "r")
x = []
counts = []
lambdax = 1.5406 # La longitud de onda de los rayos X usados (Cu)
fig, ax = plt.subplots(2)


with open(fileopen) as f:
     data = f.readlines()
     for line in data:
         values = line.split()
         x.append(float(line.split(' ')[0]))
         counts.append(float(line.split(' ')[1]))

file.close()

thetacrit_int = thetacrit(counts, x)
x_cutoff = x.index(1.705)
plt.draw()

fraclowess = 0.05
peaks_index, ques, c = smoothcounts(counts, x, fraclowess, x_cutoff, fig)
promptloop = False

while promptloop == False:
    if ques == "y":
        break
    elif ques == "n":
        fraclowess = input("Specify a new LOWESS fraction. Current one is %f. : " % fraclowess)
        peaks_index, ques, c = smoothcounts(counts, x, float(fraclowess), x_cutoff, fig)

#Ahora ya tenemos la ubicación de los picos y los valores de éstos. Ya 
#podríamos encontrar el valor del ancho de la película y del ángulo crítico.


#Calculamos el valor del theta crítico y lo comparamos con el obtenido anterior
#mente.

n = 1
slope, intercept, r_value, prompt1, prompt2 = get_slope(c, n, fig)

while promptloop == False:
    if prompt1 == 'y' and prompt2 == 'n' :
        n += 1/2
        get_slope(c,n, fig)
    elif prompt1 == 'y' and prompt2 == 'y':
        n += 3/2
        get_slope(c,n, fig)
    elif prompt1 == 'n' and prompt2 == 'y':
        n += 1
        get_slope(c,n, fig)
    elif prompt1 == 'n' and prompt2 == 'n':
        break
    break
# Calculamos el espesor de la capa especificada.
width = pow(1 / slope, 1/2)
width /= 10

print("Layer thickness has been calculated at %f nm." % width)
plt.show()
# plt.scatter(npeaks, angpeaks)
#plt.ylim(min(angpeaks), max(angpeaks))






