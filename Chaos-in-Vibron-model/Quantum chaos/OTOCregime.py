import matplotlib.pyplot as plt
import numpy as np
import palettable
import math
import PlotVMfinal
import cycler
import statistics
from matplotlib.gridspec import GridSpec
import matplotlib.pylab as pl
from math import floor, log10
plt.rcParams.update({'font.size': 12})

plt.rcParams['text.usetex'] = True

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits

    return r"${0:.{2}f}\times 10^{{{1:d}}}$".format(coeff, exponent, precision)

def ReadData(path):
    Data = []
    file = open(path,"r")
    lines = file.readlines()

    DataAllLines = []
    for line in lines:

    #print(line)
        line = line.replace("[","")   
        line = line.replace("]","")
        line = line.replace(",","")
        line = line.replace(";","")
        line = line.replace("Any","")
            
        elements = line.split()
        Data = [float(element) for element in elements]
        #print(Data)
        DataAllLines.append(Data)
    
    return DataAllLines

#fig = plt.figure(figsize=(5, 5))
f,(ax0,ax1) = plt.subplots(1,2,sharey=True, facecolor='w',figsize=(5, 2))
#gs = GridSpec(nrows=1, ncols=2)
# For Sine Function

N = 15
NumLevels = int((N+1)*(N+2)/2)
xi,epsilon,OTOCname = 0.1,0.4,"[Dx(t),Dx]SHORT"
length = 1000

firstlevel = 1
secondlevel = 50

pars = " N=%1.0i,chi=%1.2f,eps=%1.2f len = %1.0i" % (N,xi,epsilon,length)

OTOC = 'OTOC'
pth = 'OTOCdata/' + OTOCname + pars
    
Spectrum = ReadData(pth+"Spectrum.txt")[0]
Values = np.array(ReadData(pth+"Values.txt")[0]).reshape((length,NumLevels))
Times = ReadData(pth+"Times.txt")[0]

#
#   long
#
xi,epsilon,OTOCname = 0.1,0.4,"[Dx(t),Dx]LONG"
length = 1000

pars = " N=%1.0i,chi=%1.2f,eps=%1.2f len = %1.0i" % (N,xi,epsilon,length)

OTOC = 'OTOC'
pth = 'OTOCdata/' + OTOCname + pars
    
SpectrumL = ReadData(pth+"Spectrum.txt")[0]
ValuesL = np.array(ReadData(pth+"Values.txt")[0]).reshape((length,NumLevels))
TimesL = ReadData(pth+"Times.txt")[0]


#ax0 = fig.add_subplot(gs[0, 0])
#ax0.plot(Times,Values[:,1], label = Spectrum[1]/N)
name1 =r'$E = %1.2f, \bar{C} = $ '% (Spectrum[firstlevel]/N) 
name1 = name1 + sci_notation(statistics.mean(ValuesL[:,1]),1)

name2 =r'$E = %1.2f, \bar{C} = $ '% (Spectrum[secondlevel]/N) 
name2 = name2 + sci_notation(statistics.mean(ValuesL[:,secondlevel]),1)

name1 = r'$E = %1.2f$'% (Spectrum[firstlevel]/N) 
name2 = r'$E = %1.2f$'% (Spectrum[secondlevel]/N) 

ax0.plot(Times[:],Values[:,firstlevel]/statistics.mean(ValuesL[:,firstlevel]), alpha = 0.8)
ax0.plot(Times[:],Values[:,secondlevel]/statistics.mean(ValuesL[:,secondlevel]), alpha = 0.8,lw = 3)



Value =  statistics.mean(ValuesL[:,secondlevel]) - np.sqrt(statistics.variance(ValuesL[:,secondlevel]))
#ax0.axhline(y= Value, color='r', linestyle='-')

listFor = ValuesL[:,secondlevel]
print(listFor)
res = next(x for x, val in enumerate(list(listFor))
                                  if val > Value)
print(res)
print(Value)
ax0.axvline(x = 14, color = 'black', linestyle = 'dashed',alpha = 0.5)

Value =  (statistics.mean(ValuesL[:,secondlevel]) - np.sqrt(statistics.variance(ValuesL[:,secondlevel])))/statistics.mean(ValuesL[:,secondlevel])
#ax0.axhline(y= Value, color='g', linestyle='-')
#ax0.set_yscale('log')

#ax1 = fig.add_subplot(gs[0, 1])
ax1.plot(TimesL,ValuesL[:,firstlevel]/statistics.mean(ValuesL[:,firstlevel]),alpha = 0.8,label = name1)
ax1.plot(TimesL,ValuesL[:,secondlevel]/statistics.mean(ValuesL[:,secondlevel]), alpha = 0.8,lw = 3,label =name2)



#ax1.plot(TimesL,ValuesL[:,650]/statistics.mean(ValuesL[:,650]), alpha = 0.8,lw = 3,label =name2)

#ax1.set_yscale('log')
#ax1.legend()

ax0.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax0.yaxis.tick_left()
ax1.tick_params(labelright='off')
ax1.yaxis.tick_right()
#ax1.axes.yaxis.set_ticklabels([])

#ax1.set_xticks([ 0.4*1e9,0.8*1e9], [r'$0.4 \times 10^9$',r'$0.8 \times 10^9$'])
#ax1.set_xlabel(r'$t$')

#ax0.text(156,-0.5, r"$t$", ha='center', va='center', size=12)
ax1.text(100,-0.5, r"$t$", ha='center', va='center', size=14)

ax1.xaxis.set_label_coords(0.95, -0.025)

#ax0.set_ylabel(r'$\frac{C(t)}{\bar{C}}$', rotation=0,size=14)
ax0.text(-25,1.5,r'$\frac{C(t)}{\bar{C}}$', ha='center', va='center', size=14)


d = .015 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax0.transAxes, color='k', clip_on=False)
ax0.plot((1-d,1+d), (-d,+d), **kwargs)
ax0.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
ax1.plot((-d,+d), (1-d,1+d), **kwargs)
ax1.plot((-d,+d), (-d,+d), **kwargs)



plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)

#plt.show()
plt.savefig("OTOC.png",
            dpi = 1000,
            bbox_inches ="tight")
plt.savefig("OTOC.pdf",
            dpi = 1000,
            bbox_inches ="tight")