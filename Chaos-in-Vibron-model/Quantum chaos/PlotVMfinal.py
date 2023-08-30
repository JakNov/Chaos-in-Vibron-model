import matplotlib.pyplot as plt
import numpy as np
import math
import csv
import palettable
import statistics
plt.rcParams['text.usetex'] = True

def ReadData(path):
    Data = []
    file = open(path,"r")
    line = file.readline()

    #print(line)
    line = line.replace("[","")   
    line = line.replace("]","")
    line = line.replace(",","")
    line = line.replace(";","")
    line = line.replace("Any","")
            
    elements = line.split()
    Data = [float(element) for element in elements]
    
    return np.array(Data)

def SmoothFull(data,n):
    length = len(data)
    #print(length)
    SmoothData = [sum(data[0:i+n])/(i+n) for i in range(0,n)]
    print(len(SmoothData))

    SmoothData = np.append(SmoothData,[sum(data[i-n:i+n])/(2*n + 1) for i in range(n,length-n)])
    print(len(SmoothData))
    
    SmoothData = np.append(SmoothData,[sum(data[length-2*n+i:length])/(2*n - i+1) for i in range(0,n)])
    print(len(SmoothData))
    return SmoothData


def SmoothBlock(xdata,ydata,n):

    newXData = []
    newYData = []
    varData = []

    print(range(math.floor(len(xdata)/n)))
    for i in range(int(math.floor(len(xdata)/n))):
        i1 = i*n +1
        i2 = (i+1)*n +1

        np.append(newYData,sum(ydata[i1:i2])/n)
        np.append(newXData,sum(xdata[i1:i2])/n)
        np.append(varData,np.sqrt(statistics.variance((ydata[i1:i2]))))
    

    return newXData, newYData, varData


def ResizeData(data):
    data = np.array(data)
    data = data - min(data)
    return data/max(data)


def ResizeDataSm(data,Sm):
    datas = SmoothFull(data,Sm)
    data = np.array(data)
    data = data - min(datas)
    return data/max(data - min(datas))

def PlotVM(N,xi,epsilon,OTOCname):
    dim = (N+1)*(N+2)/2

    #pars = " N=%1.0i,ξ=%1.2f,ϵ=%1.2f" % (N,xi,epsilon)
    pars = " N=%1.0i,chi=%1.2f,eps=%1.2f" % (N,xi,epsilon)

    OTOC = 'OTOC'
    pth = 'vysledky/' + OTOC + OTOCname + pars
    
    Spectrum = ReadData(pth+"Spectrum.txt")
    Mean = ReadData(pth+"Mean.txt")
    Var = ReadData(pth+"Var.txt")

    sm = 25
    nameFreg = "vysledky/freg_%1.1f_%1.1f.txt" %(xi,epsilon)
    nameLya = "vysledky/lyapunov_%1.1f_%1.1f.txt"%(xi,epsilon)


    DataFreg = np.loadtxt(nameFreg, delimiter = "\t")
    DataLya = np.loadtxt(nameLya, delimiter = "\t")

    fig, ax = plt.subplots(figsize=(5,3))

    im = ax.scatter(Spectrum/N,Var/Mean,s = 10,alpha = 0.3,c = 'black')
    im.set_edgecolor("none")
    ax.grid()
    ax.set_xlabel(r"$E$", color="k", size=10)
    ax.set_ylabel(r"$\nu$", color="k", size=10)


    ax.set_prop_cycle('color', palettable.cartocolors.qualitative.Safe_10.mpl_colors)
    plt.plot(DataLya[:,0],DataLya[:,1],lw=3)#,c = 'green')
    plt.plot(DataFreg[:,0],DataFreg[:,1]/2,lw=3)#,c = 'royalblue',lw = 3)
    plt.plot(Spectrum/N,SmoothFull(Var/Mean,sm),lw=3)#,c='gold',lw = 3)

    

    plt.show()

    return ax

#PlotVM(100,0.8,0.4,'[W2(t),n]')
#PlotVM(100,0.8,0.4,'[Dx(t),Dx]')
#PlotVM(100,0.4,0.4,'[Dx(t),Dx]')

#a1,a2,a3=SmoothBlock([1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8],4)
#print(a1)