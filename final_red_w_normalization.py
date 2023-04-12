# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:53:04 2023

@author: weisz
"""
import pylab as pl
from matplotlib import pyplot as plt 
import numpy as np 
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import math
from mpmath import linspace
from pylab import figure, cm
from matplotlib.colors import LogNorm
import pickle
#%%%
def loading(plotting = False):
    abc = np.loadtxt('E:/nano/2 peak fitting/anyagmegmondo/bulk_ABC_red.txt' )
    abab = np.loadtxt('E:/nano/2 peak fitting/anyagmegmondo/red_bulk_AB_sample6.txt')
    
    if plotting == True:
        fig, ax1 = plt.subplots(figsize=(6.5, 4))
        
        ax1.plot(abab[:, 0], abab[:, 1], 'b')
        ax1.plot(abc[:, 0], abc[:, 1], 'r')
      
        ax1.set_xlim([2400, 2800])
        ax1.set_xlabel('Red Raman shift [$\mathrm{cm^{-1}}$]', size=15)
        ax1.set_ylabel('intensity [a.u.]', size=15)
        ax1.tick_params(axis='x', labelsize=15)
        ax1.tick_params(axis='y', labelsize=15)
        plt.grid()
    
        plt.show()

    zoomabx=[]; zoomaby =[]; zoomabcx=[]; zoomabcy=[];
    
    for i in range(len(abab)):
        if abab[i][0] > 2500 and abab[i][0] < 3000:
            zoomabx.append(abab[i][0])
            zoomaby.append(abab[i][1])
                                 
    zoomabx = np.array(zoomabx);zoomaby = np.array(zoomaby)
    
    
    
    for i in range(len(abc)):
        if abc[i][0] > 2500 and abc[i][0] < 3000:
            zoomabcx.append(abc[i][0])
            zoomabcy.append(abc[i][1])
            
    zoomabcx = np.array(zoomabcx);zoomabcy = np.array(zoomabcy)
    #zoomabx,zoomaby,zoomabcx,zoomabcy=loading()
    #zoomabcy=zoomabcy-sum(zoomabcy)/len(zoomabcy);zoomaby=zoomaby-sum(zoomaby)/len(zoomaby);
    zoomabcy=zoomabcy-min(zoomabcy);zoomaby=zoomaby-min(zoomaby);
    na=max(zoomaby);naby=zoomaby/na
    nb=max(zoomabcy);nabcy=zoomabcy/nb
    nabcy=np.array(nabcy);naby=np.array(naby)
    mb=(zoomabx[np.unravel_index(np.argmax(zoomaby, axis=None),zoomaby.shape)[0]])
    mc=(zoomabcx[np.unravel_index(np.argmax(zoomabcy, axis=None),zoomabcy.shape)[0]])
    delta=mc-mb
    
    #avgabc = sum(zoomabcy) / len(zoomabcy)
    return zoomabx, zoomaby , zoomabcx, zoomabcy,na,nb,naby,nabcy,delta

def lorentzian2( x,xb0 ,yb0,w2,q,a1, b1, g1, a2, b2, g2):
    return ((yb0+((a1*q) *w2* g1**2) /  ( w2* g1**2 +( x - b1 + xb0 )**2))+ 
            (yb0+((a2*q) *w2* g2**2) / ( w2*g2**2 +  ( x - b2 + xb0 )**2)))

def abfugg(x, x0ab,y0ab,w2):
        return lorentzian2(x, x0ab, y0ab,w2,*poptab[3::])
    
def l2fit():
    yb0 = 0;xb0 = 0;w2=1;a1 = 1;b1 = 2678;g1 = 15;a2 = 0.2;b2 = 2635;g2 = 20;q=0.1
    poptab, pcovab = curve_fit(lorentzian2, zoomabx,naby, p0=[xb0,yb0,w2,q, a1, b1, g1, a2, b2, g2])
    #perrab = np.sqrt(np.diag(pcovab)[::])
    abfit = lorentzian2(zoomabx, *poptab)
    return abfit,poptab

def lorentzian4(x, x0,y0,q,  a1, b1, g1, a2, b2, g2, a3, b3, g3, a4, b4, g4):
    return ((y0 + ((a1 * q) * g1 ** 2) / (g1 ** 2 +  (x - b1 + x0) ** 2)) + 
                (y0 + ((a2 * q) * g2 ** 2) / (g2 ** 2 + (x - b2 + x0) ** 2)) +
                +(y0 + ((a3 * q) * g3 ** 2) / (g3 ** 2 + (x - b3 + x0) ** 2)) + 
                (y0 + ((a4 * q) * g4 ** 2) / (g4 ** 2 + (x - b4 + x0) ** 2)))

def abcfugg(x, x0,y0):
        return lorentzian4(x, x0,y0,*poptabc[2::])
       
def abcfugg2(x, x0):
        return lorentzian4(x, x0, *poptabc[1::])
def l4fit():
     # fit
    x0 = 0; y0 = 0;q=0.1
    a1 = .96; b1 = 2670; g1 = 17;
    a2 = .4;b2 = 2620;g2 = 35;
    a3 = .2 ;b3 = 2590;g3 = 25;
    a4 = .89;b4 = 2650;g4 = 18
    
    poptabc, pcovabc = curve_fit(lorentzian4, zoomabcx, nabcy,
                                 p0=[x0,y0, q,a1, b1, g1, a2, b2, g2, a3, b3, g3, a4, b4, g4])
    
    #perrabc =  np.sqrt(np.diag(pcovabc)[::])
    abcfit = lorentzian4(zoomabcx, *poptabc)
    
    return abcfit,poptabc

def readmap(map_name):
    mp = np.loadtxt(map_name)
    ini=[];X=[]
    for i in range(len(mp)):
        ini.append(mp[i][0])
        if i > 1150 and i < 1550:
            X.append(mp[i][0])
    mp=np.delete(mp, 0, axis=1)
    nagymatrix=[];matrix=[]
    for u in range(0,len(mp[0])):
        pixel=[];szurt=[]
        for i in range(0,len(mp)):
            pixel.append(mp[i][u])
            if i > 1150 and i < 1550:
                szurt.append(mp[i][u])
        nagymatrix.append(pixel)
        matrix.append(szurt)
    return nagymatrix, matrix,ini,X

def finder(nagymatrix,g_peak): 
    "Find graphite with G peak (g_peak: avg size of g_peak) "
    array_nagy=np.array(nagymatrix)
    non_graphene = array_nagy < g_peak 
    return non_graphene

def find_peaks(Ab_map_index,Abc_map_index,xb0= -9.65382706e+02 ,yb0=1 ,w2=1,x0=-10,y0=0):
    "find peaks in x scale "
    "xb0,yb0- params for ab"
    "xb,yb- params for abc"
    
    ya=(matrix[Ab_map_index])#/max(matrix[Ab_map_index])  #-min(matrix[Ab_map_index])
    pl.plot(x,ya,label='ab spect')
    def abxslide(x,ya,xb0,yb0):
        ba, bba =curve_fit(abfugg,x,y,p0=[xb0,yb0,w2],bounds=((-3e+5,-1000,0),(3e+5,1000,20)))
        pl.plot(x,abfugg(x,*ba),label='ab fit')
        AB =abfugg(x,*ba)
        print(ba)
        return AB, ba
    
    y=(matrix[Abc_map_index])#/max(matrix[Abc_map_index]) #-min(matrix[Abc_map_index])
    pl.plot(x,y,label='abc spect')
    def abcxslide(x,y,x0,y0):
        cba, ccba =curve_fit(abcfugg,x,y,p0=[x0,y0],bounds=((-3e+05,-1000),(3e+05,1000)))
        pl.plot(x,abcfugg(x,*cba),label='abc fit')
        plt.legend()
        print(cba)
        ABC = abcfugg(x,*cba)
        return ABC, cba
    
    plt.legend()
    AB,ba= abxslide(x,y0,xb0,yb0)
    ABC,cba= abcxslide(x,y,x0,y0)
    
    return ABC,cba,AB,ba

def mapping(x,y,egylista,errAB,errABC,find_non_graph):
    def tog(x,h1,t1,t2,zoomsize):
        return((h1*(ABC+t1) + (1-h1)*(AB+t2))*zoomsize)
    
    popt,pcov = curve_fit(tog,x,y,bounds=((0,-1000,-1000,0),
                                      (1,1000,1000,10)))
    errAB= np.sqrt((AB - tog(x,*popt))**2)
    errABC = np.sqrt((ABC-tog(x,*popt))**2)
    
    errbar = np.sum(np.abs(errAB[::]))/len(errAB[::])+np.sum(np.abs(errABC[::]))/len(errABC[::])
    
    if find_non_graph == True and (np.sum(tog(x,*popt)) /(len(tog(x,*popt)) / errbar)) > 0.008:
        popt=[0,0,0,0]
    egylista.append(popt)
    
    return(egylista,errAB,errABC)  

def mapcalc(picklename):
    for i in range (len(matrix)):
        errAB=[];errABC=[]
        egylista=[];y= matrix[i]; 
        find_non_graph=np.all(non_graphene[i][400:900]==True)
        try:
            egylista,errAB,errABC= mapping(x,y,egylista,errAB,errABC,find_non_graph)
        except(RuntimeError):
            egylista=[-1,0]
            errAB=[-10,0]
            errABC=[-10,0]
            
        kiiertekeltmap.append(egylista)
        ABerrmap.append(errAB)
        ABCerrmap.append(errABC)

        
    szazalek=[]
    a = np.array([1, 2, 3])
    for i in range(len(kiiertekeltmap)):
        if type(kiiertekeltmap[i][0]) == type(a):
            #szazalek.append((kiiertekeltmap[i][0][0]/((kiiertekeltmap[i][0][0])+(kiiertekeltmap[i][0][1])))*100)
            szazalek.append(kiiertekeltmap[i][0][0]*100)
        else:
            szazalek.append(-1)
    with open(picklename +'.pickle', 'wb') as handle:
        pickle.dump(kiiertekeltmap, handle)
    with open(picklename+'AB_error' +'.pickle', 'wb') as handle:
        pickle.dump(ABerrmap, handle)
    with open(picklename+'ABC_error' +'.pickle', 'wb') as handle:
        pickle.dump(ABCerrmap, handle)    
    return szazalek,ABerrmap,ABCerrmap

def picture(width,szazalek,ABerrmap,ABCerrmap,AB_pos,ABC_pos,er_min,er_max,log=False):
    def pic_map(datas):
        ekep=[];o=0;esor=[]
        while o < len(datas):
            esor.append(datas[o])
            if o !=0 :
                if width != 0 : 
                    if (o+1) % width == 0:
                        ekep.append(esor)
                        esor=[]
                elif (o+1) % round(np.sqrt(len(datas)))== 0:
                    ekep.append(esor)
                    esor=[]
            o += 1
        return np.array(ekep)     
        
    def heatmap2d(arr: np.ndarray):
        pl.imshow(arr, cmap='gist_ncar',interpolation='none') #'gist_ncar'
        pl.colorbar()
        pl.show()
        
    kep=pic_map(szazalek)
    heatmap2d(kep)
    
    ablocerr=[]
    for i in range(len(ABerrmap)):
        ablocerr.append(np.sum(np.abs(ABerrmap[i][::]))/len(ABerrmap[i])+np.sum(np.abs(ABerrmap[i][::]))/len(ABerrmap[i]))
    abclocerr=[]
    for i in range(len(ABerrmap)):
        abclocerr.append(np.sum(np.abs(ABCerrmap[i][::]))/len(ABCerrmap[i])+np.sum(np.abs(ABCerrmap[i][::]))/len(ABCerrmap[i])) 
    #err_min=((ablocerr[AB_pos]+abclocerr[ABC_pos])/2+(ablocerr[AB_pos]+abclocerr[ABC_pos])/5),
    #err_max=((ablocerr[AB_pos]+abclocerr[ABC_pos])/2-(ablocerr[AB_pos]+abclocerr[ABC_pos])/5),
    def heatmap2d(arr: np.ndarray,tit):
        pl.imshow(arr, cmap='inferno',interpolation='none', vmin=er_min, vmax=er_max) #'gist_ncar'
        pl.title(tit)
        pl.colorbar()
        pl.show()
    #for i in ABCerrmap:
     #   ABC_locerr.append(np.sum(np.abs(i[::]))/len(i))
   
    abcerror=pic_map(abclocerr);aberror=pic_map(ablocerr)
    heatmap2d(abcerror,'abc error');heatmap2d(aberror,'ab error')
    
    if log == True:
        def logkep(ekep,tit):
            f = figure(figsize=(6.2, 5.6))
            ax = f.add_axes([0.17, 0.02, 0.72, 0.79])
            axcolor = f.add_axes([0.90, 0.02, 0.03, 0.79])
            
            im = ax.matshow(ekep, cmap='inferno',interpolation='none', norm=LogNorm(vmin=30, vmax=70))    
            t = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0,1.5,2,3,5,10]
            f.colorbar(im,cax=axcolor, ticks=t, format="$%.2f$")
            f.title('log'+tit)
            
            f.show()
        logkep(abcerror,'abc error')
        logkep(aberror,'ab error')
#%%
zoomabx,zoomaby,zoomabcx,zoomabcy,na,nb,naby,nabcy,delta=loading()
abfit,poptab=l2fit() #,abfit0_errors
abcfit,poptabc=l4fit() #,abcfit0_errors
nagymatrix, matrix,ini, x = readmap("E:/nano/2 peak fitting/anyagmegmondo/red_long_avg_correct.txt")  
non_graphene= finder(nagymatrix,11)
matrix=np.array(matrix)
for i, e in enumerate(matrix):
    matrix[i]= e/max(e)
#%%
ABC,cba,AB,ba=find_peaks(2200,2700,xb0= -1e3 ,yb0=1 ,w2=6,x0=-10,y0=0)
ABC = np.array(ABC, dtype=float);AB = np.array(AB, dtype=float);
#%%
kiiertekeltmap=[];ABerrmap=[];ABCerrmap=[]
szazalek,ABerrmap,ABCerrmap=mapcalc('E:/nano/2 peak fitting/anyagmegmondo/mp')
    #%%
picture(64,szazalek,ABerrmap,ABCerrmap,2200,2600,0,2)  










