"""
Example of the code used to analyze results from the coupled-cluster 
calculations and plots generated in the paper.

"""

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as patches
from scipy.interpolate import UnivariateSpline

plt.style.use('seaborn-talk')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#This function calculates error bars as a result of truncating the 
#effective field theory expansion.  Formula derived using bayesian analysis
# observable ('X') values, order-by-order LO - o - NLO - NNLO - N3LO - ...
#
#  X(p) = X0 \sum_{n=0}^{\infty} [c_n*Q^n]
#

def eft_error(X, kF, Lambda, conservative):
 
    # physical constants
    avg_pion_mass=138.039
 
    # insert 0.0 at index 1 in X list
    X.insert(1,0.0)
    order=['N0LO','----','N1LO','N2LO','N3LO','N4LO','N5LO']
     
    # expansion parameter
    Q=max(kF/Lambda,avg_pion_mass/Lambda)
    scale = X[0]
    if scale == 0.0:
        scale = 0.0000001
    conservative_error=[]
     
    for nu1 in range(0,len(X)):
        cons_max = 0.0
        if nu1==1:
            conservative_error.append(cons_max)
            continue
        for nu2 in range(nu1+1,len(X)):
            if nu2==1:
                continue
            cons_diff = abs(X[nu1]-X[nu2])
            if cons_diff > cons_max:
                cons_max = cons_diff
        conservative_error.append(cons_max)
         
    # extract expansion coefficients 'c'
    # initialize with constant nu=0,1 values. 
    c=[1.0,0.0]
    c_abs = [0,0.0]
    # exoected values from simple Q^nu scaling
    c_expected = [Q**0,0.0]
    # terms up to nu=1
    sub_terms = c[0]*Q**0 + c[1]*Q**1
 
    # compute errors
    sigma=[]
    #error on LO result
    sigma.append(Q**2*max(c))
    #append vanishing Q**1 order
    sigma.append(0)
 
    for nu in range(2,len(X)):
        c_nu = (X[nu]/scale - sub_terms)/Q**nu
        c.append(c_nu)
        c_abs.append(abs(c_nu))
        sub_terms = sub_terms + c[nu]*Q**nu
        c_expected.append(Q**nu)
        sigma.append(Q**(nu+1)*(max(c_abs)))
 
    error=[]
    for nu in range(0,len(X)):
        if conservative=='yes':
            error.append(max(conservative_error[nu],abs(sigma[nu]*X[0])))
        else:
            error.append(abs(sigma[nu]*X[0]))
             
    return error

#%%
# For the EFT-truncation error: this version requires to use LO-NLO-NNLO calcs with idential
# kF values. Otherwise, modify the code to use splines or other interpolations.
 
def plotting_with_error(files, figurefile,matter):
    plt.style.use('seaborn-talk')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
   # data_dir = "/home/pschwar3/deltafull-data/"
    titles = [r'LO', r'$\Delta$NLO',r'$\Delta$NNLO']
     
    Xcol=3 # rho
    kcol=2 # kF
    Ycol=7 # E/A
     
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)#, aspect='equal')
     
    # spline the data
    smooth = 'yes'
    # nof points
    Npts=100
    eps = 0.005
    for i, fn in enumerate(files):
        datafile = fn
        print("Reading {}".format(datafile))
        data = np.loadtxt(datafile,comments='#')
     
        X = data[:,Xcol]
        Y = data[:,Ycol]
        kF= data[:,kcol]
         
        if "NNLO" in datafile:
            lcolor = "blue"
            style = '-'
            X_NNLO = X
            Y_NNLO = Y
           # kF_NNLO = kF
        elif "NLO" in datafile:
            lcolor = "green"
            style = '-'
            X_NLO = X
            Y_NLO = Y
          #  kF_NLO = kF
        elif "LO" in datafile:
            lcolor = "red"
            style = '-'
            X_LO = X
            Y_LO = Y
            kF_LO = kF
     
        if smooth == 'no':
            plt.plot(X,Y,label=titles[i],ls=style,color=lcolor,linewidth=3.5)
        elif smooth == 'yes':
            #print('using splines')
            Xsmooth = np.linspace(X.min()-eps,X.max()+eps,Npts)
            spl = UnivariateSpline(X,Y)
            plt.plot(Xsmooth,spl(Xsmooth),ls=style,color=lcolor,linewidth=3.5)
 
    ndata=len(Y_LO)
    print('ndata: %d'%ndata)
    hbarc = 197.32697
    Lambda=500.0
    conservative='yes'
    kF_MeV=[]
    eLO=[]
    eNLO=[]
    eNNLO=[]
    print('Lambda=%d'%Lambda)
    for nd in range(0,ndata):
        kF_MeV.append(kF_LO[nd]*hbarc)
        EA_order_by_order=[Y_LO[nd],Y_NLO[nd],Y_NNLO[nd]]
        #print(EA_order_by_order)
        errors = eft_error(EA_order_by_order,kF_MeV[nd],Lambda,conservative)
        eLO.append(errors[0])
        eNLO.append(errors[2])
        eNNLO.append(errors[3])
        #print (kF[nd], eLO[nd], eNLO[nd], eNNLO[nd])
 
    if smooth == 'yes': 
        Xsmooth_LO = np.linspace(X_LO.min()-eps,X_LO.max()+eps,Npts)
        spl_Y_LO = UnivariateSpline(X_LO,Y_LO)
        spl_eLO = UnivariateSpline(X_LO,eLO)
         
        Xsmooth_NLO = np.linspace(X_NLO.min()-eps,X_NLO.max()+eps,Npts)
        spl_Y_NLO = UnivariateSpline(X_NLO,Y_NLO)
        spl_eNLO = UnivariateSpline(X_NLO,eNLO)
         
        Xsmooth_NNLO = np.linspace(X_NNLO.min()-eps,X_NNLO.max()+eps,Npts)
        spl_Y_NNLO = UnivariateSpline(X_NNLO,Y_NNLO,k=4)
        spl_eNNLO = UnivariateSpline(X_NNLO,eNNLO)
         
        plt.fill_between(Xsmooth_LO, spl_Y_LO(Xsmooth_LO)-spl_eLO(Xsmooth_LO), spl_Y_LO(Xsmooth_LO)+spl_eLO(Xsmooth_LO),facecolor='red',alpha=0.25,zorder=1,edgecolor='black',lw=1.0)
        plt.fill_between(Xsmooth_NLO, spl_Y_NLO(Xsmooth_NLO)-spl_eNLO(Xsmooth_NLO), spl_Y_NLO(Xsmooth_NLO)+spl_eNLO(Xsmooth_NLO),facecolor='green',alpha=0.25,zorder=2,edgecolor='black',lw=1.0)
        plt.fill_between(Xsmooth_NNLO, spl_Y_NNLO(Xsmooth_NNLO)-spl_eNNLO(Xsmooth_NNLO), spl_Y_NNLO(Xsmooth_NNLO)+spl_eNNLO(Xsmooth_NNLO),facecolor='blue',alpha=0.25,zorder=3,edgecolor='black',lw=1.0)
        plt.plot(X_LO,Y_LO,'-^r',label=titles[0])
        plt.plot(X_LO,Y_NLO,'-sg',label=titles[1])
        plt.plot(X_LO,Y_NNLO,'-ob',label=titles[2])
        
    else:
        plt.fill_between(X_LO, Y_LO-eLO, Y_LO+eLO,facecolor='red',alpha=0.25,zorder=1,edgecolor='black',lw=1.0)
        plt.fill_between(X_LO, Y_NLO-eNLO, Y_NLO+eNLO,facecolor='green',alpha=0.25,zorder=2,edgecolor='black',lw=1.0)
        plt.fill_between(X_LO, Y_NNLO-eNNLO, Y_NNLO+eNNLO,facecolor='blue',alpha=0.25,zorder=3,edgecolor='black',lw=1.0)
    
    
    ax1.add_patch(patches.Rectangle(
        (0.15, -16.5),   # (x,y) lower left corner
        0.02,            # width
        1.0,             # height
        color='black', fill=False, hatch='////', lw=2.0))
 
    plt.xlabel(r'$\rho$ [fm$^{-3}$]',fontsize=20)
    plt.ylabel(r'$E/A$ [MeV]',fontsize=20)
     
    plt.legend(loc="best", fancybox=False,frameon=False, fontsize=20,numpoints=2)
     
     
    axes = plt.gca()
    if matter == 'snm':
        spl_Y_deriv1 = spl_Y_NNLO.derivative()
        sat_pt = spl_Y_deriv1.roots()
        print(sat_pt)
        #plt.scatter(sat_pt,spl_Y_NNLO(sat_pt),s=200,c='b',marker='*')
        axes.set_xlim([0.05-eps,0.20+eps])
        axes.set_ylim([-18.0,-0.0])
    elif matter == 'pnm':
        axes.set_xlim([0.05-eps,0.20+eps])
        axes.set_ylim([0.0,25.0])
        
    axes.xaxis.set_tick_params(labelsize=20)
    axes.yaxis.set_tick_params(labelsize=20)
    fig1 = plt.gcf()
    
    plt.show()
    plt.draw()
    print("Saving to file:",figurefile)
    fig1.savefig(figurefile, format='pdf',bbox_inches='tight')
    return eLO, eNLO, eNNLO

#%%
#

inputfiles=["LO_nmax4_snm_450.dat", 
       "NLO_nmax4_snm_450_tnf3.dat",
       "NNLO_nmax4_snm_450_tnf3.dat"]

figfile = "nmax4_snm_450_tnf3.pdf"
eLO_snm_450 = []; eNLO_snm_450_tnf3 = []; eNNLO_snm_450_tnf3 = [];
eLO_snm_450, eNLO_snm_450_tnf3, eNNLO_snm_450_tnf3 = plotting_with_error(inputfiles,figfile,'snm')

inputfiles=["LO_nmax4_pnm_450.dat", 
       "NLO_nmax4_pnm_450_tnf3.dat",
       "NNLO_nmax4_pnm_450_tnf3.dat"]

figfile = "nmax4_pnm_450_tnf3.pdf"
eLO_pnm_450 = []; eNLO_pnm_450_tnf3 = []; eNNLO_pnm_450_tnf3 = [];
eLO_pnm_450, eNLO_pnm_450_tnf3, eNNLO_pnm_450_tnf3 = plotting_with_error(inputfiles,figfile,'pnm')




