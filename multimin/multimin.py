import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools
import math
import string
from time import time,strftime
from scipy.stats import multivariate_normal as multinorm
from scipy.optimize import minimize
import pickle
from hashlib import md5
import os

class Util(object):
    """
    This abstract class contains useful methods for the package.
    """ 
    
    #Mathematical functions
    """
    #Interesting but it can be problematic
    sin=math.sin
    cos=math.cos
    log=math.log
    exp=math.exp
    """
    sin=np.sin
    cos=np.cos
    log=np.log
    exp=np.exp
    
    #Stores the time of start of the script when gravray is imported
    TIMESTART=time()
    #Stores the time of the last call of elTime
    TIME=time()
    #Stores the duration between elTime consecutive calls 
    DTIME=-1
    DUTIME=[]
    
    #Lambda methods
    """
    This set of routines allows the conversion from a finite interval [0,s] to an unbound one [-inf,inf]
    
    Example:
    
        scales=[0,0,10,10,1]
        minparams=[0.0,0.0,1,1,0.7]
        uparams=Util.tIF(minparams,scales,Util.f2u)
        
        Result:
        [0.0, 0.0, -2.197224577336219, -2.197224577336219, 0.8472978603872034]
    """
    f2u=lambda x,s:Util.log((x/s)/(1-(x/s)))
    u2f=lambda t,s:s/(1+Util.exp(-t))
    tIF=lambda p,s,f:[f(p[i],s[i]) if s[i]>0 else p[i] for i in range(len(p))]

    def errorMsg(error,msg):
        """
        Add a custom message msg to an error handle.

        Parameters:
            error: error handle, handle (eg. except ValueError as error)
            msg: message to add to error, string.

        Return: None.
        """
        error.args=(error.args if error.args else tuple())+(msg,)

    def _tUnit(t):
        for unit,base in dict(d=86400,h=3600,min=60,s=1e0,ms=1e-3,us=1e-6,ns=1e-9).items():
            tu=t/base
            if tu>1:break
        return tu,unit,base
    
    def elTime(verbose=1,start=False):
        """
        Compute the time elapsed since last call of this routine.  The displayed time 
        is preseneted in the more convenient unit, ns (nano seconds), us (micro seconds), 
        ms (miliseconds), s (seconds), min (minutes), h (hours), d (days)

        Parameters: None.

        Optional:
            verbose: show the time in screen (default 1), integer or boolean.
            start: compute time from program start (deault 0), integer or boolean.

        Return: None.

        Examples:
            elTime(), basic usage (show output)
            elTime(0), no output
            elTime(start=True), measure elapsed time since program 
            print(DTIME,DUTIME), show values of elapsed time
        """
        t=time()
        dt=t-Util.TIME
        if start:
            dt=t-Util.TIMESTART    
            msg="since script start"
        else:
            msg="since last call"
        dtu,unit,base=Util._tUnit(dt)
        if verbose:print("Elapsed time %s: %g %s"%(msg,dtu,unit))
        Util.DTIME=dt
        Util.DUTIME=[dtu,unit]
        Util.TIME=time()
        return dt,[dtu,unit] 

    def mantisaExp(x):
        """Calculate the mantisa and exponent of a number.
        
        Parameters:
            x: number, float.
            
        Return:
            man: mantisa, float
            exp: exponent, float.
            
        Examples:
            m,e=mantisaExp(234.5), returns m=2.345, e=2
            m,e=mantisaExp(-0.000023213), return m=-2.3213, e=-5
        """
        xa=np.abs(x)
        s=np.sign(x)
        try:
            exp=int(np.floor(np.log10(xa)))
            man=s*xa/10**(exp)
        except OverflowError as e:
            man=exp=0
        return man,exp

class PlotGrid(object):
    """Class PlotGrid
    
    Create a grid of plots showing the projection of a N-dimensional
    
    Initialization attributes:
        dproperties: list of properties to be shown, dictionary of dictionaries (N entries)
            keys: label of attribute, ex. "q"
            dictionary: 
                label: label used in axis, string
                range: range for property, tuple (2)
        
    Optional initialization attributes:
        figsize=3 : base size for panels (the size of figure will be M x figsize), integer
        fontsize=10 : base fontsize, int
        direction='out' : direction of ticks in panels.
    
    Other attributes:
        N: number of properties, int
        M: size of grid matrix (M=N-1), int
        fw: figsize
        fs: fontsize
        fig: figure handle, figure
        axs: matrix with subplots, axes handles (MxM)
        axp: matrix with subplots, dictionary of dictionaries
        properties: list of properties labels, list of strings (N)
    
    Methods:
        tightLayout
        setLabels
        setRanges
        setTickParams
        
        plotHist
        scatterPlot
    """
    
    def __init__(self,properties,figsize=3,fontsize=10,direction='out'):

        #Basic attributes
        self.dproperties=properties
        self.properties=list(properties.keys())

        #Secondary attributes
        self.N=len(properties)
        self.M=self.N-1
        
        #Optional properties
        self.fw=figsize
        self.fs=fontsize

        #Create figure and axes: it works
        try:
            self.fig,self.axs=plt.subplots(
                self.M,self.M,
                constrained_layout=True,
                figsize=(self.M*self.fw,self.M*self.fw),
                sharex="col",sharey="row"
            )
            self.constrained=True
        except:
            self.fig,self.axs=plt.subplots(
                self.M,self.M,
                figsize=(self.M*self.fw,self.M*self.fw),
                sharex="col",sharey="row"
            )
            self.constrained=False

        if not isinstance(self.axs,np.ndarray):
            self.axs=np.array([[self.axs]])
            self.single = True
        else:
            self.single = False

        #Create named axis
        self.axp=dict()
        for j in range(self.N):
            propj=self.properties[j]
            if propj not in self.axp.keys():
                self.axp[propj]=dict()
            for i in range(self.N):
                propi=self.properties[i]
                if i==j:
                    continue
                if propi not in self.axp.keys():
                    self.axp[propi]=dict()
                if i<j:
                    self.axp[propj][propi]=self.axp[propi][propj]
                    continue
                self.axp[propj][propi]=self.axs[i-1][j]
    
        #Deactivate unused panels
        for i in range(self.M):
            for j in range(i+1,self.M):
                self.axs[i][j].axis("off")
        
        #Place ticks
        for i in range(self.M):
            for j in range(i+1):
                if not self.single:
                    self.axs[i,j].tick_params(axis='both',direction=direction)
                else:
                    self.axs[i,i].tick_params(axis='both',direction=direction)
        for i in range(self.M):
            self.axs[i,0].tick_params(axis='y',direction="out")
            self.axs[self.M-1,i].tick_params(axis='x',direction="out")
        
        #Set properties of panels
        self.setLabels()
        self.setRanges()
        self.setTickParams()
        self.tightLayout()
    
    def tightLayout(self):
        """
        Tight layout if no constrained_layout was used.
        
        Parameters: None
        
        Return: None
        """
        if self.constrained==False:
            self.fig.subplots_adjust(wspace=self.fw/100.,hspace=self.fw/100.)
        self.fig.tight_layout()
        
    def setTickParams(self,**args):
        """
        Set tick parameters.
        
        Parameters: 
            **args: same arguments as tick_params method, dictionary
        
        Return: None
        """
        opts=dict(axis='both',which='major',labelsize=0.8*self.fs)
        opts.update(args)
        for i in range(self.M):
            for j in range(self.M):
                self.axs[i][j].tick_params(**opts)
        
    def setRanges(self):
        """
        Set ranges in panels according to ranges defined in dparameters.
        
        Parameters: None
        
        Return: None
        """
        for i,propi in enumerate(self.properties):
            for j,propj in enumerate(self.properties):
                if j<=i:continue
                if self.dproperties[propi]["range"] is not None:
                    self.axp[propi][propj].set_xlim(self.dproperties[propi]["range"])
                if self.dproperties[propj]["range"] is not None:
                    self.axp[propi][propj].set_ylim(self.dproperties[propj]["range"])
    
    def setLabels(self,**args):
        """
        Set labels parameters.
        
        Parameters: 
            **args: common arguments of set_xlabel, set_ylabel and text, dictionary
        
        Return: None
        """
        opts=dict(fontsize=self.fs)
        opts.update(args)
        for i,prop in enumerate(self.properties[:-1]):
            label=self.dproperties[prop]["label"]
            self.axs[self.M-1][i].set_xlabel(label,**opts)
        for i,prop in enumerate(self.properties[1:]):
            label=self.dproperties[prop]["label"]
            self.axs[i][0].set_ylabel(label,rotation=0,labelpad=10,**opts)

        for i in range(1,self.M):
            label=self.dproperties[self.properties[i]]["label"]
            self.axs[i-1][i].text(0.5,0.0,label,ha='center',
                                  transform=self.axs[i-1][i].transAxes,**opts)
            #270 if you want rotation
            self.axs[i-1][i].text(0.0,0.5,label,rotation=0,va='center',
                                  transform=self.axs[i-1][i].transAxes,**opts)

        label=self.dproperties[self.properties[0]]["label"]
        if not self.single:
            self.axs[0][1].text(0.0,1.0,label,rotation=0,ha='left',va='top',
                                transform=self.axs[0][1].transAxes,**opts)
        
        label=self.dproperties[self.properties[-1]]["label"]
        #270 if you want rotation
        self.axs[-1][-1].text(1.05,0.5,label,rotation=0,ha='left',va='center',
                              transform=self.axs[-1][-1].transAxes,**opts)

        self.tightLayout()
        
    def plotHist(self,data,colorbar=False,**args):
        """
        Create a 2d-histograms of data on all panels of the PlotGrid.
        
        Parameters: 
            data: data to be histogramed (n=len(data)), numpy array (nxN)
            
        Optional parameters:
            colorbar=False: include a colorbar?, boolean or int (0/1)
            **args: all arguments of hist2d method, dictionary
        
        Return: 
            hist: list of histogram instances.
        """
        opts=dict()
        opts.update(args)
            
        hist=[]
        for i,propi in enumerate(self.properties):
            if self.dproperties[propi]["range"] is not None:
                xmin,xmax=self.dproperties[propi]["range"]
            else:
                xmin=data[:,i].min()
                xmax=data[:,i].max()
            for j,propj in enumerate(self.properties):
                if j<=i:continue
                    
                if self.dproperties[propj]["range"] is not None:
                    ymin,ymax=self.dproperties[propj]["range"]
                else:
                    ymin=data[:,j].min()
                    ymax=data[:,j].max()                
                
                opts["range"]=[[xmin,xmax],[ymin,ymax]]
                h,xe,ye,im=self.axp[propi][propj].hist2d(data[:,i],data[:,j],**opts)
                
                hist+=[im]
                if colorbar:
                    #Create color bar
                    divider=make_axes_locatable(self.axp[propi][propj])
                    cax=divider.append_axes("top",size="9%",pad=0.1)
                    self.fig.add_axes(cax)
                    cticks=np.linspace(h.min(),h.max(),10)[2:-1]
                    self.fig.colorbar(im,
                                      ax=self.axp[propi][propj],
                                      cax=cax,
                                      orientation="horizontal",
                                      ticks=cticks)
                    cax.xaxis.set_tick_params(labelsize=0.5*self.fs,direction="in",pad=-0.8*self.fs)
                    xt=cax.get_xticks()
                    xm=xt.mean()
                    m,e=Util.mantisaExp(xm)
                    xtl=[]
                    for x in xt:
                        xtl+=["%.1f"%(x/10**e)]
                    cax.set_xticklabels(xtl)
                    cax.text(0,0.5,r"$\times 10^{%d}$"%e,ha="left",va="center",
                             transform=cax.transAxes,fontsize=6,color='w')

        self.setLabels()
        self.setRanges()
        self.setTickParams()
        self.tightLayout()
        return hist
                    
    def scatterPlot(self,data,**args):
        """
        Scatter plot on all panels of the PlotGrid.
        
        Parameters: 
            data: data to be histogramed (n=len(data)), numpy array (nxN)
            
        Optional parameters:
            **args: all arguments of scatter method, dictionary
        
        Return: 
            scatter: list of scatter instances.
        """
        scatter=[]
        for i,propi in enumerate(self.properties):
            for j,propj in enumerate(self.properties):
                if j<=i:continue
                scatter+=[self.axp[propi][propj].scatter(data[:,i],data[:,j],**args)]

        self.setLabels()
        self.setRanges()
        self.setTickParams()
        self.tightLayout()
        return scatter
    
    def SaveFig(self, file_name=False):
        if file_name:
            self.fig.savefig(file_name)
        else:
            file_name = 'multiminfig.png'
            self.fig.savefig(file_name)

class UtilStats(object):
    """
    Abstract class with useful routines
    """
    #Golden ratio: required for golden gaussian.
    phi=(1+5**0.5)/2
    
    def genIndex(probs):
        """
        Given a set of (normalized) probabilities, randomly generate an index n following the 
        probabilities.

        For instance if we have 3 events with probabilities 0.1, 0.7, 0.2, genIndex will generate
        a number in the set (0,1,2) having those probabilities, ie. 1 will have 70% of probability.
        
        Parameters:
            probs: Probabilities, numpy array (N), adimensional
                NOTE: It should be normalized, ie. sum(probs)=1
            
        Return:
            n: Index in the set [0,1,2,... len(probs)-1], integer
            
        Example:
            genIndex([0.1,0.7,0.2])
            
        Used in:
            - ComposedMultiVariateNormal.rvs
        """
        cums=np.cumsum(probs)
        if not math.isclose(cums[-1],1,rel_tol=1e-5):
            raise ValueError("Probabilities must be normalized, ie. sum(probs) = 1")
        cond=(np.random.rand()-cums)<0
        isort=np.arange(len(probs))
        n=isort[cond][0] if sum(cond)>0 else isort[0]
        return n
        
    def setMatrixOffDiagonal(M,off):
        """
        Set a matrix with the terms of the off diagonal

        Parameters:
            M: Matrix, array, n x n 
            off: Terms off diagonal, array, n x (n-1) / 2

        Returns:
            Implicitly the matrix M has now the off diagonal terms

        Example:
            M=np.eye(3)
            off=[0.1,0.2,0.3]
            setMatrixOffDiagonal(M,off)

            Result:
            M=array([[1. , 0.1, 0.2],
                       [0.1, 1. , 0.3],
                       [0.2, 0.3, 1. ]])
                       
        Sources:
            https://newbedev.com/how-to-get-indices-of-non-diagonal-elements-of-a-numpy-array
        """
        I,J=np.where(~np.eye(M.shape[0],dtype=bool))
        ffo=list(off[::-1])
        for i,j in zip(I,J):M[i,j]=ffo.pop() if j>i else 0
        M[:,:]=np.triu(M)+np.tril(M.T,-1)
        
    def calcCovarianceFromCorrelations(sigmas,rhos):
        """
        Compute covariance matrices from the stds and correlations (rho)

        Parameters:

            sigmas: Array of values of standard deviation for variables, 
                    array, Ngauss x Nvars

            rhos: Array with correlations, array Ngauss x Nvars x Nvars

        Returns:

            Sigmas: Array with covariance matrices corresponding to these sigmas 
                    and rhos, array, Ngauss x Nvars x Nvars

        Examples:
            sigmas=[[1,2,3]]
            # rho_12, rho_13, rho_23
            rhos=[[0.1,0.2,0.3]] #Size of rhos is Nvars x (Nvars-1)/2
            UtilStats.calcCovarianceFromCorrelations(sigmas,rhos)
            
            Results:
            array([[[1. , 0.2, 0.6],
                    [0.2, 4. , 1.8],
                    [0.6, 1.8, 9. ]]])
                    
            This is equivalent to:
            
            rho=rhos[0]
            sigma=sigmas[0]
            R=np.eye(3)
            UtilStats.setMatrixOffDiagonal(R,rho)
            M=np.zeros((3,3))
            for i in range(3):
                for j in range(3):
                    M[i,j]=R[i,j]*sigma[i]*sigma[j]
            
        Sources: 

            Based on: 
            https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/    
        """
        try:
            Nvars=len(sigmas[0])        
        except:
            raise AssertionError("Array of sigmas must be an array of arrays")
        try:
            Nrhos=len(rhos[0])        
        except:
            raise AssertionError("Array of rhos must be an array of arrays")
            
        if Nrhos!=int(Nvars*(Nvars-1)/2):
            raise AssertionError(f"Size of rhos ({Nrhos}) are incompatible with Nvars={Nvars}.  It should be Nvars(Nvars-1)/2={int(Nvars*(Nvars-1)/2)}.")
            
        Sigmas=np.array(len(sigmas)*[np.eye(Nvars)])
        for Sigma,sigma,rho in zip(Sigmas,sigmas,rhos):
            UtilStats.setMatrixOffDiagonal(Sigma,rho)
            Sigma*=np.outer(sigma,sigma)
        return Sigmas

    def calcCorrelationsFromCovariances(Sigmas):
        """
        Compute the standard deviations and corresponding correlation coefficients given a set of 
        covariance matrices.

        Parameters:
            Sigmas: Array of covariance matrices, array, Ngauss x Nvars x Nvars

        Returns:
            sigmas: Array of standard deviarions, array, Ngauss x Nvars
            rhos: Array of correlation coefficiones, array, Ngauss x Nvars (Nvars-1) / 2

        Example: 
            Sigmas=[
                    [[1. , 0.2, 0.6],
                    [0.2, 4. , 1.8],
                    [0.6, 1.8, 9. ]]
                   ]
            calcCorrelationsFromCovariances(Sigmas)

            Result:

            sigmas=array([1., 2., 3.])
            rhos=[[0.1, 0.20000000000000004, 0.3]]        
        """
        if len(np.array(Sigmas).shape)!=3:
            raise AssertionError(f"Array of Sigmas (shape {np.array(Sigmas).shape}) must be an array of matrices")

        sigmas=[]
        rhos=[]
        for n,Sigma in enumerate(np.array(Sigmas)):
            sigmas+=[(np.diag(Sigma))**0.5]
            R=Sigma/np.outer(sigmas[n],sigmas[n])
            I,J=np.where(~np.eye(R.shape[0],dtype=bool))
            rhos+=[[]]
            for i,j in zip(I,J):rhos[n]+=[R[i,j]] if j>i else []
        return np.array(sigmas),np.array(rhos)    

    def calcCovarianceFromRotation(sigmas,angles):
        """
        Compute covariance matrices from the stds and the angles.
    
        Parameters:
        
            sigmas: Array of values of standard deviation for variables, 
                    array Ngauss x 3
                    
            angles: Euler angles expressing the directions of the principal 
                    axes of the distribution, array, Ngauss x 3

        Returns:
        
            Sigmas: Array with covariance matrices corresponding to these sigmas 
                    and angles, array, Ngauss x 3 x 3

        Examples:
            

        Sources: 
        
            Based on: 
            https://www.visiondummy.com/2014/04/geometric-interpretation-covariance-matrix/    
        """
        try:
            Nvars=len(sigmas[0])        
        except:
            raise AssertionError("Sigmas must be an array of arrays")
        Sigmas=[]
        for scale,angle in zip(sigmas,angles):
            L=np.identity(Nvars)*np.outer(np.ones(Nvars),scale)
            Rot=spy.eul2m(-angle[0],-angle[1],-angle[2],3,1,3) if Nvars==3 else spy.rotate(-angle[0],3)[:2,:2]
            Sigmas+=[np.matmul(np.matmul(Rot,np.matmul(L,L)),np.linalg.inv(Rot))]

        return np.array(Sigmas)

    def flattenSymmetricMatrix(M):
        """
        Given a symmetric matrix the routine returns the flatten version of the Matrix.

        Parameters:
            M: Matrix, array, n x n 

        Returns:
            F: Flatten array, array, nx(n+1)/2

        Example:
            M=[[1,0.2],[0.2,3]]
            F=flattenSymmetricMatrix(M)

            Result:
            F=[1,0.2,3]
        """
        return M[np.triu_indices(M.shape[0], k = 0)]

    def unflattenSymmetricMatrix(F,M):
        """
        Given a flatten version of a matrix, returns the symmetric matrix.

        Parameters:
            F: Flatten array, array, n x (n+1)/2
            M: Matrix where the result will be stored, array, n x n 

        Returns:
            It return the results in matrix M.

        Example:
            F = [1,0.2,3]
            M=np.zeros((2,2))
            unflattenSymmetricMatrix(F,M)

            Results:
            M=[[1,0.2],[0.2,3]]
            
        Notes:
            Number of components:

                Given f (size of flatten matrix) what is the dimension (N) of the symmetric matrix:

                $$
                N(N+1)/2 = f
                $$

                $$
                N^2+N-2f=0
                $$

                $$
                N=(\sqrt{1+8f}-1)/2
                $$
        """
        M[np.triu_indices(M.shape[0],k=0)]=np.array(F)
        M[:,:]=np.triu(M)+np.tril(M.T,-1)

class ComposedMultiVariateNormal(object):
    """
    A linear combination of multivariate normal distribution (MND) with speccial methods 
    for specifying the parameters of the distributions.
    
    Basic attributes:
        Ngauss: Number of composed MND, int
        
        Nvars: Number of random variables, int
        
        mus: Array with average (mu) of random variables, array, Ngauss x Nvars

        weigths: Array with weights of each MND, array, Ngauss
            
            NOTE: This weights are normalized at the end.
             
        sigmas: Standard deviation of each variable, array, Ngauss x Nvars
        
        rhos: Elements of the upper triangle of the correlation matrix, array, Ngauss x Nvars x (Nvars-1)/2
        
        Sigmas: Array with covariance matrices for each MND, array, Ngauss x Nvars x Nvars
        
        params: Parameters of the distribution in flatten form including symmetric elements of the covariance
                matrix, array, Ngauss*(1+Nvars+Nvars*(Nvars+1)/2)
        
        stdcorr: Parameters of the distribution in flatten form, including upper triangle of the correlation
                 matrix, array, Ngauss*(1+Nvars+Nvars*(Nvars+1)/2)
        
    Initialization:
        There are several ways of initialize a CMND:
        
        Providing: Ngauss and Nvars
            In this case the class is instantiated with zero means, unitary dispersion and 
            covariance matrix equal to Ngasus identity matrices Nvars x Nvars.
            
        Providing: params, Nvars
            In this case you have a flatted version of the parametes (weights, mus, Sigmas)
            and want to instantiate the system.  All parameters are set and no other action
            is required.

        Providing: stdcorr, Nvars
            In this case you have a flatted version of the parametes (weights, mus, sigmas, rhos)
            and want to instantiate the system.  All parameters are set and no other action
            is required.
            
        Providing: weights, mus, Sigmas (optional)
            In this case the basic properties of the CMND are set.
            
    Methods:
        pdf(X):
            Parameters:
                X: Set of values of the random variable, array Nvars

            Return:
                pdf: Value of the composed pdf at X.
                
        rvs(Nsam):
            Parameters:
                Nsam: Number of samples drawn from CMND
            
            Return:
                Xs: Array with Nsam samples, array Nsam x Nvars
            
    Examples:
    
        mus=[[0,0],[1,1]]
        weights=[0.1,0]
        MND1=ComposedMultiVariateNormal(mus=mus,weights=weights)
        MND1.setSigmas(
            [
                [[1,0.2],[0,1]],
                [[1,0],[0,1]]
            ]
        )
        print(MND1)

        params=[0.1,0.9,0,0,1,1,1,0.2,0.2,1,1,0,0,1]
        MND2=ComposedMultiVariateNormal(params=params,2)        
        print(MND2)
    
    """
    
    #Control behavior
    _ignoreWarnings=True
        
    def __init__(self,
                 Ngauss=0,
                 Nvars=0,
                 params=None,
                 stdcorr=None,
                 weights=None,mus=None,Sigmas=None):
        
        #Method 1: initialize a simple instance
        if Ngauss>0:
            mus=[[0]*Nvars]*Ngauss
            weights=[1/Ngauss]*Ngauss
            Sigmas=[np.eye(Nvars)]*Ngauss
            self.__init__(mus=mus,weights=weights,Sigmas=Sigmas)
        
        #Method 2: initialize from flatten parameters
        elif params is not None:
            self.setParams(params,Nvars)

        #Method 3: initialize from flatten parameters
        elif stdcorr is not None:
            self.setStdcorr(stdcorr,Nvars)

        #Method 4: initialize from explicit arrays
        else:
            #Basic attributes
            mus=np.array(mus)
            try:
                mus[0,0]
            except Exception as e:
                Util.errorMsg(e,"Parameter 'mus' must be a matrix, eg. mus=[[0,0]]")
                raise
            self.mus=mus

            #Number of variables 
            self.Ngauss=len(mus)
            self.Nvars=len(mus[0])

            #Weights and normalization
            if weights is None:
                self.weights=[1]+(self.Ngauss-1)*[0]
            elif len(weights)!=self.Ngauss:
                raise ValueError(f"Length of weights array ({len(weights)}) must be equal to number of MND ({self.Ngauss})")
            else:
                self._normalizeWeights(weights)

            #Secondary attributes
            if Sigmas is None:
                self.Sigmas=None
                self.params=None
            else:
                self.setSigmas(Sigmas)
                
        self._nerror=0

    def setSigmas(self,Sigmas):
        """
        Set the value of list of covariance matrices.
        
        After setting Sigmas it update params and stdcorr.
        """
        self.Sigmas=np.array(Sigmas)
        self._checkSigmas()
        self._flattenParams()
        self._flattenStdcorr()
        
    def setParams(self,params,Nvars):
        """
        Set the properties of the CMND from flatten params.
        
        After setting it generate flattend stdcorr and normalize weights.
        """
        if Nvars==0 or len(params)==0:
            raise ValueError(f"When setting from flat params, Nvars ({Nvars}) cannot be zero")
        self._unflattenParams(params,Nvars)
        self._normalizeWeights(self.weights)
        return 

    def setStdcorr(self,stdcorr,Nvars):
        """
        Set the properties of the CMND from flatten stdcorr.
        
        After setting it generate flattened params and normalize weights.
        """
        if Nvars==0 or len(stdcorr)==0:
            raise ValueError(f"When setting from flat params, Nvars ({Nvars}) cannot be zero")
        self._unflattenStdcorr(stdcorr,Nvars)
        self._normalizeWeights(self.weights)
        return
    
    def _normalizeWeights(self,weights):
        """
        Normalize weights in such a way that sum(weights)=1
        """
        self.weights=np.array(weights)/sum(np.array(weights))
            
    def _flattenParams(self):
        """
        Flatten params
        """
        self._checkParams(self.Sigmas)
        
        #Flatten covariance matrix
        SF=[UtilStats.flattenSymmetricMatrix(self.Sigmas[i]).tolist() for i in range(self.Ngauss)]
        self.params=np.concatenate((self.weights.flatten(),
                                    self.mus.flatten(),
                                    list(itertools.chain(*SF))))
        self.Npars=len(self.params) #Ngauss*(1+Nvars+Nvar*(Nvars+1)/2)
        
    def _flattenStdcorr(self):
        """
        Flatten stdcorr
        """
        self._checkParams(self.sigmas)
        
        #Flatten stds. and correlations
        self.stdcorr=np.concatenate((self.weights.flatten(),
                                     self.mus.flatten(),
                                     self.sigmas.flatten(),
                                     self.rhos.flatten()
                                    ))
        self.Ncor=len(self.stdcorr)
        
    def _unflattenParams(self,params,Nvars):
        """
        Unflatten properties from params
        """

        self.params=np.array(params)
        self.Npars=len(self.params)
        
        factor=int(1+Nvars+Nvars*(Nvars+1)/2)
        
        if (self.Npars%factor)!=0:
            raise AssertionError(f"The number of parameters {self.Npars} is incompatible with the provided number of variables ({Nvars})")

        #Number of gaussian functions
        Ngauss=int(self.Npars/factor)

        #Get the weights
        i=0
        weights=self.params[i:Ngauss]
        i+=Ngauss
        
        #Get the mus
        mus=self.params[i:i+Ngauss*Nvars].reshape(Ngauss,Nvars)
        i+=Ngauss*Nvars
        
        #Get the sigmas
        Nsym=int(Nvars*(Nvars+1)/2)
        Sigmas=np.zeros((Ngauss,Nvars,Nvars))
        [UtilStats.unflattenSymmetricMatrix(F,Sigmas[i]) for i,F in enumerate(self.params[i:i+Ngauss*Nsym].reshape(Ngauss,Nsym))]

        #Normalize weights
        self._normalizeWeights(weights)

        #Check Sigmas
        self.Nvars=Nvars
        self.Ngauss=Ngauss
        self.weights=weights
        self.mus=mus
        self.Sigmas=Sigmas
        self._checkSigmas()

        #Flatten correlations
        self._flattenStdcorr()
        
    def _unflattenStdcorr(self,stdcorr,Nvars):
        """
        Unflatten properties from stdcorr
        """

        self.stdcorr=np.array(stdcorr)
        self.Ncor=len(self.stdcorr)

        factor=int(1+Nvars+Nvars*(Nvars+1)/2)

        if (self.Ncor%factor)!=0:
            raise AssertionError(f"The number of parameters {self.Ncor} is incompatible with the provided number of variables ({Nvars})")

        #Number of gaussian functions
        Ngauss=int(self.Ncor/factor)

        #Get the weights
        i=0
        weights=self.stdcorr[i:Ngauss]
        i+=Ngauss

        #Get the mus
        mus=self.stdcorr[i:i+Ngauss*Nvars].reshape(Ngauss,Nvars)
        i+=Ngauss*Nvars

        #Get the sigmas
        sigmas=self.stdcorr[i:i+Ngauss*Nvars].reshape(Ngauss,Nvars)
        i+=Ngauss*Nvars

        #Get the rhos
        Noff=int(Nvars*(Nvars-1)/2)
        rhos=self.stdcorr[i:i+Ngauss*Noff].reshape(Ngauss,Noff)

        #Normalize weights
        self._normalizeWeights(weights)

        #Set properties
        self.Nvars=Nvars
        self.Ngauss=Ngauss
        self.weights=weights
        self.mus=mus
        self.sigmas=sigmas
        self.rhos=rhos

        #Generate Sigma
        self.Sigmas=UtilStats.calcCovarianceFromCorrelations(self.sigmas,self.rhos)
        self._checkSigmas()

        #Flatten params
        self._flattenParams()

    def _checkSigmas(self):
        """
        Check value of sigmas
        """
        self._checkParams(self.Sigmas)
        
        #Check matrix
        if len(self.Sigmas)!=self.Ngauss:
            raise ValueError(f"You provided {len(self.Sigmas)} matrix, but Ngauss={self.Ngauss} are required")
        
        elif self.Sigmas[0].shape!=(self.Nvars,self.Nvars):
            raise ValueError(f"Matrices have wrong dimensions ({self.Sigmas[0].shape}). It should be {self.Nvars}x{self.Nvars}")
                
        #Symmetrize
        for i in range(self.Ngauss):
            self.Sigmas[i]=np.triu(self.Sigmas[i])+np.tril(self.Sigmas[i].T,-1)
            """
            #This check can be done, but it can be a problem when fitting
            if not np.all(np.linalg.eigvals(self.Sigmas[i])>0):
                raise ValueError(f"Matrix {i+1}, {self.Sigmas[i].tolist()} is not positive semidefinite.")
            """

        #Get sigmas and correlations
        self.sigmas,self.rhos=UtilStats.calcCorrelationsFromCovariances(self.Sigmas)

    def _checkParams(self,checkvar=None):
        if checkvar is None:
            raise AssertionError("You must first set the parameters (Sigmas, mus, etc.)")

    def pdf(self,X):
        """
        Compute the PDF.
        
        Parameter:
            X: point in the Nvar-dimensional space, numpy array (Nvar)
        
        Return:
            p: pdf.
        """
        self._checkParams(self.params)
        self._nerror=0
        value=0
        
        for w,muvec,Sigma in zip(self.weights,self.mus,self.Sigmas):
            try:
                value+=w*multinorm.pdf(X,muvec,Sigma)
            except Exception as error:
                if not self._ignoreWarnings:
                    print(f"Error: {error}, params = {self.params.tolist()}, stdcorr = {self.params.tolist()}")
                    self._nerror+=1
                value+=0
        return value
            
    def rvs(self,Nsam=1):
        """
        Generate a random sample of points following this Multivariate distribution.
        
        Parameter:
            Nsam: number of samples.
            
        Return:
            rs: samples, numpy array (Nsam x Nvars)
        """
        self._checkParams(self.params)
        
        Xs=np.zeros((Nsam,self.Nvars))
        for i in range(Nsam):
            n=UtilStats.genIndex(self.weights)
            Xs[i]=multinorm.rvs(self.mus[n],self.Sigmas[n])
        return Xs

    def sampleCMNDLikelihood(self,uparams,data=None,pmap=None,tset="stdcorr",scales=[],verbose=0):
        """
        Compute the negative value of the logarithm of the likelihood of a sample.
        
        Parameters:
            uparams: Minimization parameters (unbound), array
        
        Optional parameters
            data = None: data for which logL is computed.
            
            pmap = None: routine to map from minparams to params or stdcorr.
            
                Example:
                    def pmap(minparams):
                        stdcorr=np.array([1]+list(minparams))
                        stdcorr[-1:]-=1
                        return stdcorr
            
            tset = "stdcorr": type of minimization parameters. Values "params", "stdcorr"
            
            scales = []: list of scales for transforming uparams (unbound) in minparams (natural scale).
            
            verbose = 0: verbosity level (0, none, 1: input parameters, 2: full definition of the CMND)
        
        """
        #Map unbound minimization parameters into their right range
        minparams=np.array(Util.tIF(uparams,scales,Util.u2f))

        #Map minimizaiton parameters into CMND parameters
        params=np.array(pmap(minparams))

        if verbose>=1:
            print("*"*80)
            print(f"Minimization parameters: {minparams.tolist()}")
            print(f"CMND parameters: {params.tolist()}")

        #Update CMND parameters according to type of minimization parameters
        if tset=="params":self.setParams(params,self.Nvars)
        else:self.setStdcorr(params,self.Nvars)

        if verbose>=2:
            print("CMND:")
            print(self)

        #Compute PDF for each point in data and sum
        logL=-np.log(self.pdf(data)).sum()
        
        if verbose>=1:
            print(f"-logL = {logL:e}")

        return logL

    def plotSample(self,data=None,
                   N=10000,
                   props=None,ranges=None,
                   figsize=2,sargs=dict(),hargs=None):
        """
        Plot a sample of the CMND.
        
        Parameters:
            data = None: Data to plot.  If None it generate a sample.

            N = 10000: Number of points to generate the sample.
            
            props = None: Array with the name of the properties.  Ex. props=["x","y"]
            
            ranges = None: Array of ranges of the properties. Ex. ranges=[-3,3],[-5,5]
            
            figsize = 2: size of each axis. 
            
            hargs = None: Dictionary with options for the hist2d function.  Ex. hargs=dict(bins=50)
            
            sargs = dict(): Dictionary with options for the scatter plot.  Ex. sargs=dict(color='r')
            
        Returns: 
            G: Graphic handle. If Nvars = 2, it is a figure object, otherwise is a PlotGrid instance.
            
            
        Examples:
        
            G=CMND.plotSample(N=10000,sargs=dict(s=1,c='r'))
            G=CMND.plotSample(N=1000,sargs=dict(s=1,c='r'),hargs=dict(bins=20))


            CMND=ComposedMultiVariateNormal(Ngauss=1,Nvars=2)
            fig=CMND.plotSample(N=1000,hargs=dict(bins=20),sargs=dict(s=1,c='r'));

            CMND=ComposedMultiVariateNormal(Ngauss=2,Nvars=3)
            print(CMND)
            mus=[[0,0],[1,1]]
            weights=[0.1,0.9]
            Sigmas=[[[1,0.2],[0,1]],[[1,0],[0,1]]]
            MND1=ComposedMultiVariateNormal(mus=mus,weights=weights,Sigmas=Sigmas)
            #MND1=ComposedMultiVariateNormal(mus=mus,weights=weights);MND1.setSigmas(Sigmas)
            print(MND1)
            print(MND1.pdf([1,1]))
            params=[0.1, 0.9, 0.0, 0.0, 1.0, 1.0, 1.0, 0.2, 1.0, 1.0, 0.0, 1.0]
            MND2=ComposedMultiVariateNormal(params=params,Nvars=2)        
            print(MND2)
            print(MND2.pdf([1,1]))

        """
        if data is None:
            self.data=self.rvs(N)
        else:
            self.data=np.copy(data)
            
        properties=dict()
        for i in range(self.Nvars):
            symbol=string.ascii_letters[i] if props is None else props[i]
            rang=None if ranges is None else ranges[i]
            properties[symbol]=dict(label=f"${symbol}$",range=rang)
            
        if self.Nvars>2:
            G=PlotGrid(properties,figsize=figsize)
            if hargs is not None:
                G.plotHist(self.data,**hargs)
            G.scatterPlot(self.data,**sargs);
            G.fig.tight_layout()
            return G
        else:
            keys=list(properties.keys())
            fig=plt.figure(figsize=(5,5))
            ax=fig.gca()
            if hargs is not None:
                ax.hist2d(self.data[:,0],self.data[:,1],**hargs)
            #Experimental
            #sns.kdeplot(x=data[:,0],y=data[:,1],shade=True,ax=ax)
            ax.scatter(self.data[:,0],self.data[:,1],**sargs)
            ax.set_xlabel(properties[keys[0]]["label"])
            ax.set_ylabel(properties[keys[1]]["label"])
            fig.tight_layout()
            return fig
    
    def _strParams(self):
        """
        Generate strings explaining which quantities are stored in the flatten arrays params and stdcorr.
        
        It also generate and aray with the bounds applicable to the stdcorr parameters for purposes
        of minimization.
        
        Returns:
            None.
            
            Set the value of:
                str_params: list of properties in params, str
                str_stdcorr: list of properties in stdcorr, str.
                bnd_stdcorr: bounds of properties in stdcorr applicable for transforming to unbound.
        """
        str_params="["
        bnd_stdcorr="["
        #Probabilities
        for n in range(self.Ngauss):
            str_params+=f"p{n+1},"
            bnd_stdcorr+=f"1,"

        #Mus
        for n in range(self.Ngauss):
            for i in range(self.Nvars):
                str_params+=f"μ{n+1}_{i+1},"
                bnd_stdcorr+=f"0,"

        str_stdcorr=str_params
        #Std. devs
        for n in range(self.Ngauss):
            for i in range(self.Nvars):
                str_stdcorr+=f"σ{n+1}_{i+1},"
                bnd_stdcorr+=f"10,"

        #Sigmas
        for n in range(self.Ngauss):
            for i in range(self.Nvars):
                for j in range(self.Nvars):
                    if j>=i:
                        str_params+=f"Σ{n+1}_{i+1}{j+1},"
                    if j>i:
                        str_stdcorr+=f"ρ{n+1}_{i+1}{j+1},"
                        bnd_stdcorr+=f"2,"

        self.str_params=str_params.strip(",")+"]"
        self.str_stdcorr=str_stdcorr.strip(",")+"]"
        self.bnd_stdcorr=bnd_stdcorr.strip(",")+"]"

    def __str__(self):
        """
        Generates a string version of the object
        """
        self.Npars=len(self.params)
        self.Ncor=len(self.stdcorr)
        
        self._strParams()
        
        msg=f"""Composition of Ngauss = {self.Ngauss} gaussian multivariates of Nvars = {self.Nvars} random variables:
    Weights: {self.weights.tolist()}
    Number of variables: {self.Nvars}
    Averages (μ): {self.mus.tolist()}
"""
        if self.Sigmas is None:
            msg+=f"""    Sigmas: (Not defined yet)
    Params: (Not defined yet)
"""
        else:
            msg+=f"""    Standard deviations (σ): {self.sigmas.tolist()}
    Correlation coefficients (ρ): {self.rhos.tolist()}

    Covariant matrices (Σ): 
        {self.Sigmas.tolist()}
    Flatten parameters: 
        With covariance matrix ({self.Npars}):
            {self.str_params}
            {self.params.tolist()}
        With std. and correlations ({self.Ncor}):
            {self.str_stdcorr}
            {self.stdcorr.tolist()}"""  
        return msg

class FitCMND():
    """
    CMND Fitting handler
    
    Attributes:
        
        Ngauss: Number of fitting MND, int.
        
        Nvars: Number of variables in each MND, int.
        
        cmnd: Fitting object of the class CMND, ComposedMultiVariateNormal
              This object will have the result of the fitting procedure.
        
        solution: Once the fitting is completed the solution object is returned.
        
            Attributes of 'solution' include:
                fun: value of the function in the minimum.
                x: value of the minimization parameters in the minimum.
                nit: number of iterations of the minimization algorithm.
                nfev: number of evaluations of the function.
                success: if True it implies that the minimization fullfills all the conditions.
        
    Secondary attributes:
    
        Ndim: Number of mus, int.    
        Ncorr: Number of correlations, int.
        minparams: Array of minimization parameters at any stage in minimization process, array.
        scales: Array of scales to convert minparams to uparams (unbound) and viceversa.
        uparams: Array of unbound minimization parameters. 

    Private attributes:
        _sigmax = 10 : Maximum value of sigma parameter.
        _ignoreWarnings = True : When matrix is singular or non-positive definite, no warning will be shown.
        
    Initialization:

        Objects of this class must be always initialized with the number of gaussians and the number
        of random variables.
        
        Example: FitCMND(1,3)
        
    Basic methods:
    
        setParams: set the value of the basic params (minparams, scales, etc.)
            
            Parameters:
                mu=0.5 : Value of all initial mus.
                sigma=1.0 : Value of all initial sigmas.
                rho=0.5 : Value of all initial rhos.
                
            Returns: None
            
                It updates minparams, scales and uprams.

    Examples:

        np.random.seed(1)
        weights=[[1.0]]
        mus=[[1.0,0.5,-0.5]]
        sigmas=[[1,1.2,2.3]]
        angles=[[10*Angle.Deg,30*Angle.Deg,20*Angle.Deg]]
        Sigmas=UtilStats.calcCovarianceFromRotation(sigmas,angles)
        CMND=ComposedMultiVariateNormal(mus=mus,weights=weights,Sigmas=Sigmas)
        data=CMND.rvs(10000)
        F=FitCMND(Ngauss=CMND.Ngauss,Nvars=CMND.Nvars)
        F.cmnd._fevfreq=200
        bounds=None
        #bounds=F.setBounds(boundw=(0.1,0.9))
        #bounds=F.setBounds(boundr=(-0.9,0.9))
        #bounds=F.setBounds(bounds=(0.1,0.9*F._sigmax))
        #bounds=F.setBounds(boundsm=((-3,3),(-2,2),(-2,2)),boundw=(0.1,0.9),bounds=(0.1,0.9*F._sigmax),boundr=(-0.9,0.9))
        print(bounds)
        Util.elTime(0)
        #F.fitData(data,verbose=0,tol=1e-3,options=dict(maxiter=100,disp=True),bounds=bounds)
        F.fitData(data,verbose=0,tol=1e-3,options=dict(maxiter=100,disp=True),method=None,bounds=bounds)
        T=Util.elTime()
        print(F.cmnd)
        G=F.plotFit(figsize=3,hargs=dict(bins=30,cmap='YlGn'),sargs=dict(s=0.5,edgecolor='None',color='r'))
        F.saveFit("/tmp/fit.pkl",useprefix=False)
        F._loadFit("/tmp/fit.pkl")
        F.saveFit("/tmp/nuevo.pkl",useprefix=True,myprefix="test")

    """
    
    
    #Constants
    #Maximum value of sigma
    _sigmax=10
    _ignoreWarnings=True
    
    def __init__(self,objfile=None,Ngauss=1,Nvars=2):
        
        if objfile is not None:
            self._loadFit(objfile)
        else:
            #Basic attributes
            self.Ngauss=Ngauss
            self.Nvars=Nvars
            self.Ndim=Ngauss*Nvars
            self.Ncorr=int(Nvars*(Nvars-1)/2)

            #Define the model cmnds
            self.cmnd=ComposedMultiVariateNormal(Ngauss=Ngauss,Nvars=Nvars)

            #Set parameters
            self.setParams()
            
        #Other
        self.fig=None
        self.prefix=""
        
    def setParams(self,mu=0.5,sigma=1.0,rho=0.5):
        """
        Set the value of the basic params (minparams, scales, etc.)
        
        Parameters:
            mu=0.5 : Value of all initial mus.
            sigma=1.0 : Value of all initial sigmas.
            rho=0.5 : Value of all initial rhos.

        Returns: None

            It updates minparams, scales and uprams.
        """
        #Define the initial parameters
        #         mus             sigmas          correlations
        minparams=[mu]*self.Ndim+[sigma]*self.Ndim+[1+rho]*self.Ngauss*self.Ncorr
        scales=[0]*self.Ndim+[self._sigmax]*self.Ndim+[2]*self.Ngauss*self.Ncorr
        if self.Ngauss>1:
            self.extrap=[]
            minparams=[1/self.Ngauss]*self.Ngauss+minparams
            scales=[1]*self.Ngauss+scales
        else:
            self.extrap=[1]
            
        self.minparams=np.array(minparams)
        self.scales=np.array(scales)
        self.uparams=Util.tIF(self.minparams,self.scales,Util.f2u)
        
    def pmap(self,minparams):
        """
        Mapping routine used in sampleCMNDLikelihood.  Mapping may change depending on the 
        complexity of the parameters to be minimized.  Here we assume that all parameters in
        the stdcorr vector is susceptible to be minimized (with the exception of weights in the 
        case of Ngauss=1 when this parameter should not be included.)
        
        Parameters:
            minparams: minimization parameters.
            
        Return:
            stdcorr: flatten parameters with correlations.
        """
        stdcorr=np.array(self.extrap+list(minparams))
        stdcorr[-self.Ngauss*self.Ncorr:]-=1
        return stdcorr
    
    def logL(self,data):
        """
        Value of the -log(Likeligood)
        
        Parameters:
            data: Array with data, array, Nsam x Nvars
            
        Return: 
            logL: value of the -log(Likelihood)
        """
        
        logL=self.cmnd.sampleCMNDLikelihood(self.uparams,
                                            data=data,
                                            pmap=self.pmap,
                                            tset="stdcorr",
                                            scales=self.scales)
        return logL
    
    def fitData(self,data,verbose=0,advance=0,**args):
        """
        Minimization procedure
        
        Parameters:
            data: Array with data, array, Nsam x Nvars
            verbose=0 : verbosity level for the sampleCMNDLikelihood routine
            advance=0 : If larger than 0 show advance each "advance" iterations.
            **args: Options of the minimize routine (eg. tol=1e-6)
                    A particularly interesting parameter is the minimization method:
                        Available methods: 
                            Slow but sure
                                Powell
                            Fast but unsure:
                                CG, BFGS, COBYLA, SLSQP
                                
        Return: None
        
            It updates the solution attribute
            
        Examples:
            F=FitCMND(1,3)
            F.fitData(data,verbose=0,tol=1e-3,options=dict(maxiter=100,disp=True))

        """
        if advance:
            advance=int(advance)
            self.neval=0
            def _advance(X,show=False):
                if self.neval==0:
                    print(f"Iterations:")
                if self.neval%advance==0 or show:
                    vars = np.array2string(X, separator=', ', precision=4, max_line_width=np.inf, formatter={'float_kind':lambda x: f"{x:.2g}"})
                    fun = self.cmnd.sampleCMNDLikelihood(X,data,self.pmap,"stdcorr",self.scales,verbose)
                    print(f"Iter {self.neval}:\n\tVars: {vars}\n\tLogL/N: {fun/len(data)}")
                self.neval+=1
        else:
            _advance = None

        self.data=np.copy(data)
        self.cmnd._ignoreWarnings=self._ignoreWarnings
        self.minargs=dict(method="Powell")
        self.minargs.update(args)
        self.solution=minimize(self.cmnd.sampleCMNDLikelihood,
                               self.minparams,
                               callback=_advance,
                               args=(data,self.pmap,"stdcorr",self.scales,verbose),
                               **self.minargs)
        if advance:
            _advance(self.solution.x,show=True)
        self.uparams=self.solution.x

        #Set the new params
        self._invParams(self.cmnd.stdcorr)
        self._updatePrefix()
            
    def _loadFit(self,objfile):
        F=pickle.load(open(objfile,"rb"))
        for k in F.__dict__.keys():
            setattr(self,k,getattr(F,k))
        self._updatePrefix()
    
    def plotFit(self,N=10000,figsize=2,props=None,ranges=None,hargs=dict(),sargs=dict()):
        """
        Plot the result of the fitting procedure
        
        Parameters:
        
            N = 10000: number of points used to build a representation of the marginal distributions.
            props = None: Array with the name of the properties.  Ex. props=["x","y"]
            ranges = None: Array of ranges of the properties. Ex. ranges=[-3,3],[-5,5]
            figsize = 2: size of each axis. 
            hargs = None: Dictionary with options for the hist2d function.  Ex. hargs=dict(bins=50)
            sargs = dict(): Dictionary with options for the scatter plot.  Ex. sargs=dict(color='r')
            
        Examples:
            F=FitCMND(1,3)
            F.fitData(data,verbose=0,tol=1e-3,options=dict(maxiter=100,disp=True))
            G=F.plotFit(figsize=3,hargs=dict(bins=30,cmap='YlGn'),sargs=dict(s=0.5,edgecolor='None',color='r'))
            
        """
        Xfits=self.cmnd.rvs(N)
        properties=dict()
        for i in range(self.Nvars):
            symbol=string.ascii_letters[i] if props is None else props[i]
            if ranges is not None:
                rang=ranges[i]
            else:
                rang=None
            properties[symbol]=dict(label=f"${symbol}$",range=rang)
            properties[symbol]=dict(label=f"${symbol}$",range=rang)
            
        if self.Nvars>2:
            G=PlotGrid(properties,figsize=figsize)
            G.plotHist(Xfits,**hargs)
            G.scatterPlot(self.data,**sargs);
            G.fig.tight_layout()
            self.fig=G.fig
            return G
        else:
            keys=list(properties.keys())
            fig=plt.figure(figsize=(5,5))
            ax=fig.gca()
            ax.hist2d(Xfits[:,0],Xfits[:,1],**hargs)
            ax.scatter(self.data[:,0],self.data[:,1],**sargs)
            ax.grid()
            ax.set_xlabel(properties[keys[0]]["label"])
            ax.set_ylabel(properties[keys[1]]["label"])
            fig.tight_layout()
            self.fig=fig
            return fig
        
    def _invParams(self,stdcorr):
        minparams=np.copy(stdcorr)
        minparams[-self.Ngauss*self.Ncorr:]+=1
        self.minparams=minparams[1:] if self.Ngauss==1 else minparams

    def _updatePrefix(self,myprefix=None):
        """
        Update prefix of fit.
        
        Prefix has two parts: the number of gaussians used and a hash computed from the object.
        
        Preffix change if:
            - Data change.
            - Initial minimization parameters change (e.g. if the fit is ran twice)
            - Minimization parameters are changed.
            - Bounds are changed.

        Alternative prefix:
            self.hash=md5(pickle.dumps([self.Ngauss,self.data])).hexdigest()[:5]
            self.hash=md5(pickle.dumps(self.__dict__)).hexdigest()[:5]
            self.hash=md5(pickle.dumps(self.minparams)).hexdigest()[:5]
            self.hash=md5(pickle.dumps(self.cmnd)).hexdigest()[:5]
        """
        self.hash=md5(pickle.dumps(self)).hexdigest()[:5]
        if myprefix is not None:
            myprefix=f"_{myprefix}"
        self.prefix=f"{self.Ngauss}cmnd{myprefix}_{self.hash}"
        
    def saveFit(self,objfile=None,useprefix=True,myprefix=None):
        """
        Pickle the result of a fit
        
        Parameteres:
            objfile=None: name of the file where the fit will be stored, string.  If none the name is set 
                          by the routine as FitCMND.pkl
                          
            useprefix=True: use a prefix in the filename of the pickle file.  The prefix is normally 
                            {Ngauss}cmnd_{hash}
                            
                            Example: If objfile="fit.pkl" the final filename will be fit-1mnd_asa33.pkl
        """
        self.fig=None
        self._updatePrefix(myprefix)
        if objfile is None:
            objfile=f"/tmp/FitCMND.pkl"
        if useprefix:
            parts=os.path.splitext(objfile)
            objfile=f"{parts[0]}-{self.prefix}{parts[1]}"
        pickle.dump(self,open(objfile,"wb"))
        
    def setBounds(self,boundw=None,bounds=None,boundr=None,boundsm=None):
        """
        Set the minimization parameters
        
        Parameters:
            boundw=None: bound of weights, tuple.  If None, boundw = (-np.inf,np.inf)
            bounds=None: bound of weights, tuple.  If None, boundw = (-np.inf,np.inf)
            boundr=None: bound of weights, tuple.  If None, boundw = (-np.inf,np.inf)
            
            boundsm=None: bounds of averages, tuple of tuples (Nvars).  
                          If None, boundsm = (-np.inf,np.inf)
                          
                          Normally the bounds on averages must be expressed in this way:
                          
                              boundsm=((-min_1,max_1),(-min_2,max_2),...)
                              
                          Example: for Nvars = 2:
                          
                              boundsm=((-2,1),(-3,0))

            bounds: a list or tuple with the bounds of weights, mus, sigmas and rhos of each variable.
                Examples: 
        """
        if boundsm is None:
            boundsm=((-np.inf,np.inf),)*self.Nvars

        # Regular bounds
        if boundw is None:
            boundw=(-np.inf,np.inf)
        else:
            boundw=tuple([Util.f2u(bw,1) for bw in boundw])

        if bounds is None:
            bounds=(-np.inf,np.inf)
        else:
            bounds=tuple([Util.f2u(bs,self._sigmax) for bs in bounds])

        if boundr is None:
            boundr=(-np.inf,np.inf)
        else:
            boundr=tuple([Util.f2u(1+br,2) for br in boundr])

        bounds=(*((boundw,)*self.Ngauss),
                *(boundsm*self.Ngauss),
                *((bounds,)*self.Nvars*self.Ngauss),
                *((boundr,)*self.Ngauss*int(self.Nvars*(self.Nvars-1)/2)))
        self.bounds=bounds
        
        if self.Ngauss==1:
            bounds=bounds[1:]
        return bounds

