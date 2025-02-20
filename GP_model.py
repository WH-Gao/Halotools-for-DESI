from halotools.empirical_models import model_defaults
from halotools.empirical_models.occupation_models import OccupationComponent
import math as mt
import random as rd
import scipy as sp
import numpy as np
from scipy.special import erf

def match(a,b):
    sorted_indices = np.argsort(a)
    sorted_halo = a[sorted_indices]

    positions_in_sorted = np.searchsorted(sorted_halo, b)
    positions = sorted_indices[positions_in_sorted]
    
    return positions

def gpd(a,b,n):
    p = np.ones(len(a))
    if np.max(a)>0:
        p[np.where(a>0)] = a[np.where(a>0)]*(a[np.where(a>0)]+b*n[np.where(a>0)])**(n[np.where(a>0)]-1)/sp.special.factorial(n[np.where(a>0)],exact=True)*np.exp(-a[np.where(a>0)]-b*n[np.where(a>0)])
    return p

class Zcens(OccupationComponent):
    def __init__(
        self,
        threshold=model_defaults.default_luminosity_threshold,
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        **kwargs
    ):
        upper_occupation_bound = 1.0
        self.sec_haloprop_key = 'halo_delta'
        super(Zcens, self).__init__(
            gal_type="centrals",
            threshold=threshold,
            upper_occupation_bound=upper_occupation_bound,
            prim_haloprop_key=prim_haloprop_key,
            **kwargs
        )
    def mean_occupation(self,**kwargs):
        mass = kwargs["table"][self.prim_haloprop_key]
        delta = kwargs["table"][self.sec_haloprop_key]
        logM = np.log10(mass)
        sigma_logM = self.param_dict["sigma_logM"]
        f_Gamma = self.param_dict["f_Gamma"]
        fenv = self.param_dict["fenv"]
        deltaenv = self.param_dict["deltaenv"]
        sigmaenv = self.param_dict["sigmaenv"]
        logMmin = self.param_dict["logMmin"]+np.log10(1+fenv*erf((delta-deltaenv)/sigmaenv))
        mean_ncen = 0.5 * f_Gamma * (1.0 + erf((logM - logMmin) / sigma_logM))
        return mean_ncen

class ZSats(OccupationComponent):
    def __init__(
        self,
        threshold=model_defaults.default_luminosity_threshold,
        prim_haloprop_key=model_defaults.prim_haloprop_key,
        modulate_with_cenocc=True,
        cenocc_model=None,
        **kwargs
    ):
        upper_occupation_bound = float("inf")
        super(ZSats, self).__init__(
            gal_type="satellites",
            threshold=threshold,
            upper_occupation_bound=upper_occupation_bound,
            prim_haloprop_key=prim_haloprop_key,
            **kwargs
        )
        self.modulate_with_cenocc = modulate_with_cenocc
        self.central_occupation_model = cenocc_model
    def mean_occupation(self,**kwargs):
        if self.modulate_with_cenocc:
            for key, value in self.param_dict.items():
                if key in self.central_occupation_model.param_dict:
                    self.central_occupation_model.param_dict[key] = value
        mass = kwargs["table"][self.prim_haloprop_key]
        Mcut = 10**self.param_dict["mcut"]
        Msat = 10**self.param_dict["msat"]
        Alphasat = self.param_dict["alphasat"]
        mean_nsat = np.zeros_like(mass)
        mean_nsat = np.exp(-1 * Mcut / mass) * (mass / Msat) ** Alphasat
        mean_ncen = self.central_occupation_model.mean_occupation(**kwargs)
        mean_nsat *= mean_ncen
        return mean_nsat
    def mc_occupation(self, **kwargs):
        b = self.param_dict["Gp"]
        if b>=0:
            table = kwargs['table']
            meanocc = self.mean_occupation(**kwargs)
            a = meanocc*(1-b)
            p = np.random.rand(len(a))
            x = np.zeros(len(a))
            temp = gpd(a,b,x)
            while np.min(temp-p)<0:
                if True in np.isinf(temp):
                    temp[np.where(np.isinf(temp))] = 2
                elif np.max(x[np.where(temp-p<0)])<170: # Make sure N_sat<170 to avoid numeric overflow.
                    x[np.where(temp-p<0)]+=1
                    temp[np.where(temp-p<0)] += gpd(a[np.where(temp-p<0)],b,x[np.where(temp-p<0)])
                else:
                    x[np.where((temp-p<0))]+=1
                    temp[np.where((temp-p<0)&(x<170))] += gpd(a[np.where((temp-p<0)&(x<170))],b,x[np.where((temp-p<0)&(x<170))])
                    temp[np.where((temp-p<0)&(x>=170))] = 2
            x=np.ceil(x).astype(int)
            table['halo_num_satellites'] = x
            return x
        else:
            table = kwargs['table']
            meanocc = self.mean_occupation(**kwargs)
            a = meanocc*(1-b)
            maxn=a/(-b)
            norm=np.zeros(len(a))
            n=np.zeros(len(a))
            for i in range(int(min(max(a/(-b)),170))):
                norm[np.where(n<maxn)]+=gpd(a[np.where(n<maxn)],b,n[np.where(n<maxn)])
                n+=1
            p = np.random.rand(len(a))*norm
            x = np.zeros(len(a))
            temp = gpd(a,b,x)
            pn = np.random.rand(len(a))
            x[np.where((maxn<2)&(pn<meanocc))]+=1 # If m<2, turn to Bernoulli distribution.
            temp[np.where((maxn<2))]+=10
            while np.min(temp-p)<0:
                if True in np.isinf(temp):
                    temp[np.where(np.isinf(temp))] = 2
                elif np.max(x[np.where(temp-p<0)])<170:
                    x[np.where(temp-p<0)]+=1
                    temp[np.where(temp-p<0)] += gpd(a[np.where(temp-p<0)],b,x[np.where(temp-p<0)])
                else:
                    x[np.where((temp-p<0))]+=1
                    temp[np.where((temp-p<0)&(x<170))] += gpd(a[np.where((temp-p<0)&(x<170))],b,x[np.where((temp-p<0)&(x<170))])
                    temp[np.where((temp-p<0)&(x>=170))] = 2
            x=np.ceil(x).astype(int)
            table['halo_num_satellites'] = x
            return x
