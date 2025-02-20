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
        logMmin = self.param_dict["logMmin"]+np.log10(1+fenv*erf((delta-deltaenv)/sigmaenv))# Caculate Mmin after considering environment bias.
        mean_ncen = 0.5 * f_Gamma * (1.0 + erf((logM - logMmin) / sigma_logM))# N_cen
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
        mean_nsat *= mean_ncen# N_sat
        return mean_nsat
