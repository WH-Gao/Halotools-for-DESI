# Generate galaxy mock with Halotools

## 1.Read halo files
To generate galaxy mock using Halotools, we need read halo files first and save as hdf5 files.

In read_halo.ipynb, I provide two examples to read Aemulus_nu simulation and UNIT simulation. 

The basical information of halos contain: 'halo_id', 'halo_m200m'/'halo_mvir', 'halo_r200m'/'halo_rvir', 'halo_x', 'halo_y', 'halo_z', 'halo_vx', 'halo_vy', 'halo_vz', 'halo_vrms', 'halo_vmax', 'halo_rs', 'halo_upid'. To be noticed, for halotools, if there is no 'halo_mvir' and 'halo_rvir' being read, there may be some error. So for Aemulus_nu I just duplicate 'halo_m200m' and 'halo_r200m' as 'halo_mvir' and 'halo_rvir'. 

When reading halo files, setting row_cut_eq_dict = {'halo_upid': -1} can remove subhalos, and Set row_cut_min_dict = {'halo_mvir': min_halo_mass} can remove halos under some mass.

In our model, to add galaxy bias related to environment, I also add 'halo_delta' in halo files. The halo files path of Aemulus_nu and UNIT are '/home/wenhao/Halotools/halo/out_xxx/C_xxx_halo_vir.list' and '/home/wenhao/UNIT_sim/xxx/xxx_halo_shape.list'.

## 2.Generate galaxy mock
Our fiducial HOD model is from Aemulus V.

$N_{cen}(M)=\frac{1}{2} f_{max}[1+erf(\frac{log_{10}M-log_{10}M_{min}}{\sigma_{logM}})]$

$N_{sat}(M)=(\frac{M}{M_{sat}})^{\alpha}exp(-\frac{M_{cut}}{M})N_{cen}(M)$

Also we introduce additional parameters $\eta_{con}$, $\eta_{vs}$ and $\eta_{vc}$ to rescale galaxy concentration, satellites velocity and centrals velocity respectively. If we use Biased NFW model providen by Halotools, there may be error when $\eta_{con}$<1, hence I just multiply halo concentration with eta_con. 

We also have $\gamma_f$ to measure the deviation with GR. There galaxy bias parameters are $f_{env}$, $\sigma_{env}$ and $\delta_{env}$:

$\bar{M_{min}}=M_{min}[1+f_{env}erf(\frac{\delta - \delta_{env}}{\sigma_{env}}]$

The way to generate galaxy mock is in halotools.ipynb. For Aemulus_nu simulations, we have to calculate halo_concentration and halo_vrms by ourselves. 

When generating galaxy mock, we get galaxy positions first and then adjust their velocities by hand.

For model with non-poissonity, we can just add parameter Gp as in halotools_Gp.ipynb.

After generating galaxy mock, we can get a galaxy catalog  and we can obtaion galaxy position and velocity through: model_instance.mock.galaxy_table['x'],model_instance.mock.galaxy_table['y'],model_instance.mock.galaxy_table['z'],model_instance.mock.galaxy_table['vx'],model_instance.mock.galaxy_table['vy'],model_instance.mock.galaxy_table['vz'].
