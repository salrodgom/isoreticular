# PLUMED reweighting template for re-weighting an OPES Expanded multithermal-multibaric
# simulation to a TARGET (T, P).
# Placeholders (filled by run_reweight_G<k>.sh):
#   __IDX__      walker index (0, 1, 2)
#   __ENE_REF__  reference energy in kJ/mol (for numerical-stability shift)
#   __VOL_REF__  reference volume in nm^3 (idem)
#   __P_SIM__    simulation pressure (the OPES centre, in 0.06022140857*P_bar units)
#   __P_TARGET__ target pressure to reweight to (same units)
#   __TEMP__     target temperature in K
energy: READ FILE=COLVARS.__IDX__ VALUES=energy IGNORE_TIME
volume: READ FILE=COLVARS.__IDX__ VALUES=vol IGNORE_TIME
cell:   READ FILE=COLVARS.__IDX__ VALUES=cell IGNORE_TIME
delta:  READ FILE=COLVARS.__IDX__ VALUES=delta IGNORE_TIME
bias:   READ FILE=COLVARS.__IDX__ VALUES=bias IGNORE_TIME

# Shift energy and volume (numerical stability)
rvol: COMBINE ARG=volume PARAMETERS=__VOL_REF__ PERIODIC=NO
rene: COMBINE ARG=energy PARAMETERS=__ENE_REF__ PERIODIC=NO

# Pure-pressure reweighting (T_target == T_sim).
# PLUMED 2.9+ does not accept ENERGY together with REWEIGHT_PRESSURE in the
# absence of REWEIGHT_TEMP: it interprets the input as ambiguous. Pure-P
# reweighting only requires VOLUME and the (P_sim, P_target) pair.
# w(P -> P') = exp(-beta (P' - P) V)
bias_weights: REWEIGHT_BIAS TEMP=__TEMP__ ARG=bias
temp_press_weights: REWEIGHT_TEMP_PRESS TEMP=__TEMP__ \
                     PRESSURE=__P_SIM__ REWEIGHT_PRESSURE=__P_TARGET__ \
                     VOLUME=rvol

avg_delta: AVERAGE ARG=delta LOGWEIGHTS=bias_weights,temp_press_weights
avg_cell:  AVERAGE ARG=cell  LOGWEIGHTS=bias_weights,temp_press_weights

PRINT ARG=cell,delta,energy,vol,bias_weights,temp_press_weights \
      FILE=COLVAR_REWEIGHT.__IDX__ STRIDE=1
