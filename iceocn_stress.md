
Ocn-ice stress is calculated in `models/ice/cice/src/source/ice_dyn_evp.F90` function `evp_prep2`.
Stress is named `strocnx`, `strocny`. It is calculated using `waterx`, `watery` which are rotated
POP grid ocean current speed `uocn`, `vocn`.

Also note that ice fraction has been multiplied in the calculation (air-ocn calculation does not 
take ice-fraction into account).


In a nutshell, code reads

`stress ~ (ice fraction) * (relative speed of ocn and ice)^2 * (ocn speed)`.

So sea-ice is treated as static when calculating its work on ocn.
