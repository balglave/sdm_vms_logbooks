# Toy simulation-estimation examples for species distribution models fitted to 'VMS x logbooks' data

- **1_ipp_discrete_toy_model** : discrete version of a species distribution model accounting for preferential sampling (Alglave et al., 2022)

- **2_ipp_st_toy_model** : spatio-temporal version of the species distribution model accounting for preferential sampling.
(n.b. the SPDE approach is included in this toy example. The random effects are estimated on a triangulated mesh and are linked to data points by linear interpolation. This allows to gain speed.)