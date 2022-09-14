*** Minimize violation of fluxes that are zero due to non-production of enzymes ***
*       Author: Hoang Dinh
***********************************************************************************

$INLINECOM /*  */
$include "./min_flux_violation_GAMS_settings.txt"
$setGlobal nscale 1000

options
	LP = soplex /*Solver selection*/
	limrow = 1000000 /*number of equations listed, 0 is suppresed*/
	limcol = 1000000 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 1000000 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets
i
$include "%species_path%"
j
$include "%rxns_path%"
rxns_enzsyn(j)
$include "%rxns_enzsyn_path%"
rxns_enzload(j)
$include "%rxns_enzload_path%"
nuc_translation(j)
$include "%nuc_trans_path%"
mito_translation(j)
$include "%mito_trans_path%"
uptake(j) /*list of uptake so that all of them are properly turned off*/
$include "%uptake_path%"
media(j) /*list of allowable uptake based on simulated media conditions*/
$include "%media_path%"
rxns_inactive(j)
$include "%rxns_inactive_path%"
prodata_set(j)
$include "%proteome_data_set_path%"
;

Parameters
S(i,j)
$include "%sij_path%"
NAA(j)
$include "%prolen_path%"
pro_val(j)
$include "%proteome_data_path%"
;

Variables
z, v(j), venzSlack
;
venzSlack.lo = 0; venzSlack.up = %venzSlackAllow%;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0; v.up(j) = 1e4 * %nscale%;

* Disable enzyme synthesis and enzyme load network
v.fx(j)$rxns_enzsyn(j) = 0;
v.fx(j)$rxns_enzload(j) = 0;

* Media
v.up(j)$uptake(j) = 0;
v.up(j)$media(j) = 1e4 * %nscale%;

* Turning off all versions of biomass dilution reaction
* You need to turn on the respective version corresponding to your growth condition
v.fx('BIOSYN-COFACTORANAEROBIC') = 0; v.fx('BIOSYN-compCERANAEROBIC') = 0;
v.fx('BIOSYN-BIODILNOGAM') = 0; v.fx('BIOSYN-BIODILAERO') = 0; v.fx('BIOSYN-BIODILBATCHANAERO') = 0;
v.fx('BIOSYN-BIODILCHEMOANAERO') = 0; v.fx('BIOSYN-BIODILSTARVE') = 0;

* Growth rate, substrate and oxygenation, and secretions
$include "%phenotype_path%"

* Carbohydrate, RNA, and protein fraction constraints
* C-lim protein fraction
*v.fx('BIOSYN-PROTTOBIO') = %mu% * (36.94 + 34.22*%mu%) / 100 * %nscale%;
* N-lim protein fraction
*v.fx('BIOSYN-PROTTOBIO') = %mu% * (10.57 + 108.56*%mu%) / 100 * %nscale%;
* RNA fraction for both C-lim and N-lim
*v.fx('BIOSYN-RNATOBIO') = %mu% * (3.79 + 13.88*%mu%) / 100 * %nscale%;

* Proteome allocation for purposes other than metabolism and ribosome
* Clim
*v.up('BIOSYN-PROTMODELED') = (1 - %nonmodeled_proteome_allocation%) * %mu% * (36.94 + 34.22*%mu%) / 100 * %nscale%;
* Nlim
*v.up('BIOSYN-PROTMODELED') = (1 - %nonmodeled_proteome_allocation%) * %mu% * (10.57 + 108.56*%mu%) / 100 * %nscale%;

* Additional constraints
*v.fx('RXN-EX_pyr_e_FWD-SPONT') = 0;
*v.fx('RXN-EX_btd_e_FWD-SPONT') = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, ProData, Nonmodel
;

Obj..				z =e= sum(j$rxns_inactive(j), v(j));
*Obj..				z =e= venzSlack;
Stoic(i)..			sum(j, S(i,j)*v(j)) =e= 0;
RiboCapacityMito.. 		v('RIBOSYN-ribomito') * %kribomito% =g= %mu% * sum(j$mito_translation(j), NAA(j) * v(j));
RiboCapacityNuc.. 		v('RIBOSYN-ribonuc') * %kribonuc% =g= %mu% * sum(j$nuc_translation(j), NAA(j) * v(j));
ProData(j)$prodata_set(j)..	v(j) =e= pro_val(j) * (1 - venzSlack);
Nonmodel..			v('BIOSYN-PROTMODELED') =l= (1 - %nonmodeled_proteome_allocation%) * v('BIOSYN-PROTTOBIO');

*** BUILD OPTIMIZATION MODEL ***
Model rba
/Obj, Stoic, RiboCapacityNuc, RiboCapacityMito, ProData, Nonmodel
/;
rba.optfile = 1;

*** SOLVE ***
Solve rba using lp minimizing z;

file ff /min_flux_violation.modelStat.txt/;
put ff;
put rba.modelStat/;
putclose ff;

file ff2 /min_flux_violation.flux_gamsscaled.txt/;
put ff2;
loop(j,
	if ( (v.l(j) gt 1e-12),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
	);
);
putclose ff2;

file ff3 /min_flux_violation.flux_essential_inactive_rxns_gamsscaled.txt/;
put ff3;
loop(j$rxns_inactive(j),
	if ( (v.l(j) gt 1e-12),
		put j.tl:0, system.tab, 'v', system.tab, v.l(j):0:15/;
	);
);
putclose ff3;

file ff4 /min_flux_violation.venzSlack.txt/;
put ff4;
put venzSlack.l:0:11/;
putclose ff4;
