******** Approximate enzyme concentration from proteome *********
*       Author: Hoang Dinh
*****************************************************************

$INLINECOM /*  */
$include "./enz_from_proteome_GAMS_settings.txt"

options
	LP = soplex /*Solver selection*/
	limrow = 0 /*number of equations listed, 0 is suppresed*/
	limcol = 0 /*number of variables listed, 0 is suppresed*/
	iterlim = 1000000 /*iteration limit of solver, for LP it is number of simplex pivots*/
	decimals = 8 /*decimal places for display statement*/
	reslim = 1000000 /*wall-clock time limit for solver in seconds*/
	sysout = on /*solver status file report option*/
	solprint = on /*solution printing option*/
        
Sets
i
$include "%pro_and_enz_path%"
j
$include "%rxns_pro_and_enz_path%"
prodata(j)
$include "%rxns_pro_data_path%"
pronodata(j)
$include "%rxns_pro_nodata_path%"
enzout(j)
$include "%rxns_enz_path%"
;

Parameters
S(i,j)
$include "%sij_pro_and_enz_path%"
vprotexpmt(j)
$include "%proteome_data_path%"
;

Variables
z, v(j)
;

*** SET FLUX LOWER AND UPPER BOUNDS ***
v.lo(j) = 0;
v.up(j)$enzout(j) = 1e4;
v.up(j)$prodata(j) = vprotexpmt(j)*1e6;
v.up(j)$pronodata(j) = 0;

*** EQUATION DEFINITIONS ***
Equations
Obj, Stoic
;

Obj..			z =e= v('%enzobj%');
*Obj..			z =e= v('ENZSYN-ATPSCPLX');
*Obj..			z =e= v('ENZSYN-SNZ3SNO1');
*Obj..			z =e= v('ENZSYN-YFR053C');
Stoic(i)..		sum(j, S(i,j)*v(j)) =e= 0;

*** BUILD OPTIMIZATION MODEL ***
Model optmodel
/Obj, Stoic/;
optmodel.optfile = 1;

*** SOLVE ***
Solve optmodel using lp maximizing z;

file ff /enz_from_proteome.modelStat.txt/;
put ff;
put optmodel.modelStat/;
putclose ff;

file ff2 /enz_from_proteome.objval.txt/;
put ff2;
put z.l:0:11;
putclose ff2;
