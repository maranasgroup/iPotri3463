*GAMS code for the knockdown analysis for poplar drought model

$INLINECOM /*  */
$ONLISTING

options 

	limrow = 1000
	optCR = 1E-9
	optCA = 1E-9
	iterlim = 1000000
	decimals = 8
	reslim = 1000000
  lp = cplex

;
        

SETS
        j   Reactions of k-ath model
$include "rxns_H3.txt"

        i   Metabolites of k-ath model
$include "mets_H3.txt"

				j_kd(j)		reaction knockdown candidates
$include "kd_candidates_H3.txt"

;

*add an alias needed for FVA
alias (j_kd1, j_kd);

PARAMETERS

  S(i,j) 				Stoichiometric matrix for ED 
$include "Sij_H3.txt"

	LB(j)					lower bounds of reactions
$include "LB_H3.txt"

	UB(j)					lower bounds of reactions
$include "UB_H3.txt"

	imp_LB(j_kd)	upper bound to impose on the candidate
$include "impose_LBs_H3.txt"

	imp_UB(j_kd)	upper bound to impose on the candidate
$include "impose_UBs_H3.txt"

	c(j)					objective vector, used to update objective function

	min(j)				stores objective values for min FVA

	max(j)				stores objective values for max FVA

;

VARIABLES

	z				    	primal objective function (rxn fluxes in each loop) for ED

  v(j)       		Flux of ED rxns

;
	

v.lo(j)=LB(j);
v.up(j)=UB(j);

*growth rates for the particular harvest
SCALAR

	leaf_g_H3		/1.1840431049994806e-05/
	stem_g_H3		/1.184043105e-05/
	root_g_H3		/0.0021022624299999997/
	beta
	stem_frac	
	root_frac
	leaf_frac
	kd_level		/0.5/

;

beta = root_g_H3 / stem_g_H3

*********************EQUATIONS NAMES**************************
EQUATIONS

	/*standard FVA constraints*/
	obj					primal objective function (rxn fluxes in each loop) for ED
	mb(i) 	 		Mass balance for ED

	/*model specific constraints*/
	eqn14				enforces that at most co2 fixed from stem in leaf is 2.7% of that uptaken
	eqn15				enforces that equal amounts of CO2 are passed to stem and exported from root
	eqn15_1			enforces that equal amounts of CO2 are passed to stem and exported from root
	eqn17				enforces that 10% of night respired CO2 is stored for the day
	eqn17_1			requires that CO2 stored at night is used in the day
	eqn18				enforces equal starch storage and use between day and night in leaf
	eqn19				enforces equal starch storage and use between day and night in leaf
	eqn21				requires photorespiration

	/*restrictions on growth rates sicne we are seeking for change in yield*/
	LeqS				leaf growth rate equals stem
	RpropS			ensures root growth rate is in same proportion to stem

;

****************************** Equations*******************************
obj..		  					z =e= v('RXN_BiomassRxn__S___DARK');
mb(i)..     				sum(j,S(i,j)*v(j)) =e= 0;
eqn14..							v('RXN_L2CP2_CO2_LIGHT') + 0.027 * v('RXN_E2L_CO2_LIGHT') =g= 0;
eqn15..							v('RXN_E2R_CO2_LIGHT') - v('RXN_R2CP1_CO2_LIGHT') =e= 0;
eqn15_1..						v('RXN_E2R_CO2_DARK') - v('RXN_R2CP1_CO2_DARK') =e= 0;
eqn17..							0.1 * v('RXN_E2S_CO2_DARK') - v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn17_1..						v('SK_MET_carbon_dioxide_c__S___LIGHT') + v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn18..							v('DM_MET_starch_p__L___LIGHT') + v('SK_MET_starch_p__L___DARK') =e= 0;
eqn19..							v('DM_MET_starch_p__S___LIGHT') + v('SK_MET_starch_p__S___DARK') =e= 0;
eqn21..							v('RXN_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p__L___LIGHT') + 2 * v('RXN_RXN_961_p__L___LIGHT') =e= 0;
LeqS..							v('RXN_BiomassRxn__L___DARK') - v('RXN_BiomassRxn__S___DARK') =e= 0;
RpropS..						v('RXN_BiomassRxn__S___DARK') * beta - v('RXN_BiomassRxn__R___DARK') =e= 0;

*****************************************************************************

*note: all the puclose and append stuff allows for real-time writing
file RESULT /knockdown_result_H3_05.txt/;
put RESULT;
put "candidate",system.tab,"imposed_LB",system.tab,"imposed_UB",system.tab,"stem_growth",system.tab,"leaf_growth",system.tab,"root_growth",system.tab,"stem_fraction",system.tab,"leaf_fraction",system.tab,"root_fraction",system.tab,"model_stat",system.tab,"solver_stat";

model KNOCKDOWN
/
	obj
	mb
	eqn14
	eqn15
	eqn15_1
	eqn17
	eqn17_1
	eqn18
	eqn19
	eqn21
	LeqS
	RpropS
/
;

PUTCLOSE;

LOOP(j_kd1,
	
	result.ap = 1;
	PUT result;

	/*reset bounds*/
	v.lo(j)=LB(j);
	v.up(j)=UB(j);

	/*set bounds of limited reaction*/

	if(imp_LB(j_kd1)<=imp_UB(j_kd1),

		/*only solve for those where the bounds make sense*/
		v.lo(j_kd1) = kd_level * imp_LB(j_kd1);
		v.up(j_kd1) = kd_level * imp_UB(j_kd1);

		/*free up biomasss reactions to allow for fractional growth*/
		v.lo('RXN_BiomassRxn__R___DARK') = 0;
		v.up('RXN_BiomassRxn__R___DARK') = 10000;
		v.lo('RXN_BiomassRxn__L___DARK') = 0;
		v.up('RXN_BiomassRxn__L___DARK') = 10000;
		v.lo('RXN_BiomassRxn__S___DARK') = 0;
		v.up('RXN_BiomassRxn__S___DARK') = 10000;
	
		KNOCKDOWN.optfile=1;
		Solve KNOCKDOWN using lp maximizing z;

		stem_frac = v.l('RXN_BiomassRxn__S___DARK') / stem_g_H3;
		leaf_frac = v.l('RXN_BiomassRxn__L___DARK') / leaf_g_H3;
		root_frac = v.l('RXN_BiomassRxn__R___DARK') / root_g_H3;
	
		put j_kd1.tl:0:65,system.tab,imp_LB(j_kd1):0:8,system.tab,imp_LB(j_kd1):0:8,system.tab,v.l('RXN_BiomassRxn__S___DARK'):0:8,system.tab,v.l('RXN_BiomassRxn__L___DARK'):0:8,system.tab,v.l('RXN_BiomassRxn__R___DARK'):0:8,system.tab,stem_frac,system.tab,leaf_frac,system.tab,root_frac,system.tab,KNOCKDOWN.Modelstat,system.tab,KNOCKDOWN.Solvestat;

	else

		put j_kd1.tl:0:65,system.tab,imp_LB(j_kd1):0:8,system.tab,imp_LB(j_kd1):0:8,system.tab,0,system.tab,0,system.tab,0,system.tab,0,system.tab,0,system.tab,0,system.tab,"problem",system.tab,"problem";	

	);
	

	PUTCLOSE;

);
