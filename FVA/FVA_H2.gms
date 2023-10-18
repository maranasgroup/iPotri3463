*GAMS code for Flux Variability Analysis for poplar drought model

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
$include "rxns_H2.txt"

        i   Metabolites of k-ath model
$include "mets_H2.txt"

;

*add an alias needed for FVA
alias (j1, j);

PARAMETERS

  S(i,j) 		Stoichiometric matrix for ED 
$include "Sij_H2.txt"

	LB(j)			lower bounds of reactions
$include "LB_H2.txt"

	UB(j)			lower bounds of reactions
$include "UB_H2.txt"

	c(j)			objective vector, used to update objective function

	min(j)		stores objective values for min FVA

	max(j)		stores objective values for max FVA

;

VARIABLES

	z				    	primal objective function (rxn fluxes in each loop) for ED

  v(j)       		Flux of ED rxns

;
	

v.lo(j)=LB(j);
v.up(j)=UB(j);


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

;

****************************** Equations*******************************
obj..		  					z =e= sum(j, c(j) * v(j));
mb(i)..     				sum(j,S(i,j)*v(j)) =e= 0;
eqn14..							v('RXN_L2CP2_CO2_LIGHT') + 0.027 * v('RXN_E2L_CO2_LIGHT') =g= 0;
eqn15..							v('RXN_E2R_CO2_LIGHT') - v('RXN_R2CP1_CO2_LIGHT') =e= 0;
eqn15_1..						v('RXN_E2R_CO2_DARK') - v('RXN_R2CP1_CO2_DARK') =e= 0;
eqn17..							0.1 * v('RXN_E2S_CO2_DARK') - v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn17_1..						v('SK_MET_carbon_dioxide_c__S___LIGHT') + v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn18..							v('DM_MET_starch_p__L___LIGHT') + v('SK_MET_starch_p__L___DARK') =e= 0;
eqn19..							v('DM_MET_starch_p__S___LIGHT') + v('SK_MET_starch_p__S___DARK') =e= 0;
eqn21..							v('RXN_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p__L___LIGHT') + 2 * v('RXN_RXN_961_p__L___LIGHT') =e= 0;

*****************************************************************************

*note: all the puclose and append stuff allows for real-time writing
file RESULT /FVA_result_H2.txt/;
put RESULT;
put "rxn",system.tab,"LB",system.tab,"min",system.tab,"max",system.tab,"UB",system.tab,"ModelStat",system.tab,"SolveStat";

model FVA
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
/
;

PUTCLOSE;

LOOP(j1,
	
	result.ap = 1;
	PUT result;
	
	C(j) = 0;
	C(j1) = 1;
	
	FVA.optfile=1;
	Solve FVA using lp maximizing z;
	max(j1) = z.l;

	Solve FVA using lp minimizing z;
	min(j1) = z.l;
	
	put j1.tl:0:65,system.tab,v.lo(j1):0:8,system.tab,min(j1):0:8,system.tab,max(j1):0:8,system.tab,v.up(j1):0:8,system.tab,FVA.ModelStat,system.tab,FVA.Solvestat;

	PUTCLOSE;

);
