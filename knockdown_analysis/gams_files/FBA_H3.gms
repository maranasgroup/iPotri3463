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
$include "rxns_H3.txt"

        i   Metabolites of k-ath model
$include "mets_H3.txt"

;

*add an alias needed for FVA
alias (j1, j);

PARAMETERS

  S(i,j) 		Stoichiometric matrix for ED
$include "Sij_H3.txt"

	LB(j)			lower bounds of reactions
$include "LB_H3.txt"

	UB(j)			lower bounds of reactions
$include "UB_H3.txt"

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
	ngam1l			constrains NGAM to be equal across tissues in the light phase
	ngam2l			constrains NGAM to be equal across tissues in the light phase
	ngam1d			constrains NGAM to be equal across tissues in the dark phase
	ngam2d			constrains NGAM to be equal across tissues in the dark phase
	ngam3				constrains NGAM to be equal across light and dark phases
;

****************************** Equations*******************************
obj..		  					z =e= v('RXN_ATPM_c__S___DARK');
mb(i)..     				sum(j,S(i,j)*v(j)) =e= 0;
eqn14..							v('RXN_L2CP2_CO2_LIGHT') + 0.027 * v('RXN_E2L_CO2_LIGHT') =g= 0;
eqn15..							v('RXN_E2R_CO2_LIGHT') - v('RXN_R2CP1_CO2_LIGHT') =e= 0;
eqn15_1..						v('RXN_E2R_CO2_DARK') - v('RXN_R2CP1_CO2_DARK') =e= 0;
eqn17..							0.1 * v('RXN_E2S_CO2_DARK') - v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn17_1..						v('SK_MET_carbon_dioxide_c__S___LIGHT') + v('DM_MET_carbon_dioxide_c__S___DARK') =e= 0;
eqn18..							v('DM_MET_starch_p__L___LIGHT') + v('SK_MET_starch_p__L___DARK') =e= 0;
eqn19..							v('DM_MET_starch_p__S___LIGHT') + v('SK_MET_starch_p__S___DARK') =e= 0;
eqn21..							v('RXN_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p__L___LIGHT') + 2 * v('RXN_RXN_961_p__L___LIGHT') =e= 0;
ngam1l..						v('RXN_ATPM_c__S___LIGHT') - v('RXN_ATPM_c__L___LIGHT') =e= 0;
ngam2l..						v('RXN_ATPM_c__L___LIGHT') - v('RXN_ATPM_c__R___LIGHT') =e= 0;
ngam1d..						v('RXN_ATPM_c__S___DARK') - v('RXN_ATPM_c__L___DARK') =e= 0;
ngam2d..						v('RXN_ATPM_c__L___DARK') - v('RXN_ATPM_c__R___DARK') =e= 0;
ngam3..							v('RXN_ATPM_c__L___DARK') - v('RXN_ATPM_c__L___LIGHT') =e= 0;
*****************************************************************************


model FBA
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
	ngam1l
	ngam2l
	ngam1d
	ngam2d
	ngam3
/
;
*note: all the puclose and append stuff allows for real-time writing
file RESULT /FBA_result_H3.txt/;
put RESULT;

FBA.optfile = 1;

SOLVE FBA USING LP MAXIMIZING z;

file outputfname /"results.txt"/;
put outputfname;
put "Rxn,Flux"/;
loop(j,
    put j.tl:0:300,",",v.l(j):0:8/;
)
putclose outputfname;
