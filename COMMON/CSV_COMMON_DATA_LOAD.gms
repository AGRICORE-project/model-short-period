*------------------------------------------------------------------------------------------------------------------------------
* DEFININITION OF GENERAL EXECUTION PARAMETERS
*------------------------------------------------------------------------------------------------------------------------------

* Disable listing of all simbols in output file
$OFFSYMLIST OFFSYMXREF

* Disable inclusion of model rows/cols in the output file
OPTION LIMROW  = 0
OPTION LIMCOL  = 0

* Set model limit solving time in seconds: 60*20=1200 seconds; 2 minutes
OPTION RESLIM  = 1200

* Set model limit iterations for solving
OPTION iterlim = 200

* Selecting solvers for lp and nlp
option lp       = cplex;
option nlp      = conopt4;
  
*------------------------------------------------------------------------------------------------------------------------------
* DECLARATION OF PARAMETER SETS
*------------------------------------------------------------------------------------------------------------------------------

* proc = included products in the data
SET whole_population_proc all the processes in the population
/
$include SET_whole_population_proc.txt
/;

* proc = included products in the data
SET proc(whole_population_proc) all the processes in the processed region
/
$include SET_proc.txt
/;

* jx = arable non-perennial crops and milk
SET jx arable non-perennial crops and milk
/
$include SET_Arable_procs.txt
/;

* jrr = crops used as fodder for livestock
$Ifi exist SET_LivestockFood_procs.txt $goto LivestockFood_present
set jrr reused fodder crops / 0 /;
$goto jrrcomplete
$label LivestockFood_present
set jrr reused fodder crops 
/
$include SET_LivestockFood_procs.txt
/;
$label jrrcomplete


* jll = set for only milk as a subset of proc

SET jll(proc) milk
/milk/;

* jvv = arable non-perennial crops (no milk)
set jvv arable non-perennial crops (no milk)
/
$include SET_vegetal_procs.txt
/;

* water = TBD - fixed
SET water water footprint /green
blue
grey/;

* emi = proc (products) for which emissions are available / calculated
SET emi  CO2 emissions
/
$include SET_emi.txt
/;

* whole_population_farm = ids of each of the whole population
SET whole_population_farm number of farms
/
$include SET_whole_population_farm.txt
/;

* farm = ids of each of the included farms
SET farm(whole_population_farm) number of farms
/
$include SET_farm.txt
/;

* jnff = crops that fix nitrogen
$Ifi exist SET_FixingNitrogen_procs.txt $goto FixingNitrogen_present
set jnnf  N-fixing crops / 0 /;
$goto jnnfcomplete
$label FixingNitrogen_present
set jnnf  N-fixing crops 
/
$include SET_FixingNitrogen_procs.txt
/;
$label jnnfcomplete

* jggr = crops that are meadows
$Ifi exist SET_MeadowsAndPastures_procs.txt $goto MeadowsAndPastures_present
set jggr  Meadows / 0 /;
$goto jggrcomplete
$label MeadowsAndPastures_present
set jggr  Meadows 
/
$include SET_MeadowsAndPastures_procs.txt
/;
$label jggrcomplete

* jr1 = rotation crops 1
$Ifi exist SET_Semipermanent_procs.txt $goto Semipermanent_present
set jrr1  rotation crops 1 / 0 /;
$goto jrr1complete
$label Semipermanent_present
set jrr1  rotation crops 1 
/
$include SET_Semipermanent_procs.txt
/;
$label jrr1complete

* jr2 = rotation crops 2
$Ifi exist SET_COPOrIndustrial_procs.txt $goto COPOrIndustrial_present
set jrr2  rotation crops 2 / 0 /;
$goto jrr2complete
$label COPOrIndustrial_present
set jrr2  rotation crops 2 
/
$include SET_COPOrIndustrial_procs.txt
/;
$label jrr2complete


* jr3 = rotation crops 3
$Ifi exist SET_Cereal_procs.txt $goto Cereal_present
set jrr3  rotation crops 3 / 0 /;
$goto jrr3complete
$label Cereal_present
set jrr3  rotation crops 3 
/
$include SET_Cereal_procs.txt
/;
$label jrr3complete

* STABLE = parameters related to livestock (just names)
set STABLE /
rebreeding_cows number of rebreeding cows (lsu) per farm
dairy_cows number of dairy cows per farm
tot_lsu capi da milk + rebreeding_cows
pr_milk milk price per ton
prod_milk production in tons
c_milk milk cost per ton
/;

* su = decoupled subsidies
SET su subsidies
/
$include SET_subsidies.txt
/;

* coupled policies
set coupled_policy/
$include coupled_subsidies_v2_names.txt
/;

* ab = agent based parameters (just names)
SET ab  agent based
/age
successors
successor_age
COD_RAGR
OTE
altitude
family_members number of family_members members
CURRENT_ASSET immediate liquidity
/;

* ra = available agrarian regions in the data
SET ra agrarian region
/
$include SET_agrarian_region.txt
/;

* tax_level = level of taxation on CO2 emission
SET tax_level level of taxation on CO2 emission
/
$include SET_tax_level.txt
/;

*------------------------------------------------------------------------------------------------------------------------------------------
* INPORT DATA (PARAMETERS)
*----------------------------------------------------------------------------------------------------------------------------------------

* HA1 = utilized agricultural area ha per each farm and product
TABLE HA1(farm,proc) utilized agricultural area ha
$onDelim
$include UAA.csv
$offDelim;
DISPLAY$(showdisplay > 0) HA1;

* XBAR1 = produced quantities tons per each farm and product
TABLE XBAR1(farm,proc) produced quantities tons
$onDelim
$include PROD.csv
$offDelim;
DISPLAY$(showdisplay > 0) XBAR1;

* PR1 = selling price per tons per each farm and product
TABLE PR1(farm,proc) selling price per tons
$onDelim
$include PR.csv
$offDelim;
DISPLAY$(showdisplay > 0) PR1;

* C1 = specific variable cost per tons per each farm and product
TABLE C1(farm,proc) specific variable cost per tons
$onDelim
$include C_VAR.csv
$offDelim;
DISPLAY$(showdisplay > 0) C1;

* LIVESTOCK = livestock information
TABLE livestock(farm,stable) livestock information
$onDelim
$include stable.csv
$offDelim;
DISPLAY$(showdisplay > 0) livestock;

* Rewrite data in Ha1, XBAR1, PR1 and C1 for milk product extracting the information from the livestock table
* As milk is hardcoded, this should be controlled for populations without milk
HA1(farm,"milk")$(proc("milk"))=livestock(farm,"dairy_cows");
DISPLAY$(showdisplay > 0) HA1;

XBAR1(farm,"milk")$(proc("milk"))=livestock(farm,"prod_milk");
DISPLAY$(showdisplay > 0) XBAR1;

PR1(farm,"milk")$(proc("milk"))=livestock(farm,"pr_milk");
DISPLAY$(showdisplay > 0) pr1;

C1(farm,"milk")$(proc("milk"))=livestock(farm,"c_milk");
DISPLAY$(showdisplay > 0) C1;

Parameter rebreeding_ratio(farm) ratio between the numer fo rebreeding cows on the dairy cows; 
rebreeding_ratio(farm)$ livestock(farm,"dairy_cows") = livestock(farm,"rebreeding_cows")/livestock(farm,"dairy_cows");
DISPLAY$(showdisplay > 0) rebreeding_ratio;


set rent_information /
area area in ha
rent_total_price rent price in €
/;

TABLE Land_rents(farm,farm,rent_information) rents information
$onDelim
$include rents.csv
$offDelim;
DISPLAY$(showdisplay > 0) Land_rents;

TABLE ha_land_rented(farm,farm) land rented in ha;
ha_land_rented(farm,farm)$(Land_rents(farm,farm,"area") gt 0)=Land_rents(farm,farm,"area");

TABLE COUPLED_SUBSIDIES_RATES(coupled_policy,proc) coupled subsidies per each policy for each product
$onDelim
$include COUPLED_SUBS_v2_rates.csv
$offDelim;
DISPLAY$(showdisplay > 0) COUPLED_SUBSIDIES_RATES;

* TABLE COUPLED_SUBSIDIES_PAYMENTS(farm,coupled_policy) coupled subsidies payments per each policy for each farm
* $onDelim
* $include COUPLED_SUBS_v2.csv
* $offDelim;
* DISPLAY$(showdisplay > 0) COUPLED_SUBSIDIES_PAYMENTS;



* (!) maybe eggs

* DATABASE ERROR CORRECTION - Don't uses, removes relevant data
* HA1(farm,proc)$((XBAR1(farm,proc) EQ 0) OR (C1(farm,proc) eq 0) OR (PR1(farm,proc) eq 0))=0;
* XBAR1(farm,proc)$((HA1(farm,proc) EQ 0) OR (C1(farm,proc) eq 0) OR (PR1(farm,proc) eq 0))=0;
* C1(farm,proc)$((XBAR1(farm,proc) EQ 0)  OR (HA1(farm,proc) eq 0) OR (PR1(farm,proc) eq 0))=0;
* PR1(farm,proc)$((XBAR1(farm,proc) EQ 0) OR (HA1(farm,proc) eq 0) OR (C1(farm,proc) eq 0))=0;
* HA1(farm,"milk") $ (LIVESTOCK(farm,"rebreeding_cows") eq 0)=0;
LIVESTOCK(farm,"rebreeding_cows")$(HA1(farm,"milk") eq 0)=0;
DISPLAY$(showdisplay > 0) c1, pr1,ha1,xbar1;

* ABM = agent based parameters
TABLE ABM(farm,ra,ab)
$onDelim
$include AGENT_BASED.csv
$offDelim;
DISPLAY$(showdisplay > 0) ABM;

* sus = subsidies received by each farm for each uncoupled subsidy
TABLE sus(farm,su)
$onDelim
$include SUBSIDIES.csv
$offDelim;
DISPLAY$(showdisplay > 0) sus;

* COD_RAGR = agrarian region code for each farm 
parameter COD_RAGR(farm,ra) agrarian region code/
$onDelim
$include COD_RAGR.csv
$offDelim
/;
DISPLAY$(showdisplay > 0) COD_RAGR;

* H20 = water footprint for each product for each type of water
TABLE H20(proc,water) crops water footprint
$onDelim
$include waterfootprint.csv
$offDelim;
DISPLAY$(showdisplay > 0) h20;

* EM = CO2 emissions for each product
PARAMETER EM(proc) co2 emission/
$onDelim
$include emission.csv
$offDelim
/;
DISPLAY$(showdisplay > 0) EM;

* nitro_manure = nitrogen emission from livestock manure for each type of cow (rebreeding / dairy)
PARAMETER nitro_manure(stable) nitrogen emission from livestock manure/
$onDelim
$include nitrogen.csv
$offDelim
/;

* rental_price = rental price of the land for each farm (price per hectare)
PARAMETER rental_price(farm) rental price of the land/
$onDelim
$include rental_price.csv
$offDelim
/;

*
PARAMETER tax_emission(tax_level) tax on emission/
$onDelim
$include tax_emission.csv
$offDelim
/;

PARAMETER greening_surface(farm) greening surface in ha
/
$onDelim
$include greening_surface.csv
$offDelim
/;


*-----------------------------------------------------------------------------------------------------------------------------
* SELECTING DATA OF THE CROPS CULTIVATED WITHIN THE REGION
* CALCULATING AVERAGE PARAMETERS AND OTHER PARAMETERS                     
*-------------------------------------------------------------------------------------------------------------------------------

* LivestockX = TBD (this is exactly as livestock, no changes, why?)
PARAMETER LIVESTOCKX(farm,stable);
LIVESTOCKX(farm,stable)=LIVESTOCK(farm,stable)
DISPLAY$(showdisplay > 0) LIVESTOCKX;

* Hax = UAA per farm per proc only for arable non-perennial crops and milk
PARAMETER HAx(farm,proc) regional farm uaa;
HAx(farm,proc)$(jx(proc)) = HA1(farm,proc);
DISPLAY$(showdisplay > 0) HAx;

* XBARx = production per farm per proc only for arable non-perennial crops and milk 
PARAMETER XBARx(farm,proc) farm regional production;
XBARx(farm,proc)$(jx(proc)) = XBAR1(farm,proc);
DISPLAY$(showdisplay > 0) XBARx;

* YIELD = yield per farm per proc only for arable non-perennial crops and milk
PARAMETER YIELD(farm,proc) farm regional yield ;
YIELD(farm,proc)$HAx(farm,proc)=XBARx(farm,proc)/HAx(farm,proc);
DISPLAY$(showdisplay > 0) YIELD;

* AVEYIELD = average yield per proc only for arable non-perennial crops and milk (excluding from the average the crops with yield = 0)
PARAMETER AVEYIELD(proc) average yield;
AVEYIELD(proc)$(SUM(farm,YIELD(farm,proc) gt 0)) =SUM(farm,YIELD(farm,proc))/SUM(farm$YIELD(farm,proc), 1);
DISPLAY$(showdisplay > 0) AVEYIELD;

* Replacing the yields = 0 with the average yield
YIELD(farm,proc)$((HAx(farm,proc) gt 0) AND (YIELD(farm,proc) EQ 0)) = AVEYIELD(proc);
DISPLAY$(showdisplay > 0) YIELD;

* PRx = prices per farm per proc only for arable non-perennial crops and milk
PARAMETER PRx(farm,proc) farm regional prices;
PRx(farm,proc)$(jx(proc)) = PR1(farm,proc);
DISPLAY$(showdisplay > 0) PRx;

* AVEPRx = average prices per proc only for arable non-perennial crops and milk (excluding from the average the crops with price = 0)
PARAMETER AVEPRx(proc) average prices;
AVEPRx(proc)$(SUM(farm,PRx(farm,proc) gt 0))=SUM(farm,PRx(farm,proc))/SUM(farm$PRx(farm,proc), 1);
DISPLAY$(showdisplay > 0) AVEPRx;

* Replacing the prices = 0 with the average price
PRx(farm,proc)$((HAx(farm,proc) gt 0) AND (PRx(farm,proc) EQ 0)) = AVEPRx(proc);
DISPLAY$(showdisplay > 0) PRx;

* Cx = costs per farm per proc only for arable non-perennial crops and milk
PARAMETER Cx(farm,proc) farm regional production;
Cx(farm,proc)$(jx(proc)) = C1(farm,proc);
DISPLAY$(showdisplay > 0) Cx;

* AVECx = average costs per proc only for arable non-perennial crops and milk (excluding from the average the crops with cost = 0)
PARAMETER AVECx(proc) average cost;
AVECx(proc)$(SUM(farm,Cx(farm,proc) gt 0))=SUM(farm,Cx(farm,proc))/SUM(farm$Cx(farm,proc), 1);
DISPLAY$(showdisplay > 0) AVECx;

* Replacing the variable costs = 0 with the average ones
Cx(farm,proc)$((HAx(farm,proc) gt 0) AND (Cx(farm,proc) EQ 0)) = AVECx(proc);
DISPLAY$(showdisplay > 0) Cx;

* SUSx = TBD (same as sus)
PARAMETER SUSx(farm,su) farm regional subsidies ;
SUSx(farm,su)=sus(farm,su);
DISPLAY$(showdisplay > 0) SUSx;

* uaa = total utilized area per farm; CUAA = classification of the farms based on the utilized area
Parameter uaa(farm) uaa, CUAA(farm) uaa class dimension;
uaa(farm)=sum(proc$(not sameas (proc,"milk")),ha1(farm,proc)); 
CUAA(farm)$ (uaa(farm) le 10)                         = 1;
CUAA(farm)$((uaa(farm) gt 10)  and (uaa(farm) le 20))  = 2;
CUAA(farm)$((uaa(farm) gt 20)  and (uaa(farm) le 50))  = 3;
CUAA(farm)$((uaa(farm) gt 50)  and (uaa(farm) le 100)) = 4;
CUAA(farm)$((uaa(farm) gt 100) and (uaa(farm) le 300)) = 5;
CUAA(farm)$ (uaa(farm) gt 300)                        = 6;
DISPLAY$(showdisplay > 0) uaa;
DISPLAY$(showdisplay > 0) CUAA;

* farmj = subset of farms with  UAA >0 and UAA <1000 AND altitude > 0
SET farmj(farm);
farmj(farm)=no;
farmj(farm)$(SUM(proc,HAx(farm,proc) gt 0) AND (uaa(farm)le 1000)AND (sum(ra,ABM(farm,ra,"altitude")gt 0)))=yes;
DISPLAY$(showdisplay > 0) farmj;

* number = number of farms in the subset farmj
parameter number;
number= card(farmj);
DISPLAY$(showdisplay > 0) number;

* TBD
Parameter selage(farm),selsucc(farm),selalti(farm),selanfi(farm), selote(farm),
         initial_liq(farm);
selage(farmj) = sum(ra,ABM(farmj,ra,"age"));
selsucc(farmj) = sum(ra,ABM(farmj,ra,"successors"));
selalti(farmj) =sum(ra,ABM(farmj,ra,"altitude"));
selanfi(farmj)= sum(ra,ABM(farmj,ra,"successor_age"));
selote(farmj)= sum(ra,ABM(farmj,ra,"ote"));
initial_liq(farm) = sum(ra,ABM(farm,ra,"CURRENT_ASSET"));
DISPLAY$(showdisplay > 0) selage,selsucc,selalti,selanfi,selote,initial_liq;

* N_milk = nitrogen produced by dairy cows per farm
PARAMETER N_milk(farm) tonnellate NITROGEN a livello aziendale prodotto da vacche da milk;
N_milk(farmj) $ (hax(farmj,"milk")) = Nitro_manure("dairy_cows") * hax(farmj,"milk");

* N_rebreeding = nitrogen produced by rebreeding cows per farm
PARAMETER N_rebreeding(farm) tonnellate NITROGEN a livello aziendale prodotto da vacche da RIMONTA;
N_rebreeding(farmj) $ (LIVESTOCKX(farmj,"rebreeding_cows")) = Nitro_manure("rebreeding_cows") * LIVESTOCKX(farmj,"rebreeding_cows");

* NITROGEN_LSU = nitrogen produced by LSU per farm
PARAMETER NITROGEN_LSU (farm);
NITROGEN_LSU(farmj) $ (hax(farmj,"milk")) = (N_milk(farmj)+N_rebreeding(farmj))/(hax(farmj,"milk"));

* TOT_NITROGEN = total nitrogen per farm
PARAMETER TOT_NITROGEN(farm);
TOT_NITROGEN(farmj) $(hax(farmj,"milk")) = NITROGEN_LSU(farmj)*(hax(farmj,"milk"));
DISPLAY$(showdisplay > 0) N_milk, N_rebreeding, NITROGEN_LSU, TOT_NITROGEN;

* TBD: Optional code to analise distribution of variables in the population
$ontext
*ANALISI STATISTICAL ANALYSIS
*execute_unload "analisi.gdx" selage
*execute 'gdxxrw.exe analisi.gdx par=selage rng=selage!a1';
*execute_unload "analisi.gdx" selsucc
*execute 'gdxxrw.exe analisi.gdx par=selsucc rng=selsucc!a1';
*execute_unload "analisi.gdx" selote
*execute 'gdxxrw.exe analisi.gdx par=selote rng=selote!a1';
*execute_unload "analisi.gdx" selalti
*execute 'gdxxrw.exe analisi.gdx par=selalti rng=selalti!a1';
DISPLAY$(showdisplay > 0) c_intens;
$offtext

*-----------------------------------------------------------------------------------------------------------------
* CREATING A SUBSET OF ALL PRODUCTION PROCESSES ACTIVATED IN THE REGION (SUBj1)
*-----------------------------------------------------------------------------------------------------------------

* j0= all processes with uaa>0
SET j0(proc);
j0(proc)=no;
j0(proc)$(SUM(farmj,HAx(farmj,proc) gt 0))=yes;

* HAxs = Hax only for activated products in the region (jx)
Parameter HAxs(farm,proc) uaa regionale per processi seminativi (subset jx);
HAxs(farmj,proc)$jx(proc)=Hax(farmj,proc);

* Subset j1: all arable crops and milk with uaa>0 (jx)
SET j1(proc);
j1(proc)=no;
j1(proc)$(SUM(farmj,HAxs(farmj,proc) gt 0))=yes;

* HAxr = Hax only for reused forages (jrr)
Parameter HAxr(farm,proc) regional uaa for reused forages (jrr);
HAxr(farmj,proc)$jrr(proc)=Hax(farmj,proc);

* Subset j2: contains reused forages (jrr)
SET j2(proc);
j2(proc)=no;
j2(proc)$(SUM(farmj,HAxr(farmj,proc) gt 0))=yes;

* Haxl = Hax only for milk production (jll)
Parameter HAxl(farm,proc) regional lsu for milk production (jll);
HAxl(farmj,proc)$jll(proc)=Hax(farmj,proc);

* Subset j3: contains milk production (jll)
SET j3(proc);
j3(proc)=no;
j3(proc)$(SUM(farmj,HAxl(farmj,proc) gt 0))=yes;

* HAxv = Hax only for marketable arable crops (jvv)
Parameter HAxv(farm,proc) regional marketable arablebcrops  (jvv);
HAxv(farmj,proc)$jvv(proc)=Hax(farmj,proc);

* Subset j4: contains marketable arable bcrops (jvv)
SET j4(proc);
j4(proc)=no;
j4(proc)$(SUM(farmj,HAxv(farmj,proc) gt 0))=yes;

* HAnf = Hax only for regional N-fixing crops (jnf)
Parameter HAnf(farm,proc) regional N-fixing crops (jnf);
HAnf(farmj,proc)$jnnf(proc)=HAxv(farmj,proc);

*Subset j5: contains Nfixing bcrops (jnf)
SET j5(proc);
j5(proc)=no;
j5(proc)$(SUM(farmj,HAnf(farmj,proc) gt 0))=yes;

* HAg = Hax only for regional grazing crops (jgr)
Parameter HAgr(farm,proc) regional grazing crops (jgr);
HAgr(farmj,proc)$jggr(proc)=HAxv(farmj,proc);

*Subset j6: contains graz crops
SET j6(proc);
j6(proc)=no;
j6(proc)$(SUM(farmj,HAgr(farmj,proc) gt 0))=yes;

* HAr1 = Hax only for regional rotation crops 1 (jr1)
Parameter HAr1(farm,proc) regional rot crops 1 (jr1);
HAr1(farmj,proc)$jrr1(proc)=HAxv(farmj,proc);

*Subset j7: contains rotation crops 1
SET j7(proc);
j7(proc)=no;
j7(proc)$(SUM(farmj,HAr1(farmj,proc) gt 0))=yes;

* HAr2 = Hax only for regional rotation crops 2 (jr2)
Parameter HAr2(farm,proc) regional rot crops 2 (jr2);
HAr2(farmj,proc)$jrr2(proc)=HAxv(farmj,proc);

*Subset j8: contains rotation crops 2
SET j8(proc);
j8(proc)=no;
j8(proc)$(SUM(farmj,HAr2(farmj,proc) gt 0))=yes;

* HAr3 = Hax only for regional rotation crops 3 (jr3)
Parameter HAr3(farm,proc) regional rot crops 3 (jr3);
HAr3(farmj,proc)$jrr3(proc)=HAxv(farmj,proc);

*Subset j9: contains rotation crops 3
SET j9(proc);
j9(proc)=no;
j9(proc)$(SUM(farmj,HAr3(farmj,proc) gt 0))=yes;

* creating a file subj1.pnr with all these subsets. This makes this code require 2 executions, creating the subj1.prn file in the first iteration
* and then reading it in the second iteration
$Ifi exist subj1.prn $goto subj1exists

File printj /subj1.prn/;
printj.ap=0;
* This is the maximum value admitted by GAMS. This imposes a limit to the number of characters that are valid for the names of the groups (20)
printj.sw=20;
printj.lw=20;
put printj;

put "SET j(proc) /"; put /;
Loop(proc$j0(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET js(j) /"; put /;
Loop(proc$j1(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jr(j) /"; put /;
Loop(proc$j2(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jl(j) /"; put /;
Loop(proc$j3(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jv(j) /"; put /;
Loop(proc$j4(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jnf(j) /"; put /;
Loop(proc$j5(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jgr(j) /"; put /;
Loop(proc$j6(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jr1(j) /"; put /;
Loop(proc$j7(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jr2(j) /"; put /;
Loop(proc$j8(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;

put "SET jr3(j) /"; put /;
Loop(proc$j9(proc),
put proc.tl; put/;
);
put "/;";
put/;put/;
put/;put/;
$EXIT

$label subj1exists
$include subj1.prn
