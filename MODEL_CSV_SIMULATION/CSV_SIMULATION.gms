*------------------------------------------------------------------------------------------------------------------------------
* LOAD EXECUTION OPTIONS FROM FILE
*------------------------------------------------------------------------------------------------------------------------------

$include options.gms

*------------------------------------------------------------------------------------------------------------------------------
* LOAD DATA FROM EXTERNAL FILES WITH COMMON SCRIPT
*------------------------------------------------------------------------------------------------------------------------------    
$include CSV_COMMON_DATA_LOAD.gms
$Ifi not exist subj1.prn $EXIT
$include CSV_PMP_COMMON_PART.gms

*-------------------------------------------------------------------------------------------------------------------------
*                              IMPORTING CALIBRATION RESULTS 
*-------------------------------------------------------------------------------------------------------------------------
* Filtering the simulation only to the included farms but yet using the whole calibration input files
ALIAS (whole_population_proc, whole_population_ccc)

* TBD
parameter qmatcal(whole_population_proc,whole_population_ccc)/
$onDelim
$include qmat.csv
$offDelim
/;

* TBD
parameter fucal(whole_population_farm,whole_population_proc)/
$onDelim
$include fu.csv
$offDelim
/;

* importing nlsolshx (ha)
PARAMETER nlsolshxcal(whole_population_farm,whole_population_proc)/
$onDelim
$include nlsolshx.csv
$offDelim
/;

*importing nsolx (xbar)
PARAMETER nsolxcal(whole_population_farm,whole_population_proc)/
$onDelim
$include nsolx.csv
$offDelim
/;

* TBD 
PARAMETER ccostqcal(whole_population_farm,whole_population_proc)/
$onDelim
$include ccostq.csv
$offDelim
/;

* TBD
PARAMETER costqacal(whole_population_farm)/
$onDelim
$include costqa.csv
$offDelim
/;

* TBD
PARAMETER nlobjcal(whole_population_farm)/
$onDelim
$include nlobj.csv
$offDelim
/;

* TBD
PARAMETER ADCCcal(whole_population_farm)/
$onDelim
$include ADCC.csv
$offDelim
/;

PARAMETER ADCC(farm);
PARAMETER costqa(farm);
PARAMETER nlobj(farm);
PARAMETER ccostq(farm, proc);
PARAMETER nsolx(farm, proc);
PARAMETER nlsolshx(farm, proc);
PARAMETER fu(farm, proc);

ADCC(farm)=ADCCcal(farm);
costqa(farm)=costqacal(farm);
nlobj(farm)=nlobjcal(farm);
ccostq(farm, proc)=ccostqcal(farm, proc);
nsolx(farm, proc)=nsolxcal(farm, proc);
nlsolshx(farm, proc)=nlsolshxcal(farm, proc);
fu(farm, proc)=fucal(farm, proc);

ALIAS (proc,ccc);
parameter qmat(proc,ccc);
qmat(proc,ccc)=qmatcal(proc,ccc);

DISPLAY$(showdisplay > 0) ccostq;
DISPLAY$(showdisplay > 0) costqa;
DISPLAY$(showdisplay > 0) nlobj;
DISPLAY$(showdisplay > 0) ADCC;
DISPLAY$(showdisplay > 0) qmat;
DISPLAY$(showdisplay > 0) fu;
DISPLAY$(showdisplay > 0) nlsolshx;
DISPLAY$(showdisplay > 0) nsolx;

*-------------------------------------------------------------------------------------------------------------------------
*                              FINANCIAL PARAMETERS 
*-------------------------------------------------------------------------------------------------------------------------

PARAMETER var_liq (farm) variation (inizial liquidity - total costs);
  var_liq(farmj)= initial_liq(farmj)- ADCC(farmj);
DISPLAY$(showdisplay > 0) var_liq;

PARAMETER LOAN (FARM) loan requested if var_liq is negative;
  LOAN(farmj) $(var_liq(farmj) LT 0) = abs(var_liq(farmj));
DISPLAY$(showdisplay > 0) loan;

*assuming interest rate of 5% and 6 months of anticipation time
PARAMETER INTEREST(farm) bank interest rate on the loan;
  INTEREST(farmj)$(var_liq(farmj) LT 0)= loan(farmj) * 0.05* 6;
DISPLAY$(showdisplay > 0) INTEREST;

parameter tot_rev(farm,proc) total revenue in one year of activity (related to proc);
  tot_rev(farmj,j)$xbar(farmj,j)= Pr(farmj,j)*xbar(farmj,j) + sum(coupled_policy, (sh(farmj, coupled_policy, j)*ha(farmj,j)));
DISPLAY$(showdisplay > 0) tot_rev;

parameter final_liq(farm) final liquidity at the end of the year;
*initial liquidity - total variable proc costs + loan - bank loan + year revenue + decoupled payment
  final_liq(farmj)= var_liq(farmj) +loan(farmj) - INTEREST(farmj) + sum(j,tot_rev(farmj,j)) + spp(farmj);
DISPLAY$(showdisplay > 0) final_liq;

*-------------------------------------------------------------------------------------------------------------------------
*                              SIMULATION PHASE 
*-------------------------------------------------------------------------------------------------------------------------

SCALAR
 BETA aid curtailment rate by modulation / 0 / ;

Parameter CO2(proc);
  CO2(j)=EM(j);

*-------------------------------------------------------------------------------------------------------------------------
* Definition of output tables
*-------------------------------------------------------------------------------------------------------------------------

* Simulation scenarios:
* s_base: base scenario
* s _cal: calibration scenario ( must be equal to s_base)
* s_land: land exchange constraints
* s_em co2 taxation
* s_csv
SET
t  simulazioni /s_base,s_cal,s_land,s_em, s_csv/;

* Tables per area
parameter gre_prezv(ra,farm,proc,t)       crop prices;
parameter gre_prodv(ra,farm,proc,t)       produced quantities ;
parameter gre_prodr(ra,farm,proc,t)       re-used quantities ;
parameter gre_prodt(ra,farm,proc,t)       total quantities ;
parameter gre_supv(ra,farm,proc,t)        activated area (ha) per sold crops;
parameter gre_supTOT (ra,farm,t)          structural change (land rented in or out);
parameter gre_h2o(ra,farm,proc,t)         water consumption ;
parameter gre_co2(ra,farm,proc,t)         co2 emission ;
parameter gre_nitrogen (ra,farm,t)        nitrogen emission from livestock manure;
parameter gre_subpa(ra,farm,t)            greening surface;
parameter gre_pombt(ra,farm,t)            shadow price of land;
parameter gre_plvv(ra,farm,t )            salable products ;
parameter gre_plra(ra,farm,t )            re-used products ;
parameter gre_plvt (ra,farm,t);
parameter gre_compv(ra,farm,t)            compensations for sold processes ;
parameter gre_comptl(ra,farm,t)           total gross compensations;
parameter gre_aiutidis_GREEN(ra,farm,t)   decoupled subs P_GREEN ;
parameter gre_aiutidis_BASIC(ra,farm,t)   decoupled subs P_BASIC ;
parameter gre_aiutidis_PENSION(ra,farm,t) decoupled subs P_PENSION;
parameter gre_aiutiacc(ra,farm,t)         coupled subsidies ;
parameter gre_comphv(ra,farm, proc,t)     subsidies per ha ;
parameter gre_coupled_subsidy_per_policy_and_proc(ra,farm,coupled_policy, proc,t) coupled subsidies per policy and proc;
parameter gre_coupled_subsidy_per_policy(ra,farm,coupled_policy, t) coupled subsidies per policy ;
parameter gre_comptn(ra,farm,t);
parameter gre_rdm(ra,farm,t) ;
parameter gre_costv(ra,farm,t)            costs of sold crops ;
parameter gre_costt(ra,farm,t)            total costs ;
parameter gre_obj(ra,farm,t)              gross margin ;
parameter gre_supt(ra,farm)               total arable area + greening;
parameter gre_uaa(ra,farm,t)              uaa total;
parameter gre_cuaa(ra,farm,t)             Classi per uaa dimension ;
parameter gre_intens(ra,farm,t)           intensità produzione milk;
parameter gre_ote (ra,farm,t)             dairy or arable type of farming;
parameter gre_age (ra,farm,t)             holders age;
parameter gre_var_liq (ra,farm,t)         initial liquidity -total costs;
parameter GRE_LOAN(RA,FARM,T);
parameter gre_fin_liq (ra,farm,t)         farm liquidity at the end of the year;
parameter gre_rented_land (ra,farm,t)     land rented in by the farm;
parameter gre_coupledsubs(ra,farm,proc,t) coupled subsidies (form pillar 1);
parameter gre_rebreeding_cows(ra,farm,t);
parameter gre_dairy_cows(ra,farm,t);
parameter gre_land_exchange_V(ra,farm,farm,t);
parameter gre_land_exchange_Z(ra,farm,farm,t);

* regional tables
parameter re_prezv(proc,t)                crop prices; ;
parameter re_prodv(proc,t)                produced quantities ;
parameter re_prodr(proc,t)                re-used quantities ;
parameter re_prodt(proc,t)                total quantities ;
parameter re_supv(proc,t)                 activated area (ha) per sold crops; 
parameter re_h2o(proc,t)                  recupero consumo acqua ;
parameter re_co2(proc,t)                  recupero emissioni co2 ;
parameter re_nitrogen (t)                 nitrogen emission;
parameter re_subpa(t)                     bpa surface;
parameter re_pombt(t)                     shadow price of land ;
parameter re_plvv(t )                     salable productions ;
parameter re_plra(t )                     salable productions ;
parameter re_plvt (t);
parameter re_compv(t)                     compensations for sold processes ;
parameter re_comptl(t)                    total gross compensatione ;
parameter re_aiutidis_GREEN(t)            decoupled subs P_GREEN;
parameter re_aiutidis_BASIC(t)            decoupled subs P_BASIC;
parameter re_aiutidis_PENSION(t)          decoupled subs PENSION;
parameter re_aiutiacc(t)                  coupled subs;
parameter re_comphv(proc,t)               subs per ha ;
parameter re_comptn(t);
parameter re_rdm(t)                       aiuti-ha proc. venduti ;
parameter re_costv(t)                     sold crops costs ;
parameter re_costt(t)                     total costs ;
parameter re_obj(t)                       gross margin ;
parameter re_rebreeding_cows(t);
parameter re_dairy_cows(t);
parameter re_land_exchange_V(farm,farm,t);
parameter re_land_exchange_Z(farm,farm,t);

*-------------------------------------------------------------------------------------------------------------------------
* TBD
*-------------------------------------------------------------------------------------------------------------------------
gre_prezv(ra,farmj,j,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0)) = Pr(farmj,j) ;
gre_prezv(ra,farmj,proc,"s_base")$(not j(proc)) = Pr(farmj,proc) ;
gre_prodv(ra,farmj,j,"s_base")$(((sum(ab, abm(farmj,ra,ab)) gt 0)) AND (((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)))= nsolx(farmj,j) ;
gre_prodr(ra,farmj,j,"s_base")$(((sum(ab, abm(farmj,ra,ab)) gt 0)) AND (     jr(j)) and ha(farmj,"milk"))                             = nsolx(farmj,j) ;
gre_prodv(ra,farmj,proc,"s_base")$(((sum(ab, abm(farmj,ra,ab)) gt 0)) AND not j(proc)) = xbar(farmj,proc) ;
gre_prodt(ra,farmj,proc,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0)) = gre_prodv(ra,farmj,proc,"s_base") +  gre_prodr(ra,farmj,proc,"s_base");
gre_supv(ra,farmj,j,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0)) = nlsolshx(farmj,j) ;                            
gre_supv(ra,farmj,proc,"s_base")$(not j(proc))  = ha(farmj,proc) ;
gre_obj(ra,farmj,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0))      = nlobj(farmj);
gre_rebreeding_cows(ra,farmj,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0)) = HA(farmj,"milk") * rebreeding_ratio(farmj);
gre_dairy_cows(ra,farmj,"s_base")$((sum(ab, abm(farmj,ra,ab)) gt 0))  = HA(farmj,"milk") ;

DISPLAY$(showdisplay > 0) gre_supv, gre_prodv, gre_prodr;

*-------------------------------------------------------------------------------------------------------------------------
* TBD
*-------------------------------------------------------------------------------------------------------------------------
PARAMETER
  aiutitot(farm) farm total compensation
  aiutitot2(farm) farm total compensation
  compALFA(farm) average comp per ha;

PARAMETER
  aidis(farm), harif(farm),bt,av(j),cbpa, prodtab(farm), av2(farm,j),al2(farm,j),bt2(farm),pv2(farm,j),shv2(farm,coupled_policy,j),
  sht1,sht2,sht3,SHvAC(farm,j),SHAC(farm,proc),Av2M(j),Pv2M(j),Mlatte,TPPRv(farm,j), avepr(j),
  ful(farm,j);

  Av2M(j)$ (sum(farmj, 1$a(farmj,j,"land"))gt 0) = sum(farmj, a(farmj,j,"land"))/sum(farmj, 1$a(farmj,j,"land"));
  ful(farmj,j)=fu(farmj,j);
  shvAC(farmj,j)=0;
  shAC(farmj,j)=0;
  aidis(farmj)  =0;
  harif(farmj)  =0;
  cbpa       =250/10000;
  Av2(farmj,j) = A(farmj,j,"land");
  Al2(farmj,jr)=al(farmj,jr);
  Bt2(farmj) = B(farmj,"land") ;
  TPPRv(farmj,j)= pr(farmj,j);
  Pv2(farmj,j)  = pr(farmj,j);
  avepr(j)$sum(farmj, 1$pr(farmj,j))= sum(farmj, pr(farmj,j))/sum(farmj, 1$pr(farmj,j));
  SHv2(farmj,coupled_policy, j) = SH(farmj, coupled_policy,j);
DISPLAY$(showdisplay > 0) bt2, pv2,av2,avepr;

PARAMETER INS(farm);
  INS(farmj)$bt2(farmj)= 30000/bt2(farmj);

SET hlegprocs hleg covered crops /"ALFA", "PROT", "ORG_ALFA", "ORG_PROT"/;

Parameter hleg(farm);
  hleg(farmj) = sum(proc$(hlegprocs(proc)), ha(farmj, proc));
*  hleg(farmj)= ha(farmj,"ALFA")+ha(farmj,"PROT")+ha(farmj,"ORG_ALFA")+ha(farmj,"ORG_PROT");
DISPLAY$(showdisplay > 0) hleg;

Parameter ras(farm,ra), rast(ra), rastt, rastn;
  ras(farmj,ra)$(sum(ab, abm(farmj,ra,ab)) gt 0)=1;
  ras(farmj,ra)$((sum(j, ha(farmj,j)) gt 0))=ras(farmj,ra);
  rast(ra)=sum(farmj,ras(farmj,ra));
  rastt=sum(ra, rast(ra));
  rastn=sum(ra,1$rast(ra));
DISPLAY$(showdisplay > 0) ras,rast,rastt,rastn;

parameter razj(ra,farm,proc);
  razj(ra,farmj,j)$ ((ras (farmj,ra) gt 0) AND (ha(farmj,j) gt 0)) = 1;

parameter raj(ra,proc);
  raj(ra,j) = sum(farmj,razj(ra,farmj,j));
DISPLAY$(showdisplay > 0) razj, raj;

*-------------------------------------------------------------------------------------------------------------------------
* TBD
*-------------------------------------------------------------------------------------------------------------------------
 ALIAS(farmj,n,m);

* inizializing permanent at 0
Parameter perm(farm);
  perm(farmj)=0;

parameter cj(j),cjs(j);

* setting variables to build scenario equations
VARIABLES   xv(farm,j)              produzioni vendute in unita
            xhv(farm,j)             superficie produzioni vendute
            xhvs(farm)
            xh2o(farm,j,water)
            xco2(farm,j)
            xnitrogen(farm)
            sgreen(farm)
            prof(farm)              profitti di area
            profa(farm)
            profb(farm)
            proftot
            Z(farm)
            V(farm)
            ZZ(farm,farm)           terra acquisita in affitto
            VV(farm,farm)           terra ceduta in affitto
;

POSITIVE VARIABLE xv,xhv,xhvs,
                  xh2o,xco2,xnitrogen,
                  SA(farm,j),SB(farm,j),
                  sgreen,
                  Z,V,ZZ,VV
                  ;

 EQUATIONS
           RESt(farm)             land constraint with set aside
           RESt1(farm)            land constraint with set aside
           REST2(farm)            non productive land constrint
           Reuse (farm,j)
           ACCOUNTv(farm,j)       activated land (product sold) constraint
           VSUP                   total surface constraint
           DIVERS0(farm)          greening
           DIVERS1(farm)
           DIVERS2(farm,j)
           DIVERS3(farm,j)
           VECO(farm)
           VROT(farm)
           VPERM1
           VPERM2
           VH2O(farm,j,water)     water consumption
           VCO2(farm,j)           co2 emission
           VNITROGEN(farm)        nitrogen emission
           PROFITTI0(farm)        profit function
           PROFITTI1(farm)       
           PROFITT
           PROFITTI2(farm)        profitti per area con v.locazione
           PROFITTI_EM(farm)      profitti con tassa sulle emissioni 
           LOC(farm)              land exchange constraints
           LOC1(farm)             
           LOCz(farm)
           LOCv(farm)
           LOCzv(farm,farm)
           LOCzv2(farm)
           LOCzv3(farm)
           LOCzv4
           ;

*only for arable crops
RESt(farmj)$Bt2(farmj)..    sum(jV$cj(jV), xhv(farmj,jV)) =E=  Bt2(farmj);
RESt1(farmj)$Bt2(farmj)..   sum(jV$cj(jV), xhv(farmj,jV)) + sgreen (farmj) =E=  Bt2(farmj);

*reused forages 
Reuse(farmj,jr)$xbar(farmj,"milk").. (al(farmj,jr)*xv(farmj,"milk"))-xv(farmj,jr) =L=0;

ACCOUNTv(farmj,j)$(Bt2(farmj) AND (ha(farmj,j) gt 0) AND cj(j))..
          Av2(farmj,j) * xv(farmj,j) - xhv(farmj,j) =E= 0 ;
          
VSUP.. sum((n,jV),xhv(n,jv)) + sum(n,sgreen(n)) =E= sum(n, Bt2(n));


* superficie totale = terra + terra affittata - terra ceduta
$offorder
 LOC(n)$Bt2(n)..  sum(jv$cj(jv), xhv(n,jv)) =E=  Bt2(n) + Z(n) - V(n);
 LOC1(n)$Bt2(n)..  sum(jv$cj(jv), xhv(n,jv)) + sgreen(n) =E=  Bt2(n)+ Z(n) - V(n);
 LOCz(n)$Bt2(n)..  Z(n) =E= sum(m,ZZ(n,m));
 LOCv(m)$Bt2(m)..  V(m) =E= sum(n,VV(n,m));

* The total land that a farm rents to another farm should be equal to the land that the second rents from the former
 LOCzv(n,m) $((Bt2(n) gt 0) AND (Bt2(m) gt 0)).. ZZ(n,m)- VV(n,m) =E= 0;

* A farm can either rent out or rent in but not both at the same time
 LOCzv2(n)$Bt2(n).. Z(n)* V(n) =E= 0;

* A land can only rents out up to 20% of its total land
 LOCzv3(n)$Bt2(n).. Z(n) =L= 0.2*bt2(n);

* Extra check for ensuring that at regional level the total rented in land is equal to the rented out
 LOCzv4..           sum(n$Bt2(n), Z(n)) - sum(n$Bt2(n), V(n)) =E= 0;

*greening
DIVERS0(farmj)$Bt2(farmj)..      xhvs(farmj) =E= sum(js$cjs(js), xhv(farmj,js));

DIVERS1(farmj)$Bt2(farmj)..   sum(js$cjs(js), xhv(farmj,js)) =G= sum(js$cjs(js), SA(farmj,js)+SB(farmj,js));

DIVERS2(farmj,js)$(Bt2(farmj) AND cjs(js))..               SA(farmj,js) =L= 0.75* xhvs(farmj);

DIVERS3(farmj,js)$(Bt2(farmj) AND cjs(js))..   SA(farmj,js) + SB(farmj,js) =L= 0.95* xhvs(farmj);

*crops rotation and permanent at agrarian region level
VPERM1 .. sum(farmj $Bt2(farmj), sum(j$jgr(j),xhv(farmj,j))) =G= sum(farmj $Bt2(farmj),perm(farmj)*(1-0.05));
VPERM2.. sum (farmj$Bt2(farmj),sum(j$jgr(j),xhv(farmj,j))) =L= sum(farmj $Bt2(farmj),perm(farmj));
VECO(farmj)$Bt2(farmj)..     sum(j$jnf(j),xhv(farmj,j))*0.7 + sgreen(farmj)
                     =G=  0.05 * [sum(jv$cj(jv), xhv(farmj,jv))-sum(j$jgr(j),xhv(farmj,j))];
VROT(farmj)$Bt2(farmj)..     sum(j$jr1(j),xhv(farmj,j))/5 + sum(j$jr2(j),xhv(farmj,j))
                     =E=
                             sum(j$jr3(j),xhv(farmj,j));

*environmental impacts
VH2O(farmj,j,water)$(Bt2(farmj) AND cj(j)).. xh2o(farmj,j,water) =E= xhv(farmj,j)*H20(j,water);
VCO2(farmj,j)$(Bt2(farmj) AND cj(j))..     xco2(farmj,j)     =E= xhv(farmj,j)*CO2(j);
VNITROGEN(farmj)$(Bt2(farmj) AND ha(farmj,"milk")).. xnitrogen(farmj) =E= NITROGEN_LSU(farmj)*xhv(farmj,"milk");

PROFITTI0(farmj)$Bt2(farmj)..
          sum(j$(cj(j)), xv(farmj,j)  * Pv2(farmj,j) )
        + sum(j$(cj(j)), xhv(farmj,j) * sum(coupled_policy, SHv2(farmj,coupled_policy, j)))
        + spp(farmj)
        - .5 * sum((j,k)$((cj(j)) AND (cj(k))), xv(farmj,j) * Qmat(j,k) * xv(farmj,k))
        - sum(j$(cj(j)), ful(farmj,j)*xv(farmj,j))
=E= prof(farmj) ;

PROFITTI1(farmj)$Bt2(farmj)..
          sum(j, xv(farmj,j)  * Pv2(farmj,j) )
        + sum(j, xhv(farmj,j) * sum(coupled_policy, SHv2(farmj,coupled_policy, j)))
        + spp(farmj)
        - sgreen(farmj)* cbpa
        - .5 * sum((j,k), xv(farmj,j) * Qmat(j,k) * xv(farmj,k))- sum(j, ful(farmj,j)*xv(farmj,j))
=E= prof(farmj) ;

PROFITTI2(n)$Bt2(n)..
          sum(j$(cj(j)), xv(n,j)  * Pv2(n,j) )
        + sum(j$(cj(j)), xhv(n,j) * sum(coupled_policy, SHv2(n,coupled_policy, j)))
        + V(n)*(rental_price(n)) - Z(n)*(rental_price(n))
        + spp(n)
        - sgreen(n)* cbpa
        - .5 * sum((j,k)$((cj(j)) AND (cj(k))), xv(n,j) * Qmat(j,k) * xv(n,k))
        - sum(j$(cj(j)), ful(n,j)*xv(n,j))
=E= prof(n) ;

PROFITTI_EM(n)$Bt2(n)..
          sum(j$(cj(j)), xv(n,j)  * Pv2(n,j) )
        + sum(j$(cj(j)), xhv(n,j) * sum(coupled_policy, SHv2(n,coupled_policy, j)))
        + V(n)*(rental_price(n)) - Z(n)*(rental_price(n))
        - sum(j$(cj(j)),xco2(n,j) * (tax_emission("tax1")/10000))
        + spp(n)
        - sgreen(n)* cbpa
        - .5 * sum((j,k)$((cj(j)) AND (cj(k))), xv(n,j) * Qmat(j,k) * xv(n,k))
        - sum(j$(cj(j)), ful(n,j)*xv(n,j))
=E= prof(n) ;

PROFitT..  sum((farmj)$bt2(farmj), prof(farmj)) =E= proftot ;
                       
*AGENT BASED RULE ON LAND EXCHANGE
Z.fx(n) $(Bt2(n) AND (selage(n) GT 65) AND (selsucc(n) EQ 0))=0;

*modello scenario s_cal
MODEL  MODQ0   /RESt, ACCOUNTv,
                vh2o,vco2,
                Reuse,
                PROFITTI0
                profitT/ ;

*modello scenario s_land
MODEL  MODQ1   /
*                RESt1,
                ACCOUNTv,
                vh2o,vco2
                VNITROGEN
                DIVERS0
                DIVERS1
                DIVERS2
                DIVERS3
                VPERM1
                VPERM2
                VROT
                VECO
                LOC1
                LOCz
                LOCv
                LOCzv
                LOCzv2
                LOCzv3
                LOCzv4
                VSUP
                Reuse
                PROFITTI2
                profitT/ ;

*modello scenario tassa emissioni
MODEL  MODQ6   /
*                RESt1,
                ACCOUNTv,
                vh2o,vco2,
                VNITROGEN,
                DIVERS0
                DIVERS1
                DIVERS2
                DIVERS3
                VPERM1
                VPERM2
                VROT
                VECO
                LOC1
                LOCz
                LOCv
                LOCzv
                LOCzv2
                LOCzv3
                LOCzv4
                VSUP
                Reuse
                PROFITTI_EM
                profitT/;

*-------------------------------------------------------------------------------------------------------------------------
* TBD
*-------------------------------------------------------------------------------------------------------------------------

Loop(ra $((rast(ra) ge 1)),
* fase 2
* Identificazione soluzione modello quadratico con set asid e implicito MODQ

    bt2(farmj)=0;
    AV2(farmj,j)=0;
    Al2(farmj,jr)=0;
    Pv2(farmj,j)=0;
    SHv2(farmj,coupled_policy, j)=0;
    spp(farmj)=0;
    FUL(farmj,j)=0;
    Bt2(farmj)   $(abm(farmj,ra,"age") gt 0) = B(farmj,"land");
    Av2(farmj,j) $(abm(farmj,ra,"age") gt 0) = A(farmj,j,"land");
    Al2(farmj,jr) $(abm(farmj,ra,"age") gt 0) = AL(farmj,jr);
    Pv2(farmj,j) $(abm(farmj,ra,"age") gt 0) = pr(farmj,j);
    SHv2(farmj,coupled_policy, j)$(abm(farmj,ra,"age") gt 0) = SH(farmj,coupled_policy, j);
    spp(farmj)   $(abm(farmj,ra,"age") gt 0) = (susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    ful(farmj,j) $(abm(farmj,ra,"age") gt 0) = fu(farmj,j);
    cj(j)=0;
    beta        =0;
    aidis(farmj)   =0;
    modq0.optfile =1;
    cj(j)=raj(ra,j);

* attivare solo i processi presenti nella agrarian region
    xhv.l(farmj,j)=0;
    xv.l(farmj,j)=0;
    xco2.l(farmj,j)=0;
    xh2o.l(farmj,j,water)=0;
    prof.l(farmj)=0;
    sgreen.l(farmj)=0;
    rest.m(farmj)=0;
    xv.fx(farmj,j)$(xbar(farmj,j) eq 0) =0 ;
    xhv.fx(farmj,j)$(ha(farmj,j) eq 0) =0 ;

    SOLVE MODQ0 using nlp maximizing proftot ;

* lettura soluzioni
    gre_prezv(ra,farmj,j,"S_cal") = Pv2(farmj,j) ;
    gre_prezv(ra,farmj,proc,"s_cal")$(not j(proc)) = Pr(farmj,proc) ;
    gre_prodv(ra,farmj,j,"s_cal")$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)) = xv.l(farmj,j) ;
    gre_prodr(ra,farmj,j,"s_cal")$((     jr(j)) and ha(farmj,"milk"))                              = xv.l(farmj,j) ;
    gre_prodv(ra,farmj,proc,"s_cal")$(not j(proc)) = xbar(farmj,proc) ;
    gre_prodt(ra,farmj,proc,"s_cal") = gre_prodv(ra,farmj,proc,"s_cal") +  gre_prodr(ra,farmj,proc,"s_cal");
    gre_supv(ra,farmj,j,"S_cal")  = xhv.l(farmj,j) ;
    gre_supv(ra,farmj,proc,"s_cal")$(not j(proc))  = ha(farmj,proc) ;
    gre_h2o(ra,farmj,j,"S_cal")   = sum(water, xh2o.l(farmj,j,water)) ;
    gre_h2o(ra,farmj,proc,"s_cal")$(not j(proc))  = sum(water, ha(farmj,proc)*H20(proc,water));  
    gre_co2(ra,farmj,j,"S_cal")   = xco2.l(farmj,j) ;
    gre_co2(ra,farmj,proc,"s_cal")$(not j(proc))  = ha(farmj,proc)*CO2(proc); ;
    gre_nitrogen(ra,farmj,"s_cal") = NITROGEN_LSU(farmj)*xhv.l(farmj,"milk");
    gre_subpa(ra,farmj,"s_cal") = 0;
* aggiungo soluzione per espHORTare cambiamento strutturale dell'azienda agricola ma senza vincoli di scambio terra 
    gre_supTOT(ra,farmj,"s_cal") = Bt2(farmj)+ sum(proc$(not j(proc)), ha(farmj,proc));
    gre_plvv(ra,farmj,"S_cal")     = sum(j$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)), xv.l(farmj,j) * Pr(farmj,j) ) + sum(proc$(not j(proc)), xbar(farmj,proc) * Pr(farmj,proc) );
    gre_plra(ra,farmj,"S_cal")     = sum(j$((     jr(j)) and ha(farmj,"milk"))                             , xv.l(farmj,j) * avepr(j) ); 
    gre_plvt (ra,farmj,"S_cal")    = gre_plvv(ra,farmj,"S_cal") +   gre_plra(ra,farmj,"S_cal");
    gre_compv(ra,farmj,"S_cal")    = sum((j,coupled_policy), SH(farmj,coupled_policy,j) *  xhv.l(farmj,j))  + sum((proc,coupled_policy)$(not j(proc)), SH(farmj, coupled_policy,proc) *  ha(farmj,proc))
                                   + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc));
    gre_comptl(ra,farmj,"S_cal")$ras(farmj,ra) = sum((j, coupled_policy), SH(farmj,coupled_policy, j) *  xhv.l(farmj,j))  + sum((proc,coupled_policy)$(not j(proc)), SH(farmj,coupled_policy, proc) *  ha(farmj,proc))
                                               + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc))
                                               + (susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    gre_aiutidis_GREEN(ra,farmj,"s_cal")$ras(farmj,ra)    =  (susx(farmj,"P_GREEN")/10000);
    gre_aiutidis_BASIC(ra,farmj,"s_cal")$ras(farmj,ra)    =  (susx(farmj,"P_BASIC"))/10000;
    gre_aiutidis_PENSION(ra,farmj,"s_cal")$ras(farmj,ra)    =  (pens(farmj))/10000;
    gre_coupledsubs(ra,farmj,proc,"s_cal")$(not j(proc))   =  ( sum(coupled_policy, SH(farmj,coupled_policy,proc)) *  ha(farmj,proc)
                                 + ShAC(farmj,proc) * xbar(farmj,proc))*10000;
    gre_comptn(ra,farmj,"S_cal")   = gre_comptl(ra,farmj,"S_cal") ;
    gre_rdm(ra,farmj,"s_cal") = 0;
    gre_comphv(ra,farmj,j,"S_cal")$ras(farmj,ra)    = sum(coupled_policy, SHv2(farmj, coupled_policy, j));
    gre_coupled_subsidy_per_policy_and_proc(ra,farmj,coupled_policy,j,"S_cal")$ras(farmj,ra) = SHv2(farmj, coupled_policy, j) * ha(farmj,j) ;
    gre_coupled_subsidy_per_policy(ra,farmj,coupled_policy,"S_cal")$ras(farmj,ra)    = sum(j, SHv2(farmj, coupled_policy, j) * ha(farmj,j)) ;
    gre_costv(ra,farmj,"S_cal") =    .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_costt(ra,farmj,"S_cal") =    .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_obj(ra,farmj,"S_cal")      = prof.l(farmj) ;
    gre_rebreeding_cows(ra,farmj,"s_cal")= HA(farmj,"milk") * rebreeding_ratio(farmj);
    gre_dairy_cows(ra,farmj,"s_cal") = HA(farmj,"milk") ;
    gre_ote(ra,farmj,"s_cal")= ABM(farmj,ra,"OTE");
    gre_age(ra,farmj,"s_cal")= ABM(farmj,ra,"AGE");
    gre_var_liq (ra,farmj,"s_cal")$ras(farmj,ra)= initial_liq(farmj)- gre_costv(ra,farmj,"S_cal");
    gre_loan (ra,farmj,"s_cal")$(gre_var_liq(ra,farmj,"s_cal") LT 0) = abs(gre_var_liq(ra,farmj,"s_cal"));
    gre_loan (ra,farmj,"s_cal")$(gre_var_liq(ra,farmj,"s_cal") GT 0) = 0;
    gre_fin_liq (ra,farmj,"s_cal")= gre_var_liq (ra,farmj,"s_cal")
                                  +  gre_loan (ra,farmj,"s_cal") - gre_loan (ra,farmj,"s_cal") * 0.04* 6
                                  + sum(j, xv.l(farmj,j) * Pr(farmj,j)) +  gre_comptl(ra,farmj,"S_cal");
    gre_uaa(ra,farmj,"s_cal")=sum(jV, gre_supv(ra,farmj,jV, "s_cal"));
    gre_cuaa(ra,farmj,"s_cal")$ (gre_uaa(ra,farmj,"s_cal") le 10)                                       = 1;
    gre_cuaa(ra,farmj,"s_cal")$((gre_uaa(ra,farmj,"s_cal") gt 10)  and (gre_uaa(ra,farmj,"s_cal") le 20))  = 2;
    gre_cuaa(ra,farmj,"s_cal")$((gre_uaa(ra,farmj,"s_cal") gt 20)  and (gre_uaa(ra,farmj,"s_cal") le 50))  = 3;
    gre_cuaa(ra,farmj,"s_cal")$((gre_uaa(ra,farmj,"s_cal") gt 50)  and (gre_uaa(ra,farmj,"s_cal") le 100)) = 4;
    gre_cuaa(ra,farmj,"s_cal")$((gre_uaa(ra,farmj,"s_cal") gt 100) and (gre_uaa(ra,farmj,"s_cal") le 300)) = 5;
    gre_cuaa(ra,farmj,"s_cal")$ (gre_uaa(ra,farmj,"s_cal") gt 300)                                      = 6;
    gre_cuaa(ra,farmj,"s_cal")$ (gre_uaa(ra,farmj,"s_cal") eq 0)                                        = 0;

* scaling (xbar*1000)
    gre_intens(ra,farmj,"s_cal")$(Bt2(farmj)AND gre_supTOT(ra,farmj,"s_cal")GT 0)=gre_prodv(ra,farmj,"milk", "s_cal")*1000/gre_supTOT(ra,farmj,"s_cal");
);
DISPLAY$(showdisplay > 0) gre_supv;


*-------------------------------------------------------------------------------------------------------------------------
* TBD
*-------------------------------------------------------------------------------------------------------------------------
*modq1 s_land

* fase 2
* RICELUZIONE MODELLI

Loop(ra $(rast(ra) ge 1),
    bt2(farmj)=0;
    AV2(farmj,j)=0;
    Al2(farmj,jr)=0;
    Pv2(farmj,j)=0;
    SHv2(farmj,coupled_policy,j)=0;
    spp(farmj)=0;
    FUL(farmj,j)=0;

    Bt2(farmj)   $(abm(farmj,ra,"age") gt 0) = B(farmj,"land");
    Av2(farmj,j) $(abm(farmj,ra,"age") gt 0) = A(farmj,j,"land");
    Al2(farmj,jr) $(abm(farmj,ra,"age") gt 0) = AL(farmj,jr);
    Pv2(farmj,j) $(abm(farmj,ra,"age") gt 0) = pr(farmj,j);
    SHv2(farmj,coupled_policy,j)$(abm(farmj,ra,"age") gt 0) = SH(farmj,coupled_policy,j);
    spp(farmj)   $(abm(farmj,ra,"age") gt 0) = (susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    ful(farmj,j) $(abm(farmj,ra,"age") gt 0) = fu(farmj,j);

    VV.l(n,m)$(abm(n,ra,"age") gt 0 and abm(m,ra,"age") gt 0 ) = ha_land_rented(n,m);
    ZZ.l(n,m)$(abm(n,ra,"age") gt 0 and abm(m,ra,"age") gt 0 ) = ha_land_rented(m,n);


    perm(farmj)$(abm(farmj,ra,"age") gt 0)=gre_supv(ra,farmj,"GRAZ","s_cal")+gre_supv(ra,farmj,"ORG_GRAZ","s_cal");

* controllo di j
    cj(j)=raj(ra,j);
    cjs(j)$js(j)=raj(ra,j);

*i processi non attivati nella agrarian region non vengono calcolati
*xv.fx(farmj,j) $  (raj(ra,j) eq 0) = 0;
*xhv.fx(farmj,j) $  (raj(ra,j) eq 0) = 0;

    SOLVE MODQ1 using nlp maximizing proftot ;

* lettura soluzioni
    gre_prezv(ra,farmj,j,"S_land") = Pv2(farmj,j) ;
    gre_prezv(ra,farmj,proc,"s_land")$(not j(proc)) = Pr(farmj,proc) ;
    gre_prodv(ra,farmj,j,"s_land")$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)) = xv.l(farmj,j) ;
    gre_prodr(ra,farmj,j,"s_land")$((     jr(j)) and ha(farmj,"milk"))                              = xv.l(farmj,j) ;
    gre_prodv(ra,farmj,proc,"s_land")$(not j(proc)) = xbar(farmj,proc) ;
    gre_prodt(ra,farmj,proc,"s_land") = gre_prodv(ra,farmj,proc,"s_land") +  gre_prodr(ra,farmj,proc,"s_land");
    gre_supv(ra,farmj,j,"S_land")  = xhv.l(farmj,j) ;
    gre_supv(ra,farmj,proc,"s_land")$(not j(proc))  = ha(farmj,proc) ;
    gre_land_exchange_V(ra,n,m,"s_land") = VV.l(n,m);
    gre_land_exchange_Z(ra,n,m,"s_land") = ZZ.l(n,m);
    gre_h2o(ra,farmj,j,"S_land")   = sum(water, xh2o.l(farmj,j,water)) ;
    gre_h2o(ra,farmj,proc,"s_land")$(not j(proc))  = sum(water, ha(farmj,proc)*H20(proc,water));  
    gre_co2(ra,farmj,j,"S_land")   = xco2.l(farmj,j) ;
    gre_co2(ra,farmj,proc,"s_land")$(not j(proc))  = ha(farmj,proc)*CO2(proc); ;
    gre_nitrogen(ra,farmj,"s_land") = NITROGEN_LSU(farmj)*xhv.l(farmj,"milk");
    gre_subpa(ra,farmj,"s_land") = 0;
*aggiungo soluzione per esportare cambiamento strutturale dell'azienda agricola in base all'affitto di terra (vloc) 
    gre_supTOT(ra,farmj,"s_land") = Bt2(farmj)+ Z.l(farmj) - V.l(farmj) + sum(proc$(not j(proc)), ha(farmj,proc));
    gre_rented_land (ra,farmj,"s_land")=+ Z.l(farmj) - V.l(farmj);
    gre_plvv(ra,farmj,"S_land")     = sum(j$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)), xv.l(farmj,j) * Pr(farmj,j) ) + sum(proc$(not j(proc)), xbar(farmj,proc) * Pr(farmj,proc) );
    gre_plra(ra,farmj,"S_land")     = sum(j$((     jr(j)) and ha(farmj,"milk"))                             , xv.l(farmj,j) * avepr(j) );
    gre_plvt (ra,farmj,"S_land")    = gre_plvv(ra,farmj,"S_land") + gre_plra(ra,farmj,"S_land");
    gre_compv(ra,farmj,"S_land")    = sum((j, coupled_policy), SH(farmj,coupled_policy, j) *  xhv.l(farmj,j))  + sum((proc, coupled_policy)$(not j(proc)), SH(farmj,coupled_policy, proc) *  ha(farmj,proc))
                                    + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc));
    gre_comptl(ra,farmj,"S_land")$ras(farmj,ra) = sum((j, coupled_policy), SH(farmj,coupled_policy,j) *  xhv.l(farmj,j))  + sum((proc, coupled_policy)$(not j(proc)), SH(farmj,coupled_policy,proc) *  ha(farmj,proc))
                                 + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc))
                                 +(susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    gre_aiutidis_GREEN(ra,farmj,"s_land")$ras(farmj,ra)    =  (susx(farmj,"P_GREEN")/10000);
    gre_aiutidis_BASIC(ra,farmj,"s_land")$ras(farmj,ra)    =  (susx(farmj,"P_BASIC"))/10000;
    gre_aiutidis_PENSION(ra,farmj,"s_land")$ras(farmj,ra)    =  (pens(farmj))/10000;
    gre_coupledsubs(ra,farmj,proc,"s_land")$(not j(proc))   =  ( sum(coupled_policy,SH(farmj,coupled_policy,proc)) *  ha(farmj,proc)
                                  + ShAC(farmj,proc) * xbar(farmj,proc))*10000;
    gre_comptn(ra,farmj,"s_land")   = gre_comptl(ra,farmj,"s_land") ;
    gre_rdm(ra,farmj,"s_land") = 0;
    gre_comphv(ra,farmj,j,"s_land") $ras(farmj,ra)   = sum(coupled_policy, SHv2(farmj,coupled_policy,j));
    gre_coupled_subsidy_per_policy_and_proc(ra,farmj,coupled_policy,j,"s_land")$ras(farmj,ra) = SHv2(farmj, coupled_policy, j) * ha(farmj,j) ;
    gre_coupled_subsidy_per_policy(ra,farmj,coupled_policy,"s_land")$ras(farmj,ra)    = sum(j, SHv2(farmj, coupled_policy, j) * ha(farmj,j)) ;

*variable cost per crop
    gre_costv(ra,farmj,"s_land") =   .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_costt(ra,farmj,"s_land") =    .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_var_liq (ra,farmj,"s_land")$ras(farmj,ra)= initial_liq(farmj)- gre_costv(ra,farmj,"s_land");
    gre_loan (ra,farmj,"s_land")$(gre_var_liq(ra,farmj,"s_land") LT 0) = abs(gre_var_liq(ra,farmj,"s_land"));
    gre_loan (ra,farmj,"s_land")$(gre_var_liq(ra,farmj,"s_land") GT 0) = 0;
    gre_fin_liq (ra,farmj,"s_land") = gre_var_liq (ra,farmj,"s_land")
                                    +  gre_loan (ra,farmj,"s_land") - gre_loan (ra,farmj,"s_land") * 0.04* 6
                                    + sum(j, xv.l(farmj,j) * Pr(farmj,j)) + gre_comptl(ra,farmj,"s_land")
                                    - Z.l(farmj)*(rental_price(farmj)) + V.l(farmj)*(rental_price(farmj));
    gre_obj(ra,farmj,"s_land")      = prof.l(farmj) ;

* TBD - Michele, this position seems strange    
    xhv.l(farmj,j)=0;
    xv.l(farmj,j)=0;
    xco2.l(farmj,j)=0;
    xnitrogen.l(farmj)=0;
    xh2o.l(farmj,j,water)=0;
    prof.l(farmj)=0;
    z.l(farmj)=0;
    v.l(farmj)=0;
    sgreen.l(farmj)=0;
  
    gre_rebreeding_cows(ra,farmj,"s_land")= HA(farmj,"milk") * rebreeding_ratio(farmj);
    gre_dairy_cows(ra,farmj,"s_land") = HA(farmj,"milk") ;
    gre_ote(ra,farmj,"s_land")= ABM(farmj,ra,"OTE");
    gre_age(ra,farmj,"s_land")= ABM(farmj,ra,"AGE");
    gre_uaa(ra,farmj,"s_land")=sum(jV, gre_supv(ra,farmj,jV, "s_land"));
    gre_cuaa(ra,farmj,"s_land")$ (gre_uaa(ra,farmj,"s_land") le 10)                                       = 1;
    gre_cuaa(ra,farmj,"s_land")$((gre_uaa(ra,farmj,"s_land") gt 10)  and (gre_uaa(ra,farmj,"s_land") le 20))  = 2;
    gre_cuaa(ra,farmj,"s_land")$((gre_uaa(ra,farmj,"s_land") gt 20)  and (gre_uaa(ra,farmj,"s_land") le 50))  = 3;
    gre_cuaa(ra,farmj,"s_land")$((gre_uaa(ra,farmj,"s_land") gt 50)  and (gre_uaa(ra,farmj,"s_land") le 100)) = 4;
    gre_cuaa(ra,farmj,"s_land")$((gre_uaa(ra,farmj,"s_land") gt 100) and (gre_uaa(ra,farmj,"s_land") le 300)) = 5;
    gre_cuaa(ra,farmj,"s_land")$ (gre_uaa(ra,farmj,"s_land") gt 300)                                      = 6;
    gre_cuaa(ra,farmj,"s_land")$ (gre_uaa(ra,farmj,"s_land") eq 0)                                        = 0;
    gre_intens(ra,farmj,"s_land")$(Bt2(farmj)AND gre_supTOT(ra,farmj,"s_land")GT 0)=gre_prodv(ra,farmj,"milk", "s_land")*1000/gre_supTOT(ra,farmj,"s_land");
);
DISPLAY$(showdisplay > 0) gre_supv,gre_obj;

*-------------------------------------------------------------------------------------------------------------------
*modq6 s_em
*xv.l(farmj,j)=nsolx(farmj,j);
* fase 2
*RISOLUZIONE MODELLI

Loop(ra $(rast(ra) ge 1),
    bt2(farmj)=0;
    AV2(farmj,j)=0;
    Al2(farmj,jr)=0;
    Pv2(farmj,j)=0;
    SHv2(farmj,coupled_policy,j)=0;
    spp(farmj)=0;
    FUL(farmj,j)=0;

    Bt2(farmj)   $(abm(farmj,ra,"age") gt 0) = B(farmj,"land");
    Av2(farmj,j) $(abm(farmj,ra,"age") gt 0) = A(farmj,j,"land");
    Al2(farmj,jr) $(abm(farmj,ra,"age") gt 0) = AL(farmj,jr);
    Pv2(farmj,j) $(abm(farmj,ra,"age") gt 0) = pr(farmj,j);
    SHv2(farmj,coupled_policy,j)$(abm(farmj,ra,"age") gt 0) = SH(farmj,coupled_policy,j);
    spp(farmj)   $(abm(farmj,ra,"age") gt 0) = (susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    ful(farmj,j) $(abm(farmj,ra,"age") gt 0) = fu(farmj,j);

    VV.l(n,m)$(abm(n,ra,"age") gt 0 and abm(m,ra,"age") gt 0 ) = ha_land_rented(n,m);
    ZZ.l(n,m)$(abm(n,ra,"age") gt 0 and abm(m,ra,"age") gt 0 ) = ha_land_rented(m,n);
    sgreen.l(n)$(abm(n,ra,"age") gt 0) = greening_surface(n);

    perm(farmj)$(abm(farmj,ra,"age") gt 0)=gre_supv(ra,farmj,"GRAZ","s_cal")+gre_supv(ra,farmj,"ORG_GRAZ","s_cal");

*controllo di j
    cj(j)=raj(ra,j);
    cjs(j)$js(j)=raj(ra,j);

*i processi non attivati nella agrarian region non vengono calcolati
*xv.fx(farmj,j) $  (raj(ra,j) eq 0) = 0;
*xhv.fx(farmj,j) $  (raj(ra,j) eq 0) = 0;
    SOLVE MODQ6 using nlp maximizing proftot ;

* lettura soluzioni
    gre_prezv(ra,farmj,j,"s_em") = Pv2(farmj,j) ;
    gre_prezv(ra,farmj,proc,"s_em")$(not j(proc)) = Pr(farmj,proc) ;
    gre_prodv(ra,farmj,j,"s_em")$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)) = xv.l(farmj,j) ;
    gre_prodr(ra,farmj,j,"s_em")$((     jr(j)) and ha(farmj,"milk"))                              = xv.l(farmj,j) ;
    gre_prodv(ra,farmj,proc,"s_em")$(not j(proc)) = xbar(farmj,proc) ;
    gre_prodt(ra,farmj,proc,"s_em") = gre_prodv(ra,farmj,proc,"s_em") +  gre_prodr(ra,farmj,proc,"s_em");
    gre_supv(ra,farmj,j,"s_em")  = xhv.l(farmj,j) ;
    gre_supv(ra,farmj,proc,"s_em")$(not j(proc))  = ha(farmj,proc) ;
    gre_land_exchange_V(ra,n,m,"s_em") = VV.l(n,m);
    gre_land_exchange_Z(ra,n,m,"s_em") = ZZ.l(n,m);

* aggiungo soluzione sull'nitrogen
    gre_h2o(ra,farmj,j,"S_em")   = sum(water, xh2o.l(farmj,j,water)) ;
    gre_h2o(ra,farmj,proc,"s_em")$(not j(proc))  = sum(water, ha(farmj,proc)*H20(proc,water));  
    gre_co2(ra,farmj,j,"S_em")   = xco2.l(farmj,j) ;
    gre_co2(ra,farmj,proc,"s_em")$(not j(proc))  = ha(farmj,proc)*CO2(proc); ;
    gre_nitrogen(ra,farmj,"s_em") = xnitrogen.l(farmj);
    gre_subpa(ra,farmj,"s_em") = sgreen.l(farmj);
  
* aggiungo soluzione per espHORTare cambiamento strutturale dell'azienda agricola in base all'affitto di terra (vloc) 
    gre_supTOT(ra,farmj,"s_em") = Bt2(farmj)+ Z.l(farmj) - V.l(farmj) + sum(proc$(not j(proc)), ha(farmj,proc));
    gre_rented_land (ra,farmj,"s_em")=+ Z.l(farmj) - V.l(farmj);

* gre_pombt (ra,farmj,"s_em")   = rest.m(farmj) ;
    gre_plvv(ra,farmj,"S_em")     = sum(j$(((not jr(j)) and ha(farmj,"milk")) OR  (ha(farmj,"milk") eq 0)), xv.l(farmj,j) * Pr(farmj,j) ) + sum(proc$(not j(proc)), xbar(farmj,proc) * Pr(farmj,proc) );
    gre_plra(ra,farmj,"S_em")     = sum(j$((     jr(j)) and ha(farmj,"milk"))                             , xv.l(farmj,j) * avepr(j)  );
    gre_plvt (ra,farmj,"S_em")    = gre_plvv(ra,farmj,"S_em") + gre_plra(ra,farmj,"S_em") ;
    gre_compv(ra,farmj,"S_em")    = sum((j,coupled_policy), SH(farmj,coupled_policy,j) *  xhv.l(farmj,j))  + sum((proc,coupled_policy)$(not j(proc)), SH(farmj, coupled_policy,proc) *  ha(farmj,proc))
                                 + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc));
    gre_comptl(ra,farmj,"S_em")$ras(farmj,ra) = sum((j, coupled_policy), SH(farmj,coupled_policy, j) *  xhv.l(farmj,j))  + sum((proc,coupled_policy)$(not j(proc)), SH(farmj,coupled_policy, proc) *  ha(farmj,proc))
                                 + sum(j, ShAC(farmj,j) *  xv.l(farmj,j))  + sum(proc$(not j(proc)), ShAC(farmj,proc) * xbar(farmj,proc))
                                 +(susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+pens(farmj))/10000;
    gre_aiutidis_GREEN(ra,farmj,"s_em")$ras(farmj,ra)    =  (susx(farmj,"P_GREEN")/10000);
    gre_aiutidis_BASIC(ra,farmj,"s_em")$ras(farmj,ra)    =  (susx(farmj,"P_BASIC"))/10000;
    gre_aiutidis_PENSION(ra,farmj,"s_em")$ras(farmj,ra)    =  (pens(farmj))/10000;
    gre_coupledsubs(ra,farmj,j,"s_em")   =  (sum(coupled_policy, SH(farmj,coupled_policy, j))   *  xhv.l(farmj,j))
                                  + (ShAC(farmj,j) *  xv.l(farmj,j));
    gre_coupledsubs(ra,farmj,proc,"s_land")$(not j(proc))   =  ( sum(coupled_policy, SH(farmj,coupled_policy, proc)) *  ha(farmj,proc)
                                  + ShAC(farmj,proc) * xbar(farmj,proc))*10000;
    gre_comptn(ra,farmj,"s_em")   = gre_comptl(ra,farmj,"s_em") ;
    gre_rdm(ra,farmj,"s_em") = 0;
    gre_comphv(ra,farmj,j,"s_em")  $ras(farmj,ra)  = sum(coupled_policy, SHv2(farmj,coupled_policy,j));
    gre_coupled_subsidy_per_policy_and_proc(ra,farmj,coupled_policy,j,"s_em")$ras(farmj,ra) = SHv2(farmj, coupled_policy, j) * ha(farmj,j) ;
    gre_coupled_subsidy_per_policy(ra,farmj,coupled_policy,"s_em")$ras(farmj,ra)    = sum(j, SHv2(farmj, coupled_policy, j) * ha(farmj,j)) ;
    gre_costv(ra,farmj,"s_em") =    .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_costt(ra,farmj,"s_em") =    .5 * sum((j,k), xv.l(farmj,j) * qmat (j,k) * xv.l(farmj,k))+ sum(j, fu(farmj,j) * xv.l(farmj,j));
    gre_var_liq (ra,farmj,"s_em")$ras(farmj,ra)= initial_liq(farmj)- gre_costv(ra,farmj,"s_em");
    gre_loan (ra,farmj,"s_em")$(gre_var_liq(ra,farmj,"s_em") LT 0) = abs(gre_var_liq(ra,farmj,"s_em"));
    gre_loan (ra,farmj,"s_em")$(gre_var_liq(ra,farmj,"s_em") GT 0) = 0;
    gre_fin_liq (ra,farmj,"s_em") = gre_var_liq (ra,farmj,"s_em")
                                  +  gre_loan (ra,farmj,"s_em") - gre_loan (ra,farmj,"s_em") * 0.04* 6
                                  + sum(j, xv.l(farmj,j) * Pr(farmj,j)) +  gre_comptl(ra,farmj,"s_em")
                                  - Z.l(farmj)*(rental_price(farmj)) + V.l(farmj)*(rental_price(farmj));
    gre_obj(ra,farmj,"s_em")      = prof.l(farmj) ;

    xhv.l(farmj,j)=0;
    xv.l(farmj,j)=0;
    xco2.l(farmj,j)=0;
    xh2o.l(farmj,j,water)=0;
    xnitrogen.l(farmj)=0;
    prof.l(farmj)=0;
    z.l(farmj)=0;
    v.l(farmj)=0;
    sgreen.l(farmj)=0;
 
    gre_rebreeding_cows(ra,farmj,"s_em")= HA(farmj,"milk") * rebreeding_ratio(farmj);
    gre_dairy_cows(ra,farmj,"s_em") = HA(farmj,"milk") ;
    gre_ote(ra,farmj,"s_em")= ABM(farmj,ra,"OTE");
    gre_age(ra,farmj,"s_em")= ABM(farmj,ra,"AGE");
    gre_uaa(ra,farmj,"s_em")=sum(jV, gre_supv(ra,farmj,jV, "s_em"));
    gre_cuaa(ra,farmj,"s_em")$ (gre_uaa(ra,farmj,"s_em") le 10)                                       = 1;
    gre_cuaa(ra,farmj,"s_em")$((gre_uaa(ra,farmj,"s_em") gt 10)  and (gre_uaa(ra,farmj,"s_em") le 20))  = 2;
    gre_cuaa(ra,farmj,"s_em")$((gre_uaa(ra,farmj,"s_em") gt 20)  and (gre_uaa(ra,farmj,"s_em") le 50))  = 3;
    gre_cuaa(ra,farmj,"s_em")$((gre_uaa(ra,farmj,"s_em") gt 50)  and (gre_uaa(ra,farmj,"s_em") le 100)) = 4;
    gre_cuaa(ra,farmj,"s_em")$((gre_uaa(ra,farmj,"s_em") gt 100) and (gre_uaa(ra,farmj,"s_em") le 300)) = 5;
    gre_cuaa(ra,farmj,"s_em")$ (gre_uaa(ra,farmj,"s_em") gt 300)                                      = 6;
    gre_cuaa(ra,farmj,"s_em")$ (gre_uaa(ra,farmj,"s_em") eq 0)                                        = 0;
    gre_intens(ra,farmj,"s_em")$(Bt2(farmj)AND gre_supTOT(ra,farmj,"s_em")GT 0)=gre_prodv(ra,farmj,"milk", "s_em")*1000/gre_supTOT(ra,farmj,"s_em");
);
DISPLAY$(showdisplay > 0) gre_supv,gre_obj,gre_rented_land;

*-------------------------------------------------------------------------------------------------------
*saving solutions for each simulation both at farm and regional level
re_prezv(proc,t)$(sum((ra,farmj),gre_prodv(ra,farmj,proc,t)) gt 0)
                = sum((ra,farmj),gre_prezv(ra,farmj,proc,t)*gre_prodv(ra,farmj,proc,t))
                / sum((ra,farmj),gre_prodv(ra,farmj,proc,t)) ;
re_prodv(proc,t) = sum((ra,farmj), gre_prodv(ra,farmj,proc,t)) ;
re_prodr(proc,t) = sum((ra,farmj), gre_prodr(ra,farmj,proc,t)) ;
re_prodt(proc,t) = sum((ra,farmj), gre_prodt(ra,farmj,proc,t)) ;
re_supv(proc,t)  = sum((ra,farmj), gre_supv(ra,farmj,proc,t)) ;
re_land_exchange_V(n,m,t) = sum(ra, gre_land_exchange_V(ra,n,m,t));
re_land_exchange_Z(n,m,t) = sum(ra, gre_land_exchange_Z(ra,n,m,t));
re_h2o(proc,t)  = sum((ra,farmj), gre_h2o(ra,farmj,proc,t)) ;
re_co2(proc,t)  = sum((ra,farmj), gre_co2(ra,farmj,proc,t)) ;
re_nitrogen(t) = sum((ra,farmj), gre_nitrogen(ra,farmj,t)); 
re_subpa(t) = sum((ra,farmj), gre_subpa(ra,farmj,t));
re_plvv(t)     = sum((ra,farmj), gre_plvv(ra,farmj,t)) ;
re_plvt(t)     = sum((ra,farmj), gre_plvt(ra,farmj,t)) ;
re_compv(t)    = sum((ra,farmj), gre_compv(ra,farmj,t)) ;
re_aiutidis_GREEN(t)   = sum((ra,farmj), gre_aiutidis_GREEN(ra,farmj,t)) ;
re_aiutidis_BASIC(t)   = sum((ra,farmj), gre_aiutidis_BASIC(ra,farmj,t)) ;
re_aiutidis_PENSION(t)   = sum((ra,farmj), gre_aiutidis_PENSION(ra,farmj,t)) ;
re_aiutiacc(t)   = sum((ra,farmj,j), gre_coupledsubs(ra,farmj,j,t)) ;
re_comptl(t)   = sum((ra,farmj), gre_comptl(ra,farmj,t)) ;
re_comptn(t)   = sum((ra,farmj), gre_comptn(ra,farmj,t)) ;
re_comphv(proc,t)$(sum((ra,farmj),gre_supv(ra,farmj,proc,t)) gt 0)
                = sum((ra,farmj),gre_comphv(ra,farmj,proc,t)*gre_supv(ra,farmj,proc,t))
                / sum((ra,farmj),gre_supv(ra,farmj,proc,t)) ;
re_costv(t)    = sum((ra,farmj), gre_costv(ra,farmj,t)) ;
re_costt(t)    = sum((ra,farmj), gre_costt(ra,farmj,t)) ;
re_obj(t)      = sum((ra,farmj), gre_obj(ra,farmj,t)) ;

*Solutions involving production processes.
set Ky /supT,prezv,prodv,proda,prodt,supv,supb, H2O, CO2, nitrogen, sbpa,
          pombt,plvv,plra,plvt,compv,p_green,p_basic,p_pension,coupled_subs,
          uaa,CUAA,intens,c_intens,age,
          ote,rented_land,
          comptn,comphv,costv,costt,GM,liq,dairy_cows,rebreeding_cows/;
*solutions that do NOT involve production processes (e.g., gross margin).
set Kx /XXX/;

*recall all solutions (Ky and Kx) at farm level level for each simulation.
parameter ris_macroaree(ra,farm,ky,*,t) soluzioni;
  ris_macroaree(ra,farmj,"prezv",j,t)= gre_prezv(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"supT","XXX",t)= gre_supTOT(ra,farmj,t);
  ris_macroaree(ra,farmj,"prodv",j,t)= gre_prodv(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"proda",j,t)= gre_prodr(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"prodt",j,t)= gre_prodt(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"supv"  ,j   ,t) = gre_supv(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"h2o"   ,j   ,t) = gre_h2o(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"co2"   ,j   ,t) = gre_co2(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"nitrogen"   ,"XXX" ,t) = gre_nitrogen(ra,farmj,t);
  ris_macroaree(ra,farmj,"sbpa"  ,"XXX"   ,t) = gre_subpa(ra,farmj,t);
  ris_macroaree(ra,farmj,"plvv","XXX",t)= gre_plvv(ra,farmj,t);
  ris_macroaree(ra,farmj,"plra","XXX",t)= gre_plra(ra,farmj,t);
  ris_macroaree(ra,farmj,"plvt","XXX",t)= gre_plvt(ra,farmj,t);
  ris_macroaree(ra,farmj,"p_green" ,"XXX",t)  = gre_aiutidis_GREEN(ra,farmj,t);
  ris_macroaree(ra,farmj,"p_basic" ,"XXX",t)  = gre_aiutidis_BASIC(ra,farmj,t);
  ris_macroaree(ra,farmj,"p_pension" ,"XXX",t)  = gre_aiutidis_PENSION(ra,farmj,t);
  ris_macroaree(ra,farmj,"coupled_subs" ,j,t)  = gre_coupledsubs(ra,farmj,j,t);
  ris_macroaree(ra,farmj,"costv","XXX",t)= gre_costv(ra,farmj,t);
  ris_macroaree(ra,farmj,"CUAA","XXX",t)= gre_cuaa(ra,farmj,t);
  ris_macroaree(ra,farmj,"uaa","XXX",t)= gre_uaa(ra,farmj,t);
  ris_macroaree(ra,farmj,"intens","XXX",t)= gre_intens(ra,farmj,t);
  ris_macroaree(ra,farmj,"ote","XXX",t)= gre_ote(ra,farmj,t);
  ris_macroaree(ra,farmj,"age","XXX",t)= gre_age(ra,farmj,t);
  ris_macroaree(ra,farmj,"GM","XXX",t)= gre_obj(ra,farmj,t);
  ris_macroaree(ra,farmj,"liq","XXX",t)= gre_fin_liq(ra,farmj,t);
  ris_macroaree(ra,farmj,"rebreeding_cows","XXX",t)= gre_rebreeding_cows(ra,farmj,t);
  ris_macroaree(ra,farmj,"dairy_cows","XXX",t)= gre_dairy_cows(ra,farmj,t);
  ris_macroaree(ra,farmj,"rented_land","XXX",t)= gre_rented_land(ra,farmj,t);
DISPLAY$(showdisplay > 0) re_supv,re_prezv;

*================================================================================================
*                                  P U T                                                        *
*================================================================================================
*creo due file di testo in cui inserire tramite loop le informazioni per ogni scenario sia a livello aziendale che a livello regionale
*modulo_macro a livello aziendale
$include farms_output.prn
*modulo_output soluzioni aggragate a livello regionale
$include regional_output.prn
*$offtext