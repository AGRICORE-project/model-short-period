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

* --------------------------------
* PMP Continuation
* --------------------------------

variables
  uv(farm,j)  INTERCEPT OF ERROR PER farm
  Pqv(j,k) P=Lsqrt(D)
  choleobj
  PDv(j,k)
  PLvv(j,k)
  ;
positive variables
  PDv
  VarC(j) variable esplicit costs per crop to be estimated
  LAMDA(j) variable implicit costs per crop to be estimated
  VARCn(farm,j) variable esplicit costs per farm
  LAMDAn(farm,j) variable implicit costs per farm
  YPSILON(farm) ;
                    
equations
  objequ objective function : error square to be minimized
  definpv(j,jj) function to create a positive simmetric semidefinite matrix (qmat)
  mcv(farm,j) variable marginal costs for realized (observed) crops
  mc1v(farm,j) variable marginal costs for  non-realized crops
  Tcn(farm)  total variable costs per farm (must be equal to Tcostn)
  ADC(farm,j)
  ADC4(farm) cost function: esplicit observed costs + implicit estimated costs > total farm costs (tcostn)
  OPTnV(farm,j) marginal costs = maginal revenue
  FOSTAR(farm) primal obj function = dual obj function
  pesum0v(j) normal distribution of errors
;

 mcv(farmj,j)$xbar(farmj,j)..           varCn(farmj,j)+LAMDAn(farmj,j) =E= uv(farmj,j)+ sum(k, sum(jj, Pqv(j,jj)*Pqv(k,jj))* xbar(farmj,k) ) ;

 mc1v(farmj,j)$(xbar(farmj,j) eq 0)..   varCn(farmj,j)+LAMDAn(farmj,j) =L=  uv(farmj,j)+ sum(k, sum(jj, Pqv(j,jj)*Pqv(k,jj))* xbar(farmj,k) ) ;

 TCn(farmj)..                       sum(j, (varCn(farmj,j))*xbar(farmj,j)) =E= Tcostn(farmj);

 ADC(farmj,j)$(xbar(farmj,j) eq 0)..   varCn(farmj,j) =E= 0;

 ADC4(farmj)..       sum(j,uv(farmj,j)*xbar(farmj,j)) + 0.5 * sum(j,xbar(farmj,j)* sum(k, sum(jj, Pqv(j,jj)*Pqv(k,jj))* xbar(farmj,k) ) )
              =G= TCOSTn(farmj);

 OPTnV(farmj,j)$xbar(farmj,j)..         varCn(farmj,j)+LAMDAn(farmj,j)+YPSILON(farmj)*[ha(farmj,j)/xbar(farmj,j)] =E= (Pr(farmj,j) + sum(coupled_policy, sh(farmj,coupled_policy,j)*[(ha(farmj,j)/xbar(farmj,j))]));

 definpv(j,jj)..                  Pqv(j,jj) =e= sum(k, PLvv(j,k)*((PDv(k,jj)+0.000000001)**0.5) );

 FOSTAR(farmj)..     sum(i, b(farmj,i))*YPSILON(farmj)+sum(j,LAMDAn(farmj,j)*(xbar(farmj,j)+0)) =E= sum(j,PR(farmj,j)*XBAR(farmj,j))+ sum((j,coupled_policy),sh(farmj,coupled_policy,j)*ha(farmj,j))
                                                                                    -sum(j,VARCn(farmj,j)*(XBAR(farmj,j)));

 pesum0v(j)..                     sum(farmj, uv(farmj,j)) =e= 0 ;

 objequ..                        choleobj =e= sum((farmj,j), sqr(uv(farmj,j)))/2 ;
 

*$ontext
$offOrder
 PDv.L(j,k)$(ord(j) eq ord(k))  = .0051 ;
 PDv.FX(j,k)$(ord(j) ne ord(k)) = 0 ;
 PLvv.L(j,k)$(ord(j) ne ord(k)) = .005 ;
 PLvv.FX(j,k)$(ord(j) eq ord(k)) = .02 ;
 PLvv.FX(j,k)$(ord(j) lt ord(k)) = 0. ;
 PDv.LO(j,k)$(ord(j) eq ord(k)) = .001 ;
 uv.l(farmj,j) = 10 ;
$onOrder

 model qcholesky2 /
       mcv,
       mc1v,
       definpv,
       TCn,
*       ADC,
       ADC4,
       OPTnV,
       FOSTAR
       objequ
*       pesum0v
 / ;

 qcholesky2.workspace=5000.0 ;

* qcholesky2.OPTFILE=1 ;
 option  savepoint=1 ;

* execute_loadpoint 'QCHOLESKY2_p.gdx' ;

solve qcholesky2 using NLP minimizing choleobj ;

* recovering of Q matrix of variable costs
parameter qmat(j,k), FU(farm,j), ADCC ;
  qmat(j,k) = sum(jj, Pqv.L(j,jj)*Pqv.L(k,jj) ) ;
  FU(farmj,j)= uv.l(farmj,j);
  ADCC(farmj) = ADC4.L (farmj);

Parameter
ccostq(farm,j) marginal costs
costqa(farm) total costs
var_liq (farm) variation (inizial liquidity - total costs);

Parameter ccostq(farm,j), costqa(farm);
ccostq(farmj,j)= fu(farmj,j) + sum(k, qmat(j,k)*xbar(farmj,k));
costqa(farmj)= sum(j,fu(farmj,j) * xbar(farmj,j)) + .5* sum((j,k), xbar(farmj,j) * qmat(j,k) * xbar(farmj,k) );
var_liq(farmj)= initial_liq(farmj)- ADCC(farmj);

DISPLAY$(showdisplay > 0) qmat,ccostq,costqa, var_liq;

*non linear model calibration
parameter uu(j), nsolx(farm,j), lpxh(j), nsoly(farm,i), nlobj(farm), nlsolshx(farm,j),nlshx(j),
           pv(j),shv(j) ;

VARIABLES  nX(j)     output in activity units
            nxH(j)    output of compensated surface
            NLINPROF   LP PROFIT  ;

POSITIVE VARIABLE nX, NXH,Ndquota;

EQUATIONS RESOURCEN(i)   CONSTRAINED RESOURCES
           NACCOUNT(j,i)    BILANCING LAND
           NLPROFIT       LP OBJECTIVE FUNCTION
           nACCOUNTr(j)  Reuse RICErse
;

 RESOURCEN(i)..           SUM(jv,AA(jv,i)*nX(jv))   =L= BB(i) ;
 NACCOUNT(j,i)..            (AA(j,i)*nX(j))-nXH(j) =E=0;
 nACCOUNTr(jr)..             (AAl(jr)*nX("milk"))-nX(jr) =L=0;

NLPROFIT..     SUM(j, nX(j)* Pv(j) )
             + SUM(j, nXh(j)* SSH(j))
             + qqp
             - sum(j,uu(j) * nx(j))
             - .5* sum((j,k), nx(j) * qmat(j,k) * nx(k) )
             =E= NLINPROF;

MODEL  PRIMAL /RESOURCEN, NACCOUNT, NLPROFIT, NACCOUNTR/;

* LOOP
lpsolx(farmj,j)=xbar(farmj,j);
lpsolshx(farmj,j)=ha(farmj,j);

LOOP((farmj)$b(farmj,"land"),

AA(j,i) = A(farmj,j,i) ;
AAL(jr)=al(farmj,jr);
DDA(j,jj,i)= DA(farmj,j,jj,i);
BB(i) = B(farmj,i) ;
Pv(j) = PR(farmj,j) ;
SSH(j) = sum(coupled_policy, SH(farmj,coupled_policy, j));
qqp=spp(farmj);
UU(j)  = fu(farmj,j) ;

nxh.l(j) = ha(farmj,j) ;
nX.L(j) = xbar(farmj,j) ;

SOLVE primal  USING NLP MAXIMIZING NLINPROF ;
nsolx(farmj,j) = nx.l(j) ;
nlsolshx(farmj,j) = nxh.l(j) ;
nsoly(farmj,i) = resourcen.m(i)  ;
nlobj(farmj) = nlinprof.l ;

);
DISPLAY$(showdisplay > 0) ha,nlsolshx,xbar,nsolx,nsoly;

*CALIBRATION

*calibration control points:
parameter contcal1(farm,j);
contcal1(farmj,j)$lpsolx(farmj,j)=(nsolx(farmj,j)/lpsolx(farmj,j)*100)-100;
DISPLAY$(showdisplay > 0) contcal1;
parameter contcal2(farm,j);
contcal2(farmj,j)$xbar(farmj,j)=(lpsolx(farmj,j)/xbar(farmj,j)*100)-100;
DISPLAY$(showdisplay > 0) contcal2;

Parameter supco0(j),supco1(j),supco2(j);
supco0(j)=sum(farmj,ha(farmj,j));
supco1(j)=sum(farmj,lpsolshx(farmj,j));
supco2(j)=sum(farmj,nlsolshx(farmj,j));
DISPLAY$(showdisplay > 0) supco0,supco1,supco2;

*---------------------------------------------------------------------------------------------------------------------
* exporting calibration results
*--------------------------------------------------------------------------------------------------------------------
*                                                           exporting calibration results

* IF platform is linux
if (platform <= 0,
  execute_unload "calibration_results.gdx" qmat;
  execute 'gdxdump calibration_results.gdx Symb=qmat Output=qmat.csv NoHeader Format=csv';
  execute_unload "calibration_results.gdx" fu;
  execute 'gdxdump calibration_results.gdx Symb=fu Output=fu.csv NoHeader Format=csv';
*nsolshx= ha
  execute_unload "calibration_results.gdx" nlsolshx;
  execute 'gdxdump calibration_results.gdx Symb=nlsolshx Output=nlsolshx.csv NoHeader Format=csv';
*nsolx =xbar
  execute_unload "calibration_results.gdx" nsolx
  execute 'gdxdump calibration_results.gdx Symb=nsolx Output=nsolx.csv NoHeader Format=csv';
  execute_unload "calibration_results.gdx" ccostq
  execute 'gdxdump calibration_results.gdx Symb=ccostq Output=ccostq.csv NoHeader Format=csv';
  execute_unload "calibration_results.gdx" costqa
  execute 'gdxdump calibration_results.gdx Symb=costqa Output=costqa.csv NoHeader Format=csv';
  execute_unload "calibration_results.gdx" nlobj
  execute 'gdxdump calibration_results.gdx Symb=nlobj Output=nlobj.csv NoHeader Format=csv';
  execute_unload "calibration_results.gdx" nsoly
  execute 'gdxdump calibration_results.gdx Symb=nsoly Output=nsoly.csv NoHeader Format=csv';
*ADCC = TOTAL COST FUNCTION (ESPLICIT + IMPLICIT)
  execute_unload "calibration_results.gdx" ADCC
  execute 'gdxdump calibration_results.gdx Symb=ADCC Output=ADCC.csv NoHeader Format=csv';
else
  execute_unload "calibration_results.gdx" qmat
  execute 'gdxviewer.exe i=calibration_results.gdx csv=qmat.csv id=qmat';
  execute_unload "calibration_results.gdx" fu
  execute 'gdxviewer.exe i=calibration_results.gdx csv=fu.csv id=fu';
*nsolshx= ha
  execute_unload "calibration_results.gdx" nlsolshx
  execute 'gdxviewer.exe i=calibration_results.gdx csv=nlsolshx.csv id=nlsolshx';
*nsolx =xbar
  execute_unload "calibration_results.gdx" nsolx
  execute 'gdxviewer.exe i=calibration_results.gdx csv=nsolx.csv id=nsolx';
  execute_unload "calibration_results.gdx" ccostq
  execute 'gdxviewer.exe i=calibration_results.gdx csv=ccostq.csv id=ccostq';
  execute_unload "calibration_results.gdx" costqa
  execute 'gdxviewer.exe i=calibration_results.gdx csv=costqa.csv id=costqa';
  execute_unload "calibration_results.gdx" nlobj
  execute 'gdxviewer.exe i=calibration_results.gdx csv=nlobj.csv id=nlobj';
  execute_unload "calibration_results.gdx" nsoly
  execute 'gdxviewer.exe i=calibration_results.gdx csv=nsoly.csv id=nsoly';
*ADCC = TOTAL COST FUNCTION (ESPLICIT + IMPLICIT)
  execute_unload "calibration_results.gdx" ADCC
  execute 'gdxviewer.exe i=calibration_results.gdx csv=ADCC.csv id=ADCC';
) ;





$exit;
