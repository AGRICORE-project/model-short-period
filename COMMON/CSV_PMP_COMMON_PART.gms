* ------------------------------------------------------------------------------------------------------------------
* PMP CALIBRATION MODULE
* ---------------------------------------------------------------------------------------------------------------

*setting probabilities for qmat construction
set pp  probability
/1*5/

*limiting factor =land
set i /land/

*creating ALIAS es. j=jj=k=kk) to create square matrix
 ALIAS (j,jj,k,kk)
 ALIAS (i,ii)
 Alias (pp,l,ll)

 PARAMETER AA(j,i) , DDA(j,jj,i), BB(i), CC(j) , PPR(j) , XXBAR(j),
           COnT(j), AVEC2(j), lambda(farm,j), lpsolx(farm,j),
           lpsolshx(farm,j),
           lpsoly(farm,i), lpobj(farm), SSH(j), spp(farm), SH(farm,coupled_policy,proc),
           ha(farm,proc), xbar(farm,proc),pr(farm,proc),c(farm,proc), qqp;


*setting scalar epsilon to avoid linear dependence
 SCALAR EPSILON / 0.000000001/ ;
ha(farmj,j)=hax(farmj,j);
xbar(farmj,j)=xbarx(farmj,j);
pr(farmj,j)=prx(farmj,j);
c(farmj,j)=cx(farmj,j);

Parameter h2otot(j) water consumption;
h2otot(j)=sum((farmj,water),xbar(farmj,j)*h20(j,water))/10000;
DISPLAY$(showdisplay > 0) h2otot;

*single farm payment=spp = decoupled payments at farm level + retirement pension if older than 65 and no children
parameter pens(farm);
pens(farmj)$(selage(farmj) gt 65 and selsucc(farmj) eq 0)=12000;
pens(farmj)$(selage(farmj) le 65 or(selage(farmj) gt 65 and selsucc(farmj) ne 0))=0;

spp(farmj)=susx(farmj,"P_GREEN")+susx(farmj,"P_BASIC")+ pens(farmj) ;

parameter hac(farm,coupled_policy,proc) ha of crops receiving subsidies;
hac(farmj,coupled_policy,proc)$(COUPLED_SUBSIDIES_RATES(coupled_policy,proc))=hax(farmj,proc);
DISPLAY$(showdisplay > 0) hac;

*from subsidies per unit of product to subsidies per hectar
sh(farmj,coupled_policy,j)=COUPLED_SUBSIDIES_RATES(coupled_policy,j);
DISPLAY$(showdisplay > 0) spp, hac, sh;

*SCALING: value downsizing by 10000
xbar(farmj,proc) = xbar(farmj,proc)/1000;
pr(farmj,proc) = pr(farmj,proc)/10;
sh(farmj,coupled_policy,proc) = sh(farmj,coupled_policy,proc)/10000;
c(farmj,proc) = c(farmj,proc)/10;
spp(farmj)=spp(farmj)/10000;
initial_liq(farmj)=initial_liq(farmj)/10000;
rental_price(farmj)=rental_price(farmj)/10000;

*LIVESTOCK CONSTRAINTS
 Parameter T_LSU(farm) totale dairy_cows,cmr(farm,proc) forage reuse per dairy_cow, fr(farm,proc) total forage reuse, al(farm,proc) reverse yield forage reuse, aal(proc);

 T_LSU(farmj) = ha(farmj,"milk");
 CMR(farmj,jr)$T_LSU(farmj) = XBAR(farmj,jr)/T_LSU(farmj);
 FR(farmj,jr) = CMR(farmj,jr) * T_LSU(farmj);
 AL(farmj,jr) $ XBAR(farmj,"milk") = FR(farmj,jr) / XBAR(farmj,"milk") ;
 
DISPLAY$(showdisplay > 0) T_LSU, cmr, fr, al;

*sellink prices of forages is 0 (only for dairy farms) as they are ALL reused for animal feeding
pr(farmj,jr)$ha(farmj,"milk")=0;
 DISPLAY$(showdisplay > 0) pr;

*Difc: to identify companies where the price is lower than the unit cost
parameter difc(farm,j);
difc(farmj,j)=pr(farmj,j)-c(farmj,j);
c(farmj,j)$(difc(farmj,j) le 0)=pr(farmj,j)*0.95;
c(farmj,j)$(ha(farmj,j) eq 0)=0;


* defining matrixes of technical coeff and available resources (land)

* A matrix of the inverse yield of the land: how much land does it take to get one unit of the product?
* B: sum of UAA = total area
* DA: diagonal matrix

PARAMETER A(farm,j,I) , B(farm,I), DA(farm,j,jj,I) ;
A(farmj,j,I) $ XBAR(farmj,j) = HA(farmj,j)/XBAR(farmj,j) ;

*EXCLUDING MILK from the summation of cultivation processes!
B(farmj,I) = SUM(j$(not sameas (j,"milk")),(HA(farmj,j)));
DA(farmj,j,jj,i) =0;
* $offOrder
DA(farmj,j,jj,i)$(ORD(j) EQ ORD (jj)) = A(farmj,j,i)
* $onOrder
DISPLAY$(showdisplay > 0) A, DA , B, pr, xbar, c, initial_liq;

PARAMETER TCOSTn(farm);
TCOSTn(farmj)=sum(j, c(farmj,j)*xbar(farmj,j));