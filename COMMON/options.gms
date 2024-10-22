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

* Define 1 for enable DISPLAY 0 for disable it
scalar showdisplay / 0 /;

* Define 0 for linux, 1 for windows
scalar platform / 0 /;

$offListing;

OPTION solprint   = off;