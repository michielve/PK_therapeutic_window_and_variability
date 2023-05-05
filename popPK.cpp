$PARAM @annotated
TVKA : 1 : Absorption rate constant (/time)
TVCL : 2 : Clearance (L/time)
TVVC : 10 : Central V (L)
TVBASE : 100 : Effect Baseline
TVSLOPE : 10 : Slope effect


$CMT GUT CENT

$MAIN
double KA = TVKA;
double CL = TVCL*exp(ETA(1));
double VC = TVVC*exp(ETA(2));

double BASE = TVBASE;
double SLOPE = TVSLOPE;


double k10 = CL/VC;

$OMEGA 0 0
$SIGMA @labels ADD
0

$ODE
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - k10*CENT;



$TABLE
// Set PK concentration
double CONCENTRATION = CENT/VC;

// SET effect formula
double EFFECT = (BASE + SLOPE*CONCENTRATION)+ADD;


while(EFFECT < 0) {
	simeps();
	EFFECT = (BASE + SLOPE*CONCENTRATION)+ADD;
}


$CAPTURE CONCENTRATION EFFECT
