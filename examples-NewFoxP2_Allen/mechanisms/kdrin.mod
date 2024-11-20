: Delayed rectifier K+ channel

NEURON {
	SUFFIX kdrin
	USEION k READ ki, ko WRITE ik
	RANGE gbar, ik, gk, scale_gkdr
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gbar= 0.0338 (mho/cm2) <0,1e9>
	scale_gkdr = 1
	
	
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
	inf
	tau (ms)
	gk (mho/cm2)
	ek (mV)
	ki (mM)
	ko (mM)

}


INITIAL {
	rate(v)
	n = inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk= scale_gkdr*gbar*n*n*n*n
	ek = 25 * log(ko/ki)
	ik = gk*(v-ek)
	
}

DERIVATIVE states {
	rate(v)
	n' = (inf-n)/tau
}

UNITSOFF

FUNCTION alf(v){ LOCAL va, varAux1, varAux2, varAux3
	varAux1=13 : 13
	varAux2=25 : 25
	varAux3=-0.018 : -0.018
		va=v-varAux1
	   :va=v-13
	if (fabs(va)<1e-04){
	   va=va+0.0001
		alf= (varAux3*va)/(-1+exp(-(va/varAux2)))
	} else {
	  	alf = (varAux3*(v-varAux1))/(-1+exp(-((v-varAux1)/varAux2)))
	}
}


FUNCTION bet(v) { LOCAL vb, varAux1, varAux2, varAux3
	varAux1=23 : 23
	varAux2=12 : 12
	varAux3=0.0054 : 0.0054
	vb=v-varAux1
	  :vb=v-23
	if (fabs(vb)<1e-04){
	  vb=vb+0.0001
		bet= (varAux3*vb)/(-1+exp(vb/varAux2))
	} else {
	  	bet = (varAux3*(v-varAux1))/(-1+exp((v-varAux1)/varAux2))
	}
}	






PROCEDURE rate(v (mV)) {LOCAL q10, sum, aa, ab
	
	aa=alf(v) ab=bet(v) 
	
	sum = aa+ab
	inf = aa/sum
	tau = 1/(sum)
	
	
}

UNITSON	



