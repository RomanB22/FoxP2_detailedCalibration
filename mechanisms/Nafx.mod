: Fast Na+ channel
: added the 's' attenuation system from hha2.mod
: Kiki Sidiropoulou
: September 27, 2007

NEURON {
	SUFFIX Nafx
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina, gna, ar2, scale_gnaf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
	gnafbar	= 0 (mho/cm2)
	scale_gnaf = 1
	:gnafbar= 0.086 (mho/cm2) <0,1e9>
	ena = 55 (mV)
	
	:PARAMETERS FOR S ATTENUATION SYSTEM
	taumin = 30 (ms)  :min activation time for "s" attenuation system 30
        vhalfr =-60 (mV)       :half potential for "s" attenuation system, -60
        vvh=-58		(mV)  : -58
 	vvs = 2 (mV)
	a0r = 0.0003 (/ms)
        b0r = 0.0003 (/ms)
       : a0r = 0.0003 (ms)
        :b0r = 0.0003 (ms)
        zetar = 12    
	zetas = 12   
        gmr = 0.2   
	ar2 = 1.0               :initialized parameter for location-dependent 1.0
                                :Na-conductance attenuation, "s", (ar=1 -> zero attenuation)
}
STATE {
	m h s
}
ASSIGNED {
	celsius (degC)
	ina (mA/cm2)
	minf 
	hinf
	sinf 
	mtau (ms)
	htau (ms)
	stau (ms)
	gna (mho/cm2)
	
}



INITIAL {
	rate(v, ar2)
	m = minf
	h = hinf
	s = sinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = scale_gnaf*gnafbar*m*m*m*h*s
	ina = gna*(v-55)
	
}

DERIVATIVE states {
	rate(v, ar2)
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
	s' = (sinf-s)/stau
}

UNITSOFF

FUNCTION malf( v){ LOCAL va, varAux1, varAux2, varAux3
	varAux1 = 28 : 28
	varAux2 = 9.3 : 9.3
	varAux3 = -0.2816 : -0.2816
	va=v+varAux1 : Relevant parameter
	:va=v+28
	if (fabs(va)<1e-04){
	   malf= varAux3*(-varAux2 + va*0.5)
	   :malf= varAux3*(-varAux2 + va*0.5)
	}else{
	   malf = varAux3*(v+varAux1)/(-1+exp(-va/varAux2))
	}
}


FUNCTION mbet(v(mV))(/ms) { LOCAL vb, varAux1, varAux2, varAux3
	varAux1 = 1 : 1
	varAux2 = 6 : 6
	varAux3 = 0.24 : 0.2464
	vb=v+varAux1
	:vb=v+1
	if (fabs(vb)<1e-04){
	    mbet = varAux3*(varAux2+vb*0.5)
	    :mbet = 0.2464*(varAux2 + vb*0.5)
	}else{
	   mbet = varAux3*(vb)/(-1+exp(vb/varAux2))	  :/(-1+exp((v+1)/6))
	}
	}	


FUNCTION half(v(mV))(/ms) { LOCAL vc, varAux1, varAux2, varAux3
	varAux1 = 40 : 40 : Does not make a big change in fI
	varAux2 = 20 : 20
	varAux3 = 0.098 : 0.098
	vc=v+varAux1
	:vc=v+15.1	:changed to 40.1 by kiki
	if (fabs(vc)<1e-04){
	   half=varAux3*(varAux2 + vc*0.5)
	}else{
	   half=varAux3/exp(vc+varAux1/varAux2)  :43.1, also spike train attenuation
}
}


FUNCTION hbet(v(mV))(/ms) { LOCAL vd, varAux1, varAux2, varAux3
	varAux1 = 13 : 13 : Not change
	varAux2 = 10 : 10
	varAux3 = 1.4 : 1.4
	vd=v+varAux1
	:vd=v+13.1  :decreasing it increases the peak current
	if (fabs(vd)<1e-04){
	   hbet=varAux3*(varAux2 + vd*0.5)
	}else{
	   hbet=varAux3/(1+exp(-(vd-varAux1)/varAux2))  :13.1 increasing it, increases the spike train attenuation and increases spike width
} 
}


:FUNCTIONS FOR S 
FUNCTION alpv(v(mV)) {
         alpv = 1/(1+exp((v-vvh)/vvs))
}


FUNCTION alpr(v(mV)) {       :used in "s" activation system tau

  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betr(v(mV)) {       :used in "s" activation system tau

  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}



PROCEDURE rate(v (mV),ar2) {LOCAL q10, msum, hsum, ma, mb, ha, hb,c
	

	ma=malf(v) mb=mbet(v) ha=half(v) hb=hbet(v)
	
	msum = ma+mb
	minf = ma/msum
	mtau = 1/(msum)
	
	
	hsum=ha+hb
	hinf=ha/hsum
	htau = 1 / (hsum)

	stau = betr(v)/(a0r*(1+alpr(v))) 
	if (stau<taumin) {stau=taumin} :s activation tau
	c = alpv(v)
	sinf = c+ar2*(1-c) 	
}

	
UNITSON


