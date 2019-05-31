/***************************************************
Jeff Fox's canine ventricular cell ionic model 
(J. J. Fox, J. L. McHarg, and R. F. Gilmour Jr,
Ionic mechanism of electrical alternans,
Am J Physiol Heart Circ Physiol 2002; 282: 516-530) 

Incorporating Timothy Syndrome mutation by deleting the f gate of L-type Ca2+ channel

Zheng  April 2005

*****************************************************/

#include <stdio.h>
#include <iostream>
#include <math.h>

#define NAME 20

using namespace std;


	const double RT_on_F = 26.7;
	const double F = 96.5;
	const double A_cap = 1.534E-4;
	const double C_sc = 1.0;
	const double V_myo = 25.84E-6;
	const double V_sr = 2.0E-6;

	const double G_Na = 12.8;
	const double G_K1 = 2.8;
	const double G_Kr = 0.0136;
	const double G_Ks = 0.0245;
	const double G_to = 0.23815;
	const double G_Kp = 0.002216;
	const double G_Nab = 0.0031;
	const double G_Cab = 0.0003842;
	const double P_bar_Ca = 0.0000226;
	const double P_bar_CaK = 5.79E-7;
	const double P_rel = 6.0;
	const double G_leak = 0.000001;
	const double I_bar_NaK = 0.693;
	const double I_Cahalf = -0.265;
	const double I_bar_pCa = 0.05;
	const double eta = 0.35;
	const double ksat = 0.2;
	const double K_NaCa = 1500.0;
	const double Km_fCa = 0.18;
	const double Km_K1 = 13.0;
	const double Km_Na = 87.5;
	const double Km_Ca = 1380.0;
	const double Km_Nai = 10.0;
	const double Km_Ko = 1.5;
	const double Km_pCa = 0.05;
	const double Km_up = 0.32;
	const double CMDN_tot = 10.0;
	const double CSQN_tot = 10000.0;
	const double Km_cmdn = 2.0;
	const double Km_csqn = 600.0;
	const double V_up = 0.1;

	const double Na_in = 10.0;
	const double K_in = 149.4;
	const double Na_out = 138.0;
	const double K_out = 4.0;
	const double Ca_out = 2000.0;

		
	double E_Na, E_K, E_Ks, E_Ca;
	
	double Ca_in, Ca_sr;
	
	double I_Na, m, h, j, alpha_m, beta_m, alpha_h, beta_h, alpha_j, beta_j, sm, sh, sj, taum, th, tj;
    double I_Ca_L, d, f, fca, I_Cal, I_Ca_K, I_Ca_full, sd, td, sf, tf, sfca, tfca;
	double I_Kr, xr, r, sxr, txr;
	double I_Ks, xs, sxs, txs;
    double I_K1, sK1;
	double I_to, xto,  yto, axto, bxto, ayto, byto, sxto, txto, syto, tyto;
    double I_Kp, Kp;
	double I_Na_Ca_Exch, I_Na_K_pump, sigma, f_NaK, I_ns_Ca, I_pCa;
	double I_Cab, I_Nab;

	double I_ns_K, I_ns_Na;
	
	double I_Total;
	
	double I_leak, I_up, gama, I_rel, beta_sr, beta_in;
	
	double V, dV;
	double dt;
		
	int counter, nn;
	int stimuli, number, cycle_length;
	double t, tt;
	double I_inj;
	
	int step, tstep;
	
	double Calcu_I_Total( );
	double Calcu_I_Na( );
	double Calcu_I_Ca_L( );
	double Calcu_I_Ks( );
	double Calcu_I_Kr( );
	double Calcu_I_K1( );
	double Calcu_I_to( );
	double Calcu_I_Kp( );
	double Calcu_I_Na_Ca_Exch( );
	double Calcu_I_Na_K_pump( );
	double Calcu_I_pCa( );
	double Calcu_I_Cab( );
	double Calcu_I_Nab( );
	double Calcu_I_leak( );
	double Calcu_I_up( );
	double Calcu_I_rel ( );
	void calcium_update( );
		
int main (){
	
	int fold;
	
	FILE *output;
	char output_file_name[ NAME ];
	printf ("\noutput file name (something.dat):" );
	scanf("%s", output_file_name);
	output=fopen(output_file_name, "w");
		
	dt=0.0025;

	fold=1.0/dt;
	
	stimuli=100;
	cycle_length=500;
	
	counter=0;
	
	V = -94.2976747654229399;
	Ca_in = 0.0534618941422192;
	Ca_sr = 322.5921809971255243;
	d = 0.0000013588414970;
	m = 0.0003221097821464;
	h = 0.9983653992115390;
	j = 0.9947244678629071;
	fca = 0.8369621264104180;
	xr = 0.4117860639678651;
	xs = 0.0159624147222210;
	xto = 0.0000436935755321;
	yto = 0.9998244518242049;

			
	E_K=RT_on_F*log(K_out/K_in);
	E_Na=RT_on_F*log(Na_out/Na_in);
	E_Ks=RT_on_F*log((K_out+0.01833*Na_out)/(K_in+0.01833*Na_in));
		
		
    dV = 0.0; 
	tt=0;
	
	number=0;
	
	while (number<=stimuli){
	
	if (number==0) {step=10000*fold;}
	else {step=cycle_length*fold;}
	
	t=0.0;
	
	for (tstep=0; tstep<=step; tstep+=1)
	{
		if(number >0 && t<=1.0) {I_inj=-80.0;}
		else {I_inj=0.0;}
		
		E_Ca = (RT_on_F/2.0)*log(Ca_out/Ca_in);
		
		I_Total=Calcu_I_Total( );
		
		dV = -1*I_Total*dt;
		V=V+dV;
		
		calcium_update ( );
		
		
		if( number> (stimuli-2) &&  (counter%(fold*5)==0 || (fabs(I_Na)>1 && fabs(I_Na)<280  && counter%(fold/4)==0) || fabs(I_Na)>280 ) )
		
		fprintf(output, "%6.4f\t%6.2f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", tt, V, I_Total, I_Na, I_Cal, I_Ks, I_Kr, I_K1, I_to, Ca_in, Ca_sr, I_Na_K_pump, I_Na_Ca_Exch, I_Ca_K);			
		
		
		counter+=1;
		t+=dt;
		tt+=dt;
		
		}
		number+=1;
		cout << "now the number is: " << number << endl;
	}
	
	fclose ( output );
    return 0;
}

double Calcu_I_Total( ){
	
	I_Na=Calcu_I_Na( );
	I_K1=Calcu_I_K1( );
	I_Kr=Calcu_I_Kr( );
	I_Ks=Calcu_I_Ks( );
	I_to=Calcu_I_to( );
	I_Kp=Calcu_I_Kp( );
	I_Na_K_pump=Calcu_I_Na_K_pump( );
	I_Nab=Calcu_I_Nab( );
	I_Na_Ca_Exch=Calcu_I_Na_Ca_Exch( );
	I_pCa=Calcu_I_pCa( );
	I_Cab=Calcu_I_Cab( );
	I_Ca_L=Calcu_I_Ca_L( );
	
	I_rel= Calcu_I_rel( );
	I_up=Calcu_I_up( );
	I_leak=Calcu_I_leak( );
	
	I_Total = I_Na+I_Nab+I_Na_Ca_Exch+I_Na_K_pump +I_Cal+I_Cab+I_pCa + I_Ks+I_Kr+I_K1+I_to+I_Ca_K+I_Kp+I_inj;
	
	return I_Total;
}



// fast Na current
double Calcu_I_Na( ){
    
			
	alpha_m = 0.32*(V+47.13)/(1.0-exp(-0.1*(V+47.13)));
	beta_m = 0.08*exp(-V/11.0);
    
	alpha_h = 0.135*exp(-(V+80.0)/6.8);
	beta_h = 7.5/(1.0+exp(-0.1*(V+11.0)));
	alpha_j = (0.175*exp(-(V+100.0)/23.0))/(1.0+exp(0.15*(V+79.0)));
	beta_j = 0.3/(1.0+exp(-0.1*(V+32.0)));
	   	
    sm=alpha_m/(alpha_m+beta_m);
	sh=alpha_h/(alpha_h+beta_h);
	sj=alpha_j/(alpha_j+beta_j);
	
	taum=1/(alpha_m+beta_m);
	th=1/(alpha_h+beta_h);
	tj=1/(alpha_j+beta_j);
	
	
	m=sm-(sm-m)*exp(-dt/taum);
	h=sh-(sh-h)*exp(-dt/th);
	j=sj-(sj-j)*exp(-dt/tj);
		
	I_Na = G_Na*(V-E_Na)*m*m*m*h*j;
	
	return I_Na;
}


// currents through L type Ca channel
double Calcu_I_Ca_L( ){

   	I_Ca_full=(P_bar_Ca/C_sc)*(4.0*V*F/RT_on_F)*((Ca_in*exp(2.0*V/RT_on_F)-0.341*Ca_out)/(exp(2.0*V/RT_on_F)-1.0));
		
	sd=1.0/(1.0+exp(-(V+10.0)/6.24));
	td=1.0/((0.25*exp(-0.01*V)/(1.0+exp(-0.07*V)))+(0.07*exp(-0.05*(V+40.0))/(1.0+exp(0.05*(V+40.0)))));
	
	sfca=1.0/(1.0+pow((Ca_in/0.18),3.0));
	tfca=30.0;
	
	d=sd-(sd-d)*exp(-dt/td);
	fca=sfca-(sfca-fca)*exp(-dt/tfca);
		
	I_Cal=d*fca*I_Ca_full;
	
	I_Ca_K=(P_bar_CaK/C_sc)*(d*fca/(1.0+(I_Ca_full/I_Cahalf)))*(1000.0*V*F/RT_on_F)*((K_in*exp(V/RT_on_F)-K_out)/(exp(V/RT_on_F)-1.0));
	
	I_Ca_L=I_Cal+I_Ca_K;
	
	return I_Ca_L;
}

//rapidly activating potassium current
double Calcu_I_Kr () {

	sxr=1.0/(1.0+exp(-2.182-0.1819*V));
	txr=43.0+(1.0/(exp(-5.495+0.1691*V)+exp(-7.677-0.0128*V)));
	
	xr=sxr-(sxr-xr)*exp(-dt/txr);
			
	r=1.0/(1.0+2.5*exp(0.1*(V+28.0)));
	
	I_Kr=G_Kr*xr*r*sqrt(K_out/4.0)*(V-E_K);
	return I_Kr;
}


	
//slowly activating K current
double Calcu_I_Ks () {

	sxs=1.0/(1.0+exp(-(V-16.0)/13.6));
	txs=1.0/((0.0000719*(V-10.0)/(1.0-exp(-0.148*(V-10.0))))+(0.000131*(V-10.0)/(exp(0.0687*(V-10.0))-1.0)));
		
	xs=sxs-(sxs-xs)*exp(-dt/txs);
				
	I_Ks=G_Ks*xs*xs*(V-E_Ks);
	return I_Ks;
}

  
// time-independent K current
double Calcu_I_K1( ){
    
   	sK1=1.0/(2.0+exp((1.62/RT_on_F)*(V-E_K)));
    
    I_K1=G_K1*sK1*(K_out/(K_out+Km_K1))*(V-E_K); 
    return I_K1;
}

//transient outward K current
double Calcu_I_to( ) {
	
	
	axto=0.04516*exp(0.03577*V);
	bxto=0.0989*exp(-0.06237*V);
	ayto=(0.005415*exp(-(V+33.5)/5.0))/(1.0+0.051335*exp(-(V+33.5)/5.0));
	byto=(0.005415*exp((V+33.5)/5.0))/(1.0+0.051335*exp((V+33.5)/5.0));
	
	sxto=axto/(axto+bxto);
	txto=1/(axto+bxto);
	syto=ayto/(ayto+byto);
	tyto=1/(ayto+byto);
	
	xto=sxto-(sxto-xto)*exp(-dt/txto);
	yto=syto-(syto-yto)*exp(-dt/tyto);
		
	I_to=G_to*xto*yto*(V-E_K);
	return I_to;
}


// plateau K current
double Calcu_I_Kp( ){
    
	Kp=1.0/(1.0+exp((7.488-V)/5.98));
    I_Kp=G_Kp*Kp*(V-E_K);
	return I_Kp;
}

// Na-Ca exchanger current
double Calcu_I_Na_Ca_Exch( ) {
	
	I_Na_Ca_Exch=(K_NaCa/(pow(Km_Na,3.0)+pow(Na_out,3.0)))*(1.0/(Km_Ca+Ca_out))*(1.0/(1.0+ksat*exp(V*(eta-1.0)/RT_on_F)))*((pow(Na_in,3.0)*Ca_out*exp(V*eta/RT_on_F))-(pow(Na_out,3.0)*Ca_in*exp(V*(eta-1.0)/RT_on_F)));			
	return I_Na_Ca_Exch;
}

// Na_K pump current
double Calcu_I_Na_K_pump( ){
    
	sigma=(1.0/7.0)*(exp(Na_out/67.3)-1.0);
	f_NaK=1.0/(1.0+0.1245*exp(-0.1*V/RT_on_F)+0.0365*sigma*exp(-V/RT_on_F));
	I_Na_K_pump=I_bar_NaK*f_NaK*(1/(1+pow((Km_Nai/Na_in), 1.5)))*(K_out/(K_out+Km_Ko));
	return I_Na_K_pump;
}

// sarcilemmal Ca pump
double Calcu_I_pCa( ){
	I_pCa=I_bar_pCa*(Ca_in/(Km_pCa+Ca_in));
	return I_pCa;
}

    
//Ca background current
double Calcu_I_Cab( ){
   
   I_Cab=G_Cab*(V-E_Ca);
   return I_Cab;
}

// Na background current
double Calcu_I_Nab ( ){

    I_Nab=G_Nab*(V-E_Na);
	
	return I_Nab;
}

//calcium flux
double Calcu_I_leak( ){
	I_leak=G_leak*(Ca_sr-Ca_in);
	return I_leak;
}

double Calcu_I_up( ){
	//I_up=V_up*Ca_in*Ca_in/(Ca_in*Ca_in+Km_up*Km_up);
	I_up=V_up/(1.0+pow((Km_up/Ca_in),2.0));
	return I_up;
}


double Calcu_I_rel( ) {
	gama=1.0/(1.0+pow((2000.0/Ca_sr), 3));
	I_rel=P_rel*d*fca*(gama*Ca_sr-Ca_in)/(1.0+1.65*exp(V/20.0));
	return I_rel;
}


//dynamic changes of Ca in myoplasm and SR
void calcium_update  ( ) {
	
	beta_sr=1.0/(1.0+CSQN_tot*Km_csqn/pow((Km_csqn+Ca_sr), 2));
	Ca_sr=Ca_sr+dt*(beta_sr*(I_up-I_leak-I_rel)*V_myo/V_sr);
		
	beta_in=1.0/(1.0+CMDN_tot*Km_cmdn/pow((Km_cmdn+Ca_in), 2));
	Ca_in=Ca_in+dt*beta_in*(I_rel+I_leak-I_up-(A_cap/(2*F*V_myo))*(I_Cal+I_Cab+I_pCa-2*I_Na_Ca_Exch));;
	}


