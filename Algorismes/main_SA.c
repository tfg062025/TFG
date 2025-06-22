#include <stdio.h>
#include <stdlib.h>
#include "FMD.h"
#include <float.h>
#include <string.h>
#include <time.h>

/********** HIPERPARAMENTRES GA ***********************************************

* Lmax   : Màxim nombre de proves a cada temperatura
* Lamax  : Nombre màxim de proves acceptades a cada temperatura
* HTsw   : Indicador de temperatura alta 
* LTsw   : Indicador de temperatura baixa
* alp    : Factor de decerixament de temperatura
* zerot  : Indicador de congelament del sistema
* tini   : Temperatura inicial de procés de recalentat
* dt     : Increment de temperatura del procés de recalentat                                                     
* KB     : Constant de Boltzmann 

*******************************************************************************/

#define Lmax  2000
#define Lamax 200
#define HTsw  0.8
#define LTsw  0.000001
#define alp   0.99
#define zerot 1.0e-12
#define tini  1000.0
#define KB    1.0

/********** DEFINICIÓ LIMITS COMPACTE K ***************************************/

#define x_0_min 12726.0
#define x_0_max 17932.0
#define phi_min 0.12
#define phi_max 0.3489494104672237
#define lambda_min 300.0
#define lambda_max 3000.0
#define mu_min 0.0
#define mu_max 10.0
#define sigma_min 0.0
#define sigma_max 50.0
#define delta_min 0.0
#define delta_max 20000.0

/********** RANG DELS PARAMETRES ***********************************************/

#define rang_x_0        (x_0_max - x_0_min)
#define rang_phi        (phi_max - phi_min)
#define rang_lambda     (lambda_max - lambda_min)
#define rang_mu         (mu_max - mu_min)
#define rang_sigma      (sigma_max - sigma_min)
#define rang_delta      (delta_max - delta_min)

/********** POTÈNCIES PRECALCULADES ********************************************/

#define POW2_64 (18446744073709551616.0) 

/********** ESTRUCTURA DELS INDIVIDUS *****************************************

* Estructura de 7 doubles:
	* 6 Parametres 
	* Valor fitness

*******************************************************************************/
typedef struct Individu{
    double x_0;
    double phi;
    double lambda;
    double mu;
    double sigma;
    double delta;
    double fitness;
} individu;

/********** GENERACIÓ DE NOMBRES ALEATÒRIS *************************************
  
 * randomize  : Selecciona una llavor per als generadors de nombres 
 				pseudo-aleatòris

 * xorshift64 : Generador de unsigned long int pseudo-aleatòri a partir de 
  				l'algorisme xorshift64

 * uniform	  : Dona un double pseudo-aleàtori uniformement distribuit en [a,b]

 *******************************************************************************/

unsigned long int state;
void randomize(void) { state = time(NULL); }
 
unsigned long int xorshift64(void) { state ^= state >> 12;
     state ^= state << 25;
     state ^= state >> 27;
     return state * 0x2545F4914F6CDD1DUL;
}
 
double uniform(double a, double b){
     return a + ((b-a)*xorshift64())/ (POW2_64-1.0);
}

/********** FUNCIONS SA ********************************************************
 
 *  eval_fitness		  : Retorna el valor del funció objectiu avaluada en  
                            un estat a partir de la funció Fit() a "FMD.h"

 *  inicialitzar_estat    : Genera un estat a partir d'una distribució uniforme 
                            al compacte K 

 *  param_bound 		  : Asigna els valors extrems de cada parametre als 
 						    valors que es trobin fora de K 

 *  mecanisme_generacio	  : Mecanisme de generació del SA a partir d'una 
                            distribució uniforme dins d'un rang en funció de t

 *  Inner_loop			  : Du a terme una transició, és a dir genera un estat
                            a partir del mecanisme de generació i selecciona 
                            el següent estat

********************************************************************************/

void eval_fitness(individu *P){
 	
    ODE_Parameters pars;
    pars.psi=P->phi;
    pars.lambda=P->lambda;
    pars.mu= P->mu;
    pars.sigma= P->sigma;
    pars.delta= P->delta;
    pars.beta=ODE_Global_Estructural_Base_Parameters.beta;

    P->fitness=Fit(P->x_0, 25, 12, &pars);
}

void inicialitzar_estat(individu * x){

    (x)->x_0    = uniform(x_0_min,x_0_max) ;
    (x)->phi    = uniform(phi_min,phi_max) ;
    (x)->lambda = uniform(lambda_min,lambda_max) ;
    (x)->mu     = uniform(mu_min,mu_max) ;
    (x)->sigma  = uniform(sigma_min,sigma_max) ;
    (x)->delta  = uniform(delta_min,delta_max) ;

    eval_fitness(x);

}

double param_bound(double param, double a, double b){	

	if(param < a) return a;
	if(param > b) return b;
	return param;

}

void mecanisme_generacio(individu ** xp,  double t) {

    double p=0.001;
   
    if(t>1e3){

        (*xp)->x_0      = uniform(x_0_min,x_0_max) ;
		(*xp)->phi      = uniform(phi_min,phi_max) ;
		(*xp)->lambda   = uniform(lambda_min,lambda_max) ;
		(*xp)->mu       = uniform(mu_min,mu_max) ;
		(*xp)->sigma    = uniform(sigma_min,sigma_max) ;
		(*xp)->delta    = uniform(delta_min,delta_max) ;
    }

    else{

    if(t<1e-4 )    t = t*1e2;

    (*xp)->x_0      += uniform(-p*t*(rang_x_0), p*t*(rang_x_0)); 
    (*xp)->phi      += uniform(-p*t*(rang_phi), p*t*(rang_phi));
    (*xp)->lambda   += uniform(-p*t*(rang_lambda), p*t*(rang_lambda));
    (*xp)->mu       += uniform(-p*t*(rang_mu), p*t*(rang_mu));
    (*xp)->sigma    += uniform(-p*t*(rang_sigma), p*t*(rang_sigma));
    (*xp)->delta    += uniform(-p*t*(rang_delta), p*t*(rang_delta));

    }

    (*xp)->x_0    = param_bound((*xp)->x_0, x_0_min, x_0_max);
    (*xp)->phi    = param_bound((*xp)->phi, phi_min, phi_max);
    (*xp)->lambda = param_bound((*xp)->lambda, lambda_min, lambda_max);
    (*xp)->mu     = param_bound((*xp)->mu, mu_min, mu_max);
    (*xp)->sigma  = param_bound((*xp)->sigma, sigma_min, sigma_max);
    (*xp)->delta  = param_bound((*xp)->delta, delta_min, delta_max);
    
    eval_fitness(*xp);
    
}

float InnerLoop(individu **x, double t){

    unsigned int L = 0, La = 0;

    individu * xp=malloc(sizeof(individu));

    if (xp == NULL) {
        fprintf(stderr, "Error: no s'ha pogut assignar la memòria\n");
        exit(EXIT_FAILURE);
    }

    do { 

        L++;  
        *xp=**x;
        mecanisme_generacio(&xp, t);  

        if(xp->fitness <= (*x)->fitness || uniform(0.0, 1.0) < exp(((*x)->fitness - xp->fitness) / (KB * t))) {
            **x = *xp;  
            La++; 
        }  
        
    } while (L <= Lmax && La <= Lamax);  

    free(xp); 
    return ((float) La) / ((float) L);
    
}

/********** VISUALITZACIÓ RESULTATS *******************************************

 * mostrar_estat 	        : Mostra per pantalla els valors dels paràmetres i fitness 
					         d'un iestat

 * escriure_resultat	   : Escriu en un fitxer els valors dels paràmetres i fitness
					         d'un estat

 * escriure_resultat_taula : Escriu en un fitxer els valors dels paràmetres i fitness
					         en forma de taula per facilitar l'extracció de dades

*******************************************************************************/
void mostrar_estat(individu *P){

	fprintf (stdout, "\nx(0) = %.16lf\n", P->x_0);
    fprintf (stdout, "φ = %.16lf\n", P->phi);
    fprintf (stdout, "λ = %.16lf\n", P->lambda);
    fprintf (stdout, "μ = %.16lf\n", P->mu);
    fprintf (stdout, "σ = %.16lf\n", P->sigma);
    fprintf (stdout, "δ = %.16lf\n\n", P->delta);
    fprintf (stdout, "fitness: %.16lf\n", P->fitness);
    
}
void escriure_resultat(FILE * f, individu Best){

	fprintf(f,"x(0): %.16lf	 ", Best.x_0);
	fprintf(f,"φ: %.16lf   ", Best.phi);
	fprintf(f,"λ: %.16lf   ",Best.lambda);
	fprintf(f,"μ: %.16lf   ", Best.mu);
	fprintf(f,"σ: %.16lf   ",  Best.sigma);
	fprintf(f,"δ: %.16lf   ", Best.delta);
	fprintf(f,"fitness: %.16lf", Best.fitness);
	fprintf(f, "\n\n");

}
void escriure_resultat_taula(FILE * f, individu Best){

	fprintf(f,"%.16lf ", Best.x_0);
	fprintf(f,"%.16lf ", Best.phi);
	fprintf(f,"%.16lf ",Best.lambda);
	fprintf(f,"%.16lf ", Best.mu);
	fprintf(f,"%.16lf ",  Best.sigma);
	fprintf(f,"%.16lf ", Best.delta);
	fprintf(f,"%.16lf ", Best.fitness);

}

/********** FUNCIONS ANÀLISI DE RESULTATS *************************************

 * min_max_scalling : Escalat min-max per igualar les escales dels paràmetres

 * calcul_norma		: Calcula la norma L2 entre els paràmetres de dos individus 

*******************************************************************************/

double min_max_scalling(double diff, double min, double max){
	return (diff-min)/(max-min);
}

double calcul_norma(individu ref, individu p){
	double dx_0 = min_max_scalling(p.x_0,x_0_min,x_0_max)-min_max_scalling(ref.x_0,x_0_min,x_0_max);
	double dphi = min_max_scalling(p.phi,phi_min,phi_max)-min_max_scalling(ref.phi,phi_min,phi_max);
	double dlambda = min_max_scalling(p.lambda,lambda_min,lambda_max)-min_max_scalling(ref.lambda,lambda_min,lambda_max);
	double dmu = min_max_scalling(p.mu,mu_min,mu_max)-min_max_scalling(ref.mu,mu_min,mu_max);
	double dsigma = min_max_scalling(p.sigma,sigma_min,sigma_max)-min_max_scalling(ref.sigma,sigma_min,sigma_max);
	double ddelta = min_max_scalling(p.delta,delta_min,delta_max)-min_max_scalling(ref.delta,delta_min,delta_max);
	double dist = sqrt(dx_0*dx_0 + dphi*dphi + dlambda*dlambda + dmu*dmu + dsigma*dsigma + ddelta*ddelta);
    return dist;
}

/********** SIMULATED ANNEALING (PSEUDOCODI ALGORISME 5.1)**************************/

int main(){
    
    randomize();

    clock_t inici,final;

	inici = clock();
    individu *x=malloc(sizeof(individu));

    if (x == NULL) {
        fprintf(stderr, "Error: no s'ha pogut assignar la memòria\n");
        exit(EXIT_FAILURE);
    }

    inicialitzar_estat(x);
    
    double t = tini;
    float r;
    
    unsigned iter=0U;
    while((r = InnerLoop(&x,t)) < HTsw){ 
        printf("Temperatura: %.3f, Ratio de aceptación (La/L): %.3f\n", t, r);
        t = 2.0*t;
    }
    printf("Maximum temperature found: %lf (La/L = %f)\n", t, r); 
   
       do {
       
       printf("Iteration %4u: temperature: %.8e (La/L = %9.6f)", iter, t, r); 
       t = t*alp;
        
        iter++;
        mostrar_estat(x);
        } while(t > zerot && (r = InnerLoop(&x,t)) >= LTsw);
        printf("Iteration %4u: temperature: %.8e (La/L = %9.6f)", iter, t, r);
        mostrar_estat(x);

        FILE *f = fopen("Resultats_SA_finals.txt", "a");

        escriure_resultat(f,*x);
	    fclose(f);
        
        final=clock();

        //Codi per el fitxer d'anàlisi de resultats:
		//Evaluar la distancia dels parametres a l'optim 
		//Adquirir el temps de computació
        /*
        double temps = ((double) (final - inici)) / CLOCKS_PER_SEC;
	    individu ref = {15670.5560275192783593, 0.2497248909716255, 1570.2313809039706030, 0.0, 0.4904756364357690, 8944.2282749675759987, 2566.999667640135158};
	    //Comparació amb GA vs SA
        char nom_fitxer[64];
        double norma = calcul_norma(ref,*x);
	    if(norma<1e-4){
		    sprintf(nom_fitxer,"Comparació_SA_globals.txt");
	    }
	    else{
		    sprintf(nom_fitxer,"Comparació_SA_locals.txt");
	    }
	
	    f=fopen(nom_fitxer, "a");
	    escriure_resultat(f, *x);
        fprintf(f,"%lf ", norma);
        fprintf(f,"%lf ", temps);
        fprintf(f, "\n\n");
	    fclose(f);*/
        
        free(x);
        return 0;
}
