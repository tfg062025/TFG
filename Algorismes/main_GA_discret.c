#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "FMD.h"
#include <float.h>
#include <stdbool.h>


/********** HIPERPARAMENTRES GA ***********************************************

* T   	  : Tamany selecció per torneig 
* P_M     : Probabilitat de mutació 
* POPSIZE : Tamany de la població (ha de ser parell)

*******************************************************************************/

#define T 75
#define P_M 0.05
#define POPSIZE 15000

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

/********** LONGITUD DELS GENS ********************************************/

#define l_x_0 19
#define l_phi 32
#define l_lambda 19 
#define l_mu 24
#define l_sigma 19
#define l_delta 16

/********** PAS DE DISCRETITZACIÓ ***************************************/

#define h_x_0 (x_0_max - x_0_min) / (double) ((1UL << l_x_0) - 1UL)
#define h_phi (phi_max - phi_min) / (double)((1UL << l_phi) - 1UL)
#define h_lambda (lambda_max - lambda_min) /(double) ((1UL << l_lambda) - 1UL)
#define h_mu (mu_max - mu_min) / (double)((1UL << l_mu) - 1UL)
#define h_sigma (sigma_max - sigma_min) / (double)((1UL << l_sigma) - 1UL)
#define h_delta (delta_max - delta_min) /(double) ((1UL << l_delta) - 1UL)


/********** GENERACIÓ DE NOMBRES ALEATÒRIS *************************************
  
 * randomize  	 : Selecciona una llavor per als generadors de nombres 
 				   pseudo-aleatòris
 * xorshift64 	 : Generador de unsigned long int pseudo-aleatòri a partir de 
  				   l'algorisme xorshift64
 * u_limited_ran : Generador de unsigned int entre [0,lim-1]
 * bernoulli  	 : Retorna un booleà amb prob true p i prob false 1-p

 *******************************************************************************/

 unsigned long int state;
 void randomize(void) { state = time(NULL); }
 
 unsigned long int xorshift64(void) { state ^= state >> 12;
	 state ^= state << 25;
	 state ^= state >> 27;
	 return state * 0x2545F4914F6CDD1DUL;
 }
 
 unsigned int u_limited_ran (const unsigned long int lim) { 
	unsigned int r;
	r = (xorshift64() % lim);
	return r; 
 }

 bool bernoulli (const double prob) { 
	 return (xorshift64() < ldexp(prob, 64)); 
 }

 
/********** ESTRUCTURA DELS INDIVIDUS *****************************************

* Estructura de 7 components:
	* 6 gens (unsigned int)
	* Valor fitness (double)

*******************************************************************************/

typedef struct Individu{

    unsigned int x_0;
    unsigned int phi;
    unsigned int lambda;
    unsigned int mu;
    unsigned int sigma;
    unsigned int delta;
    double fitness;

} individu;


/********** FUNCIONS GA ********************************************************
 * fenotip				 : A partir de la funció morfogènica obté els valors 
						   dels parametres que formen el fenotip d'un individu

 * inicialitzar_poblacio : Genera una població P(0) a partir de la generació 
 						   de POPSIZE nombres pseudo-aleatòris distribuits 
						   uniformement per a cada parametre dins del compacte K

 * eval_fitness			 : Retorna el valor del fitness d'un individu a partir 
 						   de la funció Fit() a "FMD.h"

 * seleccio_torneig		 : Selecciona POPSIZE individus a partir de la selecció 
						   per torneig amb tamany T

 * param_crossover		 : Aplica el crossover d'un punt a dos gens pares per 
						   obtenir dos gens fills 

 * crossover			 : A partir de param_crossover genera dos cromosomes 
 						   fills

 * param_mut			 : Aplica una perturbació uniforme a cada 
 						   parametre amb una probabilitat P_M i amb 
						   param_bound per mantenir els valors dins de K

********************************************************************************/

double fenotip(unsigned int param, double h , double min){
	return min + param * h;
}

void inicialitzar_poblacio(individu * P){

    for(int i=0;i<POPSIZE;i++){

        P[i].x_0	=	u_limited_ran(1UL << l_x_0);
        P[i].phi	=	u_limited_ran(1UL << l_phi);
        P[i].lambda	=	u_limited_ran(1UL << l_lambda);
        P[i].mu		=	u_limited_ran(1UL << l_mu);
        P[i].sigma	=	u_limited_ran(1UL << l_sigma);
        P[i].delta	=	u_limited_ran(1UL << l_delta);

    }
}

double eval_fitness(individu P){
 	
    ODE_Parameters  pars;

    double x0		=	fenotip(P.x_0, h_x_0, x_0_min);
    pars.psi		=	fenotip(P.phi, h_phi, phi_min);
    pars.lambda		=	fenotip(P.lambda, h_lambda, lambda_min);
    pars.mu			= 	fenotip(P.mu, h_mu, mu_min);
    pars.sigma		= 	fenotip(P.sigma, h_sigma, sigma_min);
    pars.delta		= 	fenotip(P.delta, h_delta, delta_min);
    pars.beta		=	ODE_Global_Estructural_Base_Parameters.beta;
	double fitness	=	Fit(x0, 25, 12, &pars);
    return fitness;

}

void seleccio_torneig(individu *P,individu *best){

	int idx_best=u_limited_ran(POPSIZE);
	int idx_p;

	*best=P[idx_best];

	for(int i=1;i<T;i++){
		idx_p=u_limited_ran(POPSIZE);
		if((best)->fitness>(P[idx_p]).fitness){
			*best=P[idx_p];
		}
	}	
}

unsigned int param_crossover(unsigned int param_a, unsigned int param_b, int l){

	
	unsigned int c = 1 + u_limited_ran(l-1);
	return (param_a >> (l - c)) << (l - c) | (param_b & ~( (~0UL) << (l - c) ));

}

void crossover(individu * c_a, individu *c_b){

	individu  aux;
	aux=*c_a;

	c_a->x_0 = param_crossover(c_a->x_0, c_b->x_0, l_x_0);
	c_b->x_0 = param_crossover(c_b->x_0, aux.x_0, l_x_0);

	c_a->phi = param_crossover(c_a->phi, c_b->phi, l_phi);
	c_b->phi = param_crossover(c_b->phi, aux.phi, l_phi);

	c_a->lambda = param_crossover(c_a->lambda, c_b->lambda, l_lambda);
	c_b->lambda = param_crossover(c_b->lambda, aux.lambda, l_lambda);

	c_a->mu = param_crossover(c_a->mu, c_b->mu, l_mu);
	c_b->mu = param_crossover(c_b->mu, aux.mu, l_mu);

	c_a->sigma = param_crossover(c_a->sigma, c_b->sigma, l_sigma);
	c_b->sigma = param_crossover(c_b->sigma, aux.sigma, l_sigma);

	c_a->delta =param_crossover(c_a->delta, c_b->delta, l_delta);
	c_b->delta =param_crossover(c_b->delta, aux.delta, l_delta);

}

unsigned int param_mut(unsigned int param, int l){

	for(int i=0; i<l; i++){
		if(bernoulli(P_M)) param ^= (1UL<<i);
	}
	return param;
}

void mutacio(individu *p){
	
	p->x_0 = param_mut(p->x_0, l_x_0);
	p->phi = param_mut(p->phi, l_phi);
	p->lambda = param_mut(p->lambda, l_lambda);
	p->mu = param_mut(p->mu, l_mu);
	p->sigma = param_mut(p->sigma, l_sigma);
	p->delta = param_mut(p->delta, l_delta);

}


/********** VISUALITZACIÓ RESULTATS *******************************************

 * mostrar_individu : Mostra per pantalla els valors dels paràmetres i fitness 
					  d'un individu 

 * escriure_resultat: Escriu en un fitxer els valors dels paràmetres i fitness
					  d'un individu

*******************************************************************************/

void mostrar_individu(individu P){

	fprintf (stdout, "x(0) = %.16lf\n", fenotip(P.x_0, h_x_0, x_0_min));
    fprintf (stdout, "φ = %.16lf\n", 	fenotip(P.phi, h_phi, phi_min));
    fprintf (stdout, "λ = %.16lf\n", 	fenotip(P.lambda, h_lambda, lambda_min));
    fprintf (stdout, "μ = %.16lf\n", 	fenotip(P.mu, h_mu, mu_min));
    fprintf (stdout, "σ = %.16lf\n", 	fenotip(P.sigma, h_sigma, sigma_min));
    fprintf (stdout, "δ = %.16lf\n\n", 	fenotip(P.delta, h_delta, delta_min));
	
}

void escriure_resultat(FILE * f, individu P){

	fprintf(f,"x(0): %.16lf	 ", fenotip(P.x_0, h_x_0, x_0_min));
	fprintf(f,"φ: %.16lf   ", fenotip(P.phi, h_phi, phi_min));
	fprintf(f,"λ: %.16lf   ",fenotip(P.lambda, h_lambda, lambda_min));
	fprintf(f,"μ: %.16lf   ", fenotip(P.mu, h_mu, mu_min));
	fprintf(f,"σ: %.16lf   ",  fenotip(P.sigma, h_sigma, sigma_min));
	fprintf(f,"δ: %.16lf   ", fenotip(P.delta, h_delta, delta_min));
	fprintf(f,"fitness: %.16lf", P.fitness);
	fprintf(f, "\n\n");

}

/******** ALGORISME GENÈTIC (PSEUDOCODI ALGORISME 4.1)************************/

void GA(individu *Best){
	
	int contador=0;
	int generacio=0;

	individu * P=malloc((2*POPSIZE) *sizeof(individu ));
	individu * Q=(P+POPSIZE);
	individu *orig = P;
	inicialitzar_poblacio(P);

	double bestfitness=DBL_MAX;
	double bestfitness_ant=DBL_MAX;

	for(int j=0;j<POPSIZE;j++){
	
		P[j].fitness=eval_fitness(P[j]);
		if(P[j].fitness<bestfitness){

			bestfitness=P[j].fitness;
			*Best=P[j];
			
		}		
	}
	
	while(contador<100){
		
		bestfitness_ant=bestfitness;

		for(int j=0;j<POPSIZE/2;j++){

			
			seleccio_torneig(P,Q+2*j);
			seleccio_torneig(P,Q+(2*j+1));

			crossover(Q+2*j,Q+(2*j+1));

			mutacio(Q+2*j);
			mutacio(Q+(2*j+1));

			Q[2*j].fitness=eval_fitness(Q[2*j]);
			Q[2*j+1].fitness=eval_fitness(Q[2*j+1]);
			
			if(Q[2*j].fitness<bestfitness){
				bestfitness=Q[2*j].fitness;
				*Best=Q[2*j];
			}	
			if(Q[2*j+1].fitness<bestfitness){
				bestfitness=Q[2*j+1].fitness;
				*Best=Q[2*j+1];
			}
	
		}
		
		individu *aux;
		aux=P;
		P=Q;
		Q=aux;

		generacio++;

		if(bestfitness==bestfitness_ant){
			contador++;
		}
		else{
			contador=0;
		}
		printf("\nFitness del millor individu obtingut a la generació %i:%f \n", generacio,(Best)->fitness);
		printf("Nombre de generacions sense millorar el fitness %i \n",contador);
		
	}

	free(orig);

}

/********** CÀLCUL DE L'ÒPTIM ***********************************************/

int main(){
	
	if(POPSIZE % 2 == 1){
		printf("ERROR: La mida de la població ha de ser un nombre parell \n" );
		exit(EXIT_FAILURE);
	}
	randomize();
	individu Best;
	GA(&Best);
	printf("\nSolució obtinguda:\n");
	mostrar_individu(Best);	

	char nom_fitxer[64];

	sprintf(nom_fitxer,"GA_discret_P%i_PM%.3lf_T%i.txt", POPSIZE, P_M,T);
	FILE * f = fopen(nom_fitxer, "a");

	if (f == NULL) {
		printf("No s'ha pogut obrir el fitxer %s\n", nom_fitxer);
		
	}

	escriure_resultat(f, Best);
	fclose(f);
	return 0;

}