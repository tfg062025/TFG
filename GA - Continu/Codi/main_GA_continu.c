#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "FMD.h"
#include <float.h>
#include <stdbool.h>

/********** HIPERPARAMETRES GA ***********************************************

* T   	  : Tamany selecció per torneig 
* P_M     : Probabilitat de mutació 
* POPSIZE : Tamany de la població (Ha de ser parell)

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

/********** POTÈNCIES PRECALCULADES ********************************************/

#define POW2_64 (18446744073709551616.0) 

/********** GENERACIÓ DE NOMBRES ALEATÒRIS *************************************
  
 * randomize  : Selecciona una llavor per als generadors de nombres 
 				pseudo-aleatòris
 * xorshift64 : Generador de unsigned long int pseudo-aleatòri a partir de 
  					  l'algorisme xorshift64
 * uniform	  : Dona un double pseudo-aleàtori uniformement distribuit en [a,b]
 * bernoulli  : Retorna un booleà amb prob true p i prob false 1-p

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

bool bernoulli (const double prob) { 
	return (xorshift64() < ldexp(prob, 64)); 
}

/********** ESTRUCTURA DELS INDIVIDUS *****************************************

* Estructura de 7 doubles:
	* 6 parametres/gens
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

/********** FUNCIONS GA ********************************************************
 
 * inicialitzar_poblacio : Genera una població P(0) a partir de la generació 
 						   de POPSIZE nombres pseudo-aleatòris distribuits 
						   uniformement per a cada parametre dins del compacte K

 * eval_fitness			 : Retorna el valor del fitness d'un individu a partir 
 						   de la funció Fit() a "FMD.h"

 * seleccio_torneig		 : Selecciona POPSIZE individus a partir de la selecció 
						   per torneig amb tamany T

 * param_crossover		 : Usa la recombinació lineal extesà amb parametre p 
 						   per generar 2 fills a partir de 2 progenitors per un 
						   parametre concret

 * param_bound 			 : Asigna els valors extrems de cada parametre als 
 						   valors que es trobin fora de K 

 * crossover			 : A partir de param_crossover i param_bound genera dos 
 						   fills a K per a dos progenitors donats

 * mutacio				 : Aplica una perturbació uniforme en [-r, r] a cada 
 						   parametre amb una probabilitat P_M i aplica 
						   param_bound per mantenir els valors a K

********************************************************************************/

void inicialitzar_poblacio(individu * P){

	for(int i=0;i<POPSIZE;i++){

		P[i].x_0		= 	uniform(x_0_min,x_0_max);
		P[i].phi 		= 	uniform(phi_min,phi_max);
		P[i].lambda 	=	uniform(lambda_min,lambda_max);
		P[i].mu 		=	uniform(mu_min,mu_max);
		P[i].sigma 		= 	uniform(sigma_min,sigma_max);
		P[i].delta 		=	uniform(delta_min,delta_max);
		
    }
}

double eval_fitness(individu P){
 	
    ODE_Parameters pars;
    pars.psi=P.phi;
    pars.lambda=P.lambda;
    pars.mu= P.mu;
    pars.sigma= P.sigma;
    pars.delta= P.delta;
    pars.beta=ODE_Global_Estructural_Base_Parameters.beta;

	double fitness=Fit(P.x_0, 25, 12, &pars);
    return fitness;
}

void seleccio_torneig(individu *P,individu *best){

	int idx_best = rand() % POPSIZE;
	int idx_p;

	*best=P[idx_best];

	for(int i=1;i<T;i++){

		idx_p =	rand() % POPSIZE;

		if(best->fitness > P[idx_p].fitness){

			*best = P[idx_p];
		}
	}	
}

double param_crossover(double param_a, double param_b, double p){

	return param_a + uniform(-p,p) * (param_b - param_a);
}

double param_bound(double param, double a, double b){	

	if(param < a) return a;
	if(param > b) return b;
	return param;

}

void crossover(individu * c_a, individu *c_b, double p){

	individu aux;
	aux=*c_a;
	
	c_a->x_0	=	param_crossover(c_a->x_0, c_b->x_0, p);
	c_b->x_0	=	param_crossover(c_b->x_0, aux.x_0, p);

	c_a->x_0	=	param_bound(c_a->x_0, x_0_min, x_0_max);
	c_b->x_0	=	param_bound(c_b->x_0, x_0_min, x_0_max);

	c_a->phi	=	param_crossover(c_a->phi, c_b->phi, p);
	c_b->phi	=	param_crossover(c_b->phi, aux.phi, p);

	c_a->phi	=	param_bound(c_a->phi, phi_min, phi_max);
	c_b->phi	=	param_bound(c_b->phi, phi_min, phi_max);

	c_a->lambda	=	param_crossover(c_a->lambda, c_b->lambda, p);
	c_b->lambda	=	param_crossover(c_b->lambda, aux.lambda, p);

	c_a->lambda	=	param_bound(c_a->lambda, lambda_min, lambda_max);
	c_b->lambda	=	param_bound(c_b->lambda, lambda_min, lambda_max);

	c_a->mu		=	param_crossover(c_a->mu, c_b->mu, p);
	c_b->mu		=	param_crossover(c_b->mu, aux.mu, p);

	c_a->mu		=	param_bound(c_a->mu, mu_min, mu_max);
	c_b->mu		=	param_bound(c_b->mu, mu_min, mu_max);

	c_a->sigma	=	param_crossover(c_a->sigma, c_b->sigma, p);
	c_b->sigma	=	param_crossover(c_b->sigma, aux.sigma, p);

	c_a->sigma	=	param_bound(c_a->sigma, sigma_min, sigma_max);
	c_b->sigma	=	param_bound(c_b->sigma, sigma_min, sigma_max);

	c_a->delta	=	param_crossover(c_a->delta, c_b->delta, p);
	c_b->delta	=	param_crossover(c_b->delta, aux.delta, p);

	c_a->delta	=	param_bound(c_a->delta, delta_min, delta_max);
	c_b->delta	=	param_bound(c_b->delta, delta_min, delta_max);

}

void mutacio(individu *p, int generacio, double r) {
	
    if (bernoulli(P_M)) {

		p->x_0 	+= 	uniform(- r*(x_0_max-x_0_min) , r*(x_0_max-x_0_min));
		p->x_0	=	param_bound(p->x_0, x_0_min,x_0_max);

	}

	if (bernoulli(P_M)) {

		p->phi 	+= 	uniform(- r*(phi_max-phi_min) , r*(phi_max-phi_min));
		p->phi	=	param_bound(p->phi, phi_min,phi_max);

	}

	if (bernoulli(P_M)) {

		p->lambda 	+= 	uniform(- r*(lambda_max-lambda_min) , r*(lambda_max-lambda_min));
		p->lambda	=	param_bound(p->lambda, lambda_min,lambda_max);

	}

	if (bernoulli(P_M)) {

		p->mu 	+= 	uniform(- r*(mu_max-mu_min) , r*(mu_max-mu_min));
		p->mu	=	param_bound(p->mu, mu_min,mu_max);

	}

	if (bernoulli(P_M)) {

		p->sigma 	+= 	uniform(- r*(sigma_max-sigma_min) , r*(sigma_max-sigma_min));
		p->sigma	=	param_bound(p->sigma, sigma_min,sigma_max);

	}

	if (bernoulli(P_M)) {

		p->delta 	+= 	uniform(- r*(delta_max-delta_min) , r*(delta_max-delta_min));
		p->delta	=	param_bound(p->delta, delta_min,delta_max);

	}

}

/********** VISUALITZACIÓ RESULTATS *******************************************

 * mostrar_individu 	   : Mostra per pantalla els valors dels paràmetres i fitness 
					         d'un individu 

 * escriure_resultat	   : Escriu en un fitxer els valors dels paràmetres i fitness
					         d'un individu

 * escriure_resultat_taula : Escriu en un fitxer els valors dels paràmetres i fitness
					         en forma de taula per facilitar l'extracció de dades

*******************************************************************************/

void mostrar_individu(individu P){

	fprintf (stdout, "x(0) = %.16lf\n", P.x_0);
    fprintf (stdout, "φ = %.16lf\n", P.phi);
    fprintf (stdout, "λ = %.16lf\n", P.lambda);
    fprintf (stdout, "μ = %.16lf\n", P.mu);
    fprintf (stdout, "σ = %.16lf\n", P.sigma);
    fprintf (stdout, "δ = %.16lf\n\n", P.delta);

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

/******** ALGORISME GENÈTIC (PSEUDOCODI ALGORISME 4.1)************************/

void GA(individu *Best){

	int contador=0;
	int generacio=0;
	
	individu * P=malloc((2*POPSIZE) *sizeof(individu ));
	individu * Q=(P+POPSIZE);

	if (P == NULL) {
        fprintf(stderr, "Error: no s'ha pogut assignar la memòria\n");
        exit(EXIT_FAILURE);
    }

	if (Q == NULL) {
        fprintf(stderr, "Error: no s'ha pogut assignar la memòria\n");
        exit(EXIT_FAILURE);
    }

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
		

		char nom_hist[64];
	
		sprintf(nom_hist,"generació_%i.txt", generacio);
		FILE * hist = fopen(nom_hist, "w");
		for(int i=0;i<POPSIZE;i++){
			escriure_resultat_taula(hist,P[i]);
			fprintf(hist,"\n");
		}
		fclose(hist);
		bestfitness_ant=bestfitness;
		
		for(int j=0;j<POPSIZE/2;j++){
			
			
			seleccio_torneig(P,Q+2*j);
			seleccio_torneig(P,Q+(2*j+1));
			
			double p=(generacio)/20.0;
			if(p>1) p=1;
			
			crossover(Q+2*j,Q+(2*j+1),p);

			mutacio(Q+2*j,generacio, 0.05);
			mutacio(Q+(2*j+1),generacio, 0.05);

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
		
		// Codi per generar l'evolució del fitness

		/*char nom_fitxer_evolucio[64];

		sprintf(nom_fitxer_evolucio,"GA_continu_evolucio_NOU_P%i_PM%.2lf_T%i.txt", POPSIZE, P_M, T);
		FILE * fe = fopen(nom_fitxer_evolucio, "w");

		if (fe == NULL) {
			printf("No s'ha pogut obrir el fitxer %s\n", nom_fitxer_evolucio);
		}

		fprintf(fe, "%lf\n", Best->fitness);

		fclose(fe);*/

		individu *aux;
		aux = P;
		P = Q;
		Q = aux;
		generacio++;

		if(fabs(bestfitness-bestfitness_ant)<1e-15) contador++;
		else contador=0;

		printf("Fitness del millor individu obtingut a la generació %i:%f \n", generacio,(Best)->fitness);
		printf("Nombre de generacions sense millorar el fitness %i \n\n",contador);
	
	}

	free(orig);
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

/********** CÀLCUL DE L'ÒPTIM ***********************************************/

int main(){
	
	if(POPSIZE % 2 == 1){
		printf("ERROR: La mida de la població ha de ser un nombre parell \n" );
		exit(EXIT_FAILURE);
	}

	clock_t inici,final;

	inici = clock();
	
	randomize();

	individu Best;
	GA(&Best);

	printf("Solució obtinguda:\n");
	mostrar_individu(Best);

	char nom_fitxer[64];
	
	sprintf(nom_fitxer,"GA_continu_P%i_PM%.2lf_T%i.txt", POPSIZE, P_M, T);
	FILE * f = fopen(nom_fitxer, "a");

	if (f == NULL) {
		printf("No s'ha pogut obrir el fitxer %s\n", nom_fitxer);
		
	}

	escriure_resultat(f, Best);
	fclose(f);

	final = clock();
	
	//Codi per el fitxer d'anàlisi de resultats:
		//Evaluar la distancia dels parametres a l'optim 
		//Adquirir el temps de computació
	
	/*double temps = ((double) (final - inici)) / CLOCKS_PER_SEC;
	individu ref = {15670.5560275192783593, 0.2497248909716255, 1570.2313809039706030, 0.0, 0.4904756364357690, 8944.2282749675759987, 2566.999667640135158};
	
	double norma = calcul_norma(ref,Best);

	if(calcul_norma(ref,Best)<1e-4){
		sprintf(nom_fitxer,"Comparació_GA_globals.txt");
	}
	else{
		sprintf(nom_fitxer,"Comparació_GA_locals.txt");
	}
	
	f=fopen(nom_fitxer, "a");
	escriure_resultat_taula(f, Best);
	fprintf(f,"%lf ", norma);

	fprintf(f,"%lf ", temps);
    fprintf(f, "\n\n");

	
	fclose(f);
	*/
	
	return 0;
	
}

