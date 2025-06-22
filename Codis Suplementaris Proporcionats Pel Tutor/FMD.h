#include <stdio.h>
#include <math.h>
#include <float.h> // For MAXDOUBLE
#include "RKF78.h"

#define MIN0SQ(a) ((a) <  -1 ? (a)*(a) : (a < 0 ? sqrt(-(a)) : 0))
#define ZeRoParsThreshold 1.e-14

typedef struct { double psi; double beta; double lambda; double mu; double sigma; double delta; } ODE_Parameters;
const struct ODE_Base_Parameters { double alpha; double beta; float epsilon; float K;} ODE_Global_Estructural_Base_Parameters = { 0.3489494104672237, 0.0000243826356697, 0.11, 18822.8 };
typedef struct { unsigned short first_data_index; unsigned short number_of_years; } Fitt_period_parameters;

#define NombreDeDades 37
struct data {
    size_t n;
    unsigned short starting_year;
    unsigned char number_of_years_before_predation;
    unsigned char number_of_years_FirstEpoch; // Definition: First epoch is until external big change of population
    unsigned short population[NombreDeDades];
    double prediction[NombreDeDades];
} Birds_Data = {NombreDeDades, 1981, 17, 24, {
/*   1981   1982   1983   1984   1985   1986   1987   1988   1989
 *    [0]    [1]    [2]    [3]    [4]    [5]    [6]    [7]    [8] */
       36,   200,   546,  1200,  1200,  2200,  1850,  2861,  4266,
/*   1990   1991   1992   1993   1994   1995   1996   1997   1998 Migració per depredadors comença el 1998 fins final
 *    [9]   [10]   [11]   [12]   [13]   [14]   [15]   [16]   [17] */
     4300,  3950,  6714,  9373, 10143, 10327, 11328, 11725, 11691,
/*   1999   2000   2001   2002   2003   2004   2005   2006   2007 Increment poblacional per causes externes el 2005
 *   [18]   [19]   [20]   [21]   [22]   [23]   [24]   [25]   [26] Fins 2004 (inclusiu) s'anomena First Epoch; a l'inici (fins 1997) no té migració */
    10189, 10537, 11666, 10122, 10355,  9168, 13988, 15329, 14177,
/*   2008   2009   2010   2011   2012   2013   2014   2015   2016  2017 *** 2018: 1225
 *   [27]   [28]   [29]   [30]   [31]   [32]   [33]   [34]   [35]  [36] *** 2019: 2148 */
    13031,  9762, 11271,  8688,  7571,  6983,  4778,  2067,  1586,  793 } };

double PredictionErrorEuclideanNorm(unsigned short first_year, unsigned short number_of_years){ number_of_years += first_year;
       double d = Birds_Data.population[first_year] - Birds_Data.prediction[first_year], EN = d * d;
       for (register unsigned short i=first_year+1; i < number_of_years; i++) {
            d = Birds_Data.population[i] - Birds_Data.prediction[i]; EN += d * d;
       }
       return sqrt(EN);
}

double RSquaredOfPrediction(unsigned short first_year, unsigned short number_of_years){
       double SumY = Birds_Data.population[first_year], SumY2 = SumY * SumY;
       double EN = (SumY - Birds_Data.prediction[first_year]) * (SumY - Birds_Data.prediction[first_year]);

       for (register unsigned short i=first_year+1; i < first_year+number_of_years; i++) {
            SumY2 += Birds_Data.population[i] * Birds_Data.population[i];
            SumY += Birds_Data.population[i];
            double d = Birds_Data.population[i] - Birds_Data.prediction[i]; EN += d * d;
       }
       return (1.0 - (number_of_years * EN) / (number_of_years * SumY2 - SumY * SumY));
}

void ODE_Prediction_NoMigration(double x0, double alpha, double beta, unsigned short first_year, unsigned short number_of_years){ register unsigned ty;
     double C = x0 * alpha, D = x0 * beta, exppsit = 1.0, exppsi = ((alpha > 0) ? exp(-alpha) : exp(alpha));

     Birds_Data.prediction[first_year] = x0;  // Storing IC x(0)
     if(fabs(alpha) < ZeRoParsThreshold) { for(ty=1; ty < number_of_years; ty++) Birds_Data.prediction[first_year + ty] = x0 / (D * ty + 1); return; }
     if(alpha > 0){ for(ty=1; ty < number_of_years; ty++) { exppsit = exppsit * exppsi; Birds_Data.prediction[first_year + ty] = C / (D + (alpha - D)* exppsit); } return; }
     for(ty=1; ty < number_of_years; ty++) { exppsit = exppsit * exppsi; Birds_Data.prediction[first_year + ty] = (C  * exppsit) / (alpha +  D * (exppsit  - 1)); }
}

void EDO_Prediction_DensoIndependent(double x0, double alpha, double beta, double nu, unsigned short first_year, unsigned short number_of_years){
/* Case nu = 0 */
if(nu < ZeRoParsThreshold) { ODE_Prediction_NoMigration(x0, alpha, beta, first_year, number_of_years); return; }

/* From now on nu > 0 */
    Birds_Data.prediction[first_year] = x0;  // Storing IC x(0)
    double X0beta = x0 * beta;
    double xi = alpha * alpha - 4 * beta * nu;
    register unsigned ty;

/* Ricatti-third case: xi = 0 ; OBSERVE that, in this case, alpha != 0 because nu > 0 (and beta > 0) */
    if(fabs(xi) < ZeRoParsThreshold) {
       double C = 2 * X0beta - alpha, D = (beta * C) / 2; C = (C * alpha) / 4;
       for(ty=1; ty < number_of_years; ty++) Birds_Data.prediction[first_year + ty] = (C * ty + X0beta) / (D * ty + beta);
       return;
    }

/* Ricatti first case: xi > 0 ; AGAIN OBSERVE that, in this case, alpha != 0 because nu > 0 (and beta > 0) */
    if(xi > 0.0){ xi = sqrt(xi);
       double expxit = 1.0, expxi = exp(-xi), A = 2 * nu, B = 2 * X0beta, C = xi + alpha, D = xi - alpha, T = x0 * C - A;
              A = A + x0 * D; D = D + B; C = C - B;
       for(ty=1; ty < number_of_years; ty++) { expxit = expxit * expxi; Birds_Data.prediction[first_year + ty] = (T + A * expxit) / (D + C * expxit); }
       return;
    }

/* And Finally Ricatti second case: xi < 0 */
    xi = sqrt(-xi); // In fact, xi is 2* xi
    double xiX0 = x0 * xi, C = x0 * alpha - 2 * nu, D = 2 * X0beta - alpha, tanxi = tan(xi/2), tanxit = 0.0;
    for(ty=1; ty < number_of_years; ty++) { tanxit = (tanxit + tanxi) / (1 - tanxit * tanxi); Birds_Data.prediction[first_year + ty] = (C * tanxit + xiX0) / (D * tanxit + xi); }
    return;
}

double Fit_DensoIndependent(double x0, double alpha, double beta, double nu, unsigned short first_year, unsigned short number_of_years){
       EDO_Prediction_DensoIndependent(x0, alpha, beta, nu, first_year, number_of_years);
       return PredictionErrorEuclideanNorm(first_year, number_of_years);
}

#define ElliotSigmoidSCALE 1000
#define TwoElliotSigmoidSCALE 2000
double ElliotSigmoid(double x, double sigma, double delta) { x = sigma*(x-delta); return x/(ElliotSigmoidSCALE + fabs(x)); }
double Psi(double x, double mu, double sigma, double delta){ if(fabs(sigma) < ZeRoParsThreshold) return 1.0;
       double ES = ElliotSigmoid(x, sigma, delta);
       sigma *= delta;
       if(x < delta) { x = x / delta; ES = ES * (x + (mu*(1.0-x)*(sigma + ElliotSigmoidSCALE)) / (sigma + TwoElliotSigmoidSCALE)); }
       return ((1 - ES)*(sigma + TwoElliotSigmoidSCALE)) / (sigma*(1+mu) + TwoElliotSigmoidSCALE);
}

/* It is assumed that par->lambda > 0 */
void MigrationODE(double t, double x, double *der, void *Params){ ODE_Parameters *par = (ODE_Parameters *) Params;
     *der = par->psi * x - par->beta*x*x - par->lambda*Psi(x, par->mu, par->sigma, par->delta);
}

#define HMAX  1.0
#define HMIN  1.e-6
#define RKTOL 1.e-8
int ODE_Prediction_WithMigration(double x0, unsigned short first_year, unsigned short number_of_years, ODE_Parameters *pars){ double t = 0.0, err, h = 1.e-3;
    Birds_Data.prediction[first_year] = x0; for (register unsigned short ty=1; ty < number_of_years; ty++) Birds_Data.prediction[first_year + ty] = -1000.0;
    for (register unsigned short ty=1; ty < number_of_years; ty++) {
         while (t+h < ty) {
               int status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
               if(status) { return status; } if(x0 < 0.0) { return 99; }
         }
         do { h = MIN(ty-t, h);
               int status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, pars, MigrationODE);
               if(status) { return status; } if(x0 < 0.0) { return 99; }
         } while(t < ty) ;
         Birds_Data.prediction[first_year + ty] = x0;
    }
    return 0;
}

double Fit(double x0, unsigned short first_year, unsigned short number_of_years, ODE_Parameters *pars){
       if(fabs(pars->lambda) > ZeRoParsThreshold){ if(ODE_Prediction_WithMigration(x0, first_year, number_of_years, pars) == 66) return DBL_MAX; }
       else ODE_Prediction_NoMigration(x0, pars->psi, pars->beta, first_year, number_of_years);
       return PredictionErrorEuclideanNorm(first_year, number_of_years);
}
