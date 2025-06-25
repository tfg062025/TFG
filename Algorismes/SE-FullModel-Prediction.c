#include <stdio.h>
#include <string.h>
#include "FMD.h"
char signe(double r, double s){ if(fabs(r-s) < 1.e-12) return '='; return ((s < r) ? '-' : '+'); }

int main (int argc, char *argv[]){ Fitt_period_parameters period = { 25, 12 }; unsigned short print_full_solution = 0;

    if (argc < 7 || argc > 8) { fprintf (stderr, "\nUSAGE: %s [-FS|--full-solution] x0 φ λ μ σ δ\n\n", argv[0]); return 1; }
    if(strcmp(argv[1], "-FS") == 0 || strcmp(argv[1], "--full-solution") == 0) print_full_solution = 1;

    double x0 = atof(argv[1+print_full_solution]);
    double psi = ODE_Global_Estructural_Base_Parameters.alpha;
    if(strcmp(argv[2+print_full_solution], "alpha") && strcmp(argv[2+print_full_solution], "α")) psi = atof(argv[2+print_full_solution]);
    double rho = ODE_Global_Estructural_Base_Parameters.alpha - psi;
    ODE_Parameters pars = { psi, ODE_Global_Estructural_Base_Parameters.beta, atof(argv[3+print_full_solution]), atof(argv[4+print_full_solution]), atof(argv[5+print_full_solution]), atof(argv[6+print_full_solution]) };
    if( x0 < 1 || psi < 0.1 || rho < 0 || pars.lambda < 0 || pars.mu < 0 || pars.sigma < 0 || pars.delta < 0) { fprintf (stderr, "\nUSAGE: %s [-FS|--full-solution] x0 ρ λ μ σ δ\n\n", argv[0]); return 1; }

    fprintf(stdout,"\n#### Evaluation of the Second Epoch (period: %u-%u) Fitness Function\n",
                           Birds_Data.starting_year + period.first_data_index, Birds_Data.starting_year + period.first_data_index + period.number_of_years - 1 );
    fprintf (stdout, "#### (Database interval data: first year = %u ; last year = %u and number of years = %u)\n",
                           period.first_data_index, period.first_data_index + period.number_of_years - 1, period.number_of_years);
    fprintf (stdout, "#### Initial condition x(0) = %.16lf\n", x0);
    fprintf (stdout, "#### β (fixed) = %e\n", pars.beta);
    fprintf (stdout, "#### φ = %.16lf (ρ = %.16lf)\n", pars.psi, rho);
    fprintf (stdout, "#### λ = %.16lf\n", pars.lambda);
    fprintf (stdout, "#### μ = %.16lf\n", pars.mu);
    fprintf (stdout, "#### σ = %.16lf\n", pars.sigma);
    fprintf (stdout, "#### δ = %.16lf\n\n", pars.delta);

    fprintf (stdout, "Quadratic error: %.16lf\n", Fit(x0, period.first_data_index, period.number_of_years, &pars));
    fprintf (stdout, "R^2 of the prediction = %lf\n", RSquaredOfPrediction(period.first_data_index, period.number_of_years));

    double SNSQ = 0.0, SPSQ = 0.0;
    for (register unsigned short i=period.first_data_index; i < period.first_data_index + period.number_of_years; i++) {
         double diff = Birds_Data.population[i] - Birds_Data.prediction[i];
         if (diff < 0.0) { SNSQ += diff*diff; } else if (diff > 0.0) { SPSQ += diff*diff; }
    }
    fprintf (stdout, "Quadratic-positive error: %.15lf\nQuadratic-negative error: %.15lf\n", sqrt(SNSQ), sqrt(SPSQ));

    fprintf(stdout, "\nPredictions\n-------------------------------------------\n Year Observed        Predicted       Dpos\n===========================================\n");
    for (register unsigned short i=period.first_data_index; i < period.first_data_index + period.number_of_years; i++) fprintf(stdout, "%5u %8d %21.15lf pos: %c\n",
                                 i+Birds_Data.starting_year, Birds_Data.population[i], Birds_Data.prediction[i], signe(Birds_Data.population[i], Birds_Data.prediction[i]) );

    if(print_full_solution) {
       fprintf (stdout, "\n\n---------------------------------------------------------\nFull ODE solution in the Second Epoch (period: %u-%u)\n",
                             Birds_Data.starting_year + period.first_data_index, Birds_Data.starting_year + period.first_data_index + period.number_of_years - 1 );
       fprintf (stdout, "       t                x(t)               Error\n---------------------------------------------------------\n");
       fprintf (stdout, "           0 %25.16lf %.12g\n", x0, 0.0);

       double t = 0.0, err, h = 1.e-3;
       FILE *f = fopen("../Resultats/Solució EDO/Solució EDO.txt", "w");
       if(f==NULL){
          printf("Error: No s'ha pogut obrir el fitxer");
       }
       for (register unsigned short ty=1; ty < period.number_of_years; ty++) {
            while (t+h < ty) {
                   int status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, &pars, MigrationODE);
                   fprintf(f, "%.16lf %.16lf\n",t,x0);
                   if(status) { fprintf (stdout, "RKF78 ERROR: Error number %d; Stopping\n", status); exit(status); }
                   fprintf (stdout, "%12.8lf %25.16lf %.12g\n", t, x0, err);
                   if(x0 < 0.0) exit(99);
            }
            do { h = MIN(ty-t, h);
                 int status = RKF78(&t, &x0, &h, &err, HMIN, HMAX, RKTOL, &pars, MigrationODE);
                 fprintf(f, "%.16lf %.16lf\n",t,x0);
                 if(status) { fprintf (stdout, "RKF78 ERROR: Error number %d; Stopping\n", status); exit(status); }
                 if(t < ty) { fprintf (stdout, "%12.8lf %25.16lf %.12g\n", t, x0, err); if(x0 < 0.0) exit(99); }
            } while(t < ty) ;
            fprintf (stdout, "%12u %25.16lf %.12g\n", ty, x0, err);
            
            if(x0 < 0.0) exit(99);
       }

    }

    fprintf (stdout, "\n");
    
    return 0;
}
