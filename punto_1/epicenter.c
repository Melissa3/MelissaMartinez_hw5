#include <stdio.h>
#include <stdlib.h>
#include <math.h>
# define NSUM 25	

double likelihood(double *obs, double *model, int n);
double fRand(double fMin, double fMax);
double gaussrand();

int main()
{
	FILE *out_file = fopen("iteraciones.csv","w"); 

	int n=6;
	int n_iterations= 20000;
	int i=0, j=0;
	double x[n], y[n], t[n];
	double beta;
	float sigma_t=0.1;
	double h_0=25.0, h_new;
	double k_0=40.0, k_new;
	double lhood_0, lhood_new;
	double r[n], r_0[n], r_new[n];
	

	double v= 5.0;

	FILE *fstream = fopen("datos.csv", "r");


	for (i = 0; i < n; i++)
	{
		fscanf(fstream, "%lf,%lf,%lf\n", &x[i], &y[i], &t[i]);
		r[i]= v*t[i];
	}	
	fclose(fstream);

	for (i=0; i<n; i++)
	{
		r_0[i]= sqrt(pow(x[i]-h_0,2) + pow(y[i]-k_0,2));
	}

	lhood_0= likelihood(r, r_0, n);

	for (i=0; i<n_iterations; i++)
	{
	h_new = h_0 + gaussrand();
	k_new = k_0 + gaussrand();

		for (j=0; j<n; j++)
		{
		r_new[j]= sqrt(pow(x[j]-h_new,2) + pow(y[j]-k_new,2));
		}

	lhood_new = likelihood(r, r_new, n);

	double alpha= lhood_new/ lhood_0; 

	if(alpha>=1.0)
	{
	h_0= h_new;
	k_0= k_new;
	lhood_0= lhood_new;
	}
	else
	{
	beta= fRand(0.0, 1.0);
		if (beta<= alpha)
		{
		h_0= h_new;
		k_0= k_new;
		lhood_0= lhood_new;		
		}

	}
fprintf(out_file, "%lf,%lf\n", h_0, k_0);
}
fclose(out_file);	
}

double likelihood(double *obs, double *model, int n)
{
   double c=0.0;
   int i;
   for (i=0; i<n; i++)
   {
	c = c + pow(obs[i]-model[i], 2);
   }
   c= c*(1.0/2.0);
   return exp(-c); 
}

double fRand(double fMin, double fMax)
{
/* Generate random double from uniform distribution
http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c */
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


double gaussrand()
{
/* Generate random double from Gaussian distribution
http://c-faq.com/lib/gaussian.html*/
	double x = 0;
	int i;
	for(i = 0; i < NSUM; i++)
		x += (double)rand() / RAND_MAX;

	x -= NSUM / 2.0;
	x /= sqrt(NSUM / 12.0);

	return x;
}


