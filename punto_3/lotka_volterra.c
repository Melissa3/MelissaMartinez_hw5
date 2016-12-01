#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

void leap_frog (double alpha, double beta, double gamma, double delta, double *x_prey, double *y_predator, double *x_0, double *y_0, int N, double *time, double dtime);

double new_likelihood(double *T, double *x, double *y, int N, double *t, double *x_obs, double *y_obs, int n, double dt);

double fRand(double fMin, double fMax);

void interp_1d(double *x_original, double *f_original, int n, double *x_interp, double *f_interp, int m, double dx);

int main()
{
	int i=0, n= 96, n_iterations=20000;
	double t[n], x_prey[n], y_predator[n];
	double alpha_0= 20.0 , beta_0= 5.0 , gamma_0= 20.0 , delta_0= 0.5;
	double alpha_new, beta_new, gamma_new, delta_new;

	double lhood_0, lhood_new;

	int N=800;
	double time[N];
	time[0]= 0.0;
	double t0=0.0, tf=0.8;
	double dtime=(tf-t0)/N;	
	double x_0[N], x_new[N];
	double y_0[N], y_new[N];

	for (i=1; i<N; i++)
	{
	time[i]= time[i-1] + dtime; 
	}	
	

	FILE *out_file = fopen("parametros.csv","w"); 

	FILE *fdatos = fopen("lotka_volterra_obs.csv", "r");

	for (i = 0; i < n; i++)
	{
		fscanf(fdatos, "%lf,%lf,%lf\n", &t[i], &x_prey[i], &y_predator[i]);
	}	
	fclose(fdatos);

	leap_frog(alpha_0, beta_0, gamma_0, delta_0, x_prey, y_predator, x_0, y_0, N, time, dtime);
	
	lhood_0= new_likelihood(time, x_0, y_0, N, t, x_prey, y_predator, n, 5*dtime);

	/*lhood_0= new_likelihood(x_prey, y_predator, x_0, y_0, n, N, t, time, dtime);*/

	for (i=0; i< n_iterations; i++)
	{

	alpha_0 = alpha_0 + fRand(-1.0, 1.0);

	beta_0 = beta_0 + fRand(-0.5, 0.5);

	gamma_0 = gamma_0 + fRand(-1.0, 1.0);

	delta_0 = delta_0 + fRand(-0.1, 0.1); 

	leap_frog(alpha_new, beta_new, gamma_new, delta_new, x_prey, y_predator, x_new, y_new, N, time, dtime);

	lhood_new= new_likelihood(time, x_new, y_new, N, t, x_prey, y_predator, n, 5*dtime);

	double ratio= lhood_new/ lhood_0; 

		if(ratio>=1.0)
		{
		alpha_0 = alpha_new;
		beta_0 = beta_new;
		gamma_0 = gamma_new;
		delta_0 = delta_new;
		lhood_0= lhood_new;
		}
		else
		{
		double opt= fRand(0.0, 1.0);
			if (opt<= ratio)
			{
			alpha_0 = alpha_new;
			beta_0 = beta_new;
			gamma_0 = gamma_new;
			delta_0 = delta_new;
			lhood_0= lhood_new;		
			}
		}
	fprintf(out_file, "%lf,%lf,%lf,%lf,%lf\n", alpha_0, beta_0, gamma_0, delta_0, lhood_0);
	}
fclose(out_file);
}

void leap_frog (double alpha, double beta, double gamma, double delta, double *x_prey, double *y_predator, double *x0, double *y0, int N, double *time, double dtime)
{
	int i;
	double dx_dt, dy_dt;

	x0[0]= x_prey[0];
	y0[0]= y_predator[0];

	x0[1]= x_prey[1];
	y0[1]= y_predator[1];
	
	for (i=1; i< N-1; i++)
	{
	dx_dt= x0[i]*(alpha - beta*y0[i]);
	dy_dt= -y0[i]*(gamma - delta*x0[i]);

	x0[i+1]= x0[i-1] + dtime*dx_dt;
	y0[i+1]= y0[i-1] + dtime*dy_dt;			
	}
}

double new_likelihood(double *T, double *x, double *y, int N, double *t, double *x_obs, double *y_obs, int n, double dt)
{
	int i = 0;
	double x_interp[n], y_interp[n];
	double l=0.0;

	interp_1d(T, x, N, t, x_interp, n, dt);
	interp_1d(T, y, N, t, y_interp, n, dt);

	for (i = 0; i < n; ++i)
	{
		l = l + pow( pow(x_obs[i]-x_interp[i],2) + pow(y_obs[i]-y_interp[i],2), 0.5);
	}
	return exp(-0.5*l/n);
}


void interp_1d(double *x_original, double *f_original, int n, double *x_interp, double *f_interp, int m, double dx)
{
	int i, i2, j=0;
	double a, b;
	double delta_x;
	double delta_x2;

	for (i = 0; i < m; ++i)
	{
		i2 = i;
		j=0;
		b = 0.0;
		while(j<n)
		{
			delta_x = x_interp[i]-x_original[j];
			delta_x2 = x_original[j+1]-x_interp[i];

			if (delta_x > 0.0 && delta_x<dx)
			{				
				b = f_original[j];
				a = (f_original[j+1]-f_original[j]) / (x_original[j+1]-x_original[j]);
				f_interp[i] = a*delta_x + b;
				break;
			}
			else if(delta_x2>0.0 && delta_x2<dx)
			{
				b = f_original[j];
				a = (f_original[j+1]-f_original[j]) / (x_original[j+1]-x_original[j]);
				f_interp[i] = a*(x_original[j+1]-x_original[j]-delta_x2) + b;
				break;
			}
			j = j+1;			
		}		
	}
}

double fRand(double fMin, double fMax)
{
/* Generate random double from uniform distribution
http://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c */
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


