#include <iostream>
#include <complex>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
using namespace std;

//simulation parameters
int N = 401;
int Ntime = 100000000;
double C = 1;
double a = 1.0;
double m = 1.0;
double k_b = 1.0;
double sigma = 30;
double dt = 0.02*M_PI*sqrt(m/C)*0.1;
double sound_vel = a*sqrt(C/m);
double k0 = 2*M_PI/(N-1)*30;
double k1 = 2*M_PI/(N-1)*30;
double wp_dist = 80;


//open output file
int output_interval = 10000;
ofstream file_data;
ofstream file_data2;
ofstream file_data3;

vector<double> ddu (N);
vector<double> du (N);
vector<double> u (N);

double pair_force(double z, double a, double C) {
    double alpha = 1;
    double beta = 0;
    double force = C*(a-z) + alpha*pow((a-z), 2.0) + beta*pow((a-z), 3.0);
    return force;
}

double initial_displacement (double n, double a) {
	double u_0 = 0;
	double u_1 = 0;
	//u_1 = sin(k1*n);
	u_0 = exp(-pow((n-(N-1)/2+wp_dist), 2)/2/pow(sigma,2))*sin(k0*n); 
	u_1 = exp(-pow((n-(N-1)/2-wp_dist), 2)/2/pow(sigma,2))*sin(k1*n); 
    return u_0+u_1;
}

double initial_velocities (double n, double a) {
	double u_0 = 0;
	double u_1 = 0;
	u_0 = exp(-pow((n-(N-1)/2+wp_dist), 2)/2/pow(sigma,2))*k0*cos(k0*n) - pow(sigma, -2)*exp(-pow((n-(N-1)/2)+wp_dist, 2)/2/pow(sigma,2))*(n-(N-1)/2+wp_dist)*sin(k0*n); 
	u_1 = exp(-pow((n-(N-1)/2)-wp_dist, 2)/2/pow(sigma,2))*k1*cos(k1*n) - pow(sigma, -2)*exp(-pow((n-(N-1)/2) -wp_dist, 2)/2/pow(sigma,2))*(n-(N-1)/2-wp_dist)*sin(k1*n); 
	u_0 *= -sound_vel;
	u_1 *= sound_vel;
    return u_0+u_1;
}

int main(int argc, const char * argv[]) {

    for (int j = 0; j<N; j++) {
		u[j] = initial_displacement(j,a);
		du[j] = initial_velocities(j,a);
	}

    //loop
    for (int t=0; t<Ntime; t++) {
        if(t % output_interval == 0) {
	        string filename;		
			filename = "data/" + to_string(t) + ".dump";
			file_data.open (filename);
            for(int j=0; j<N; j++) file_data << j << " " << fixed << u[j] <<"\n";
			file_data.close();

		    // Calculate DFT of x using brute force
			int k, j ;
			float Xre[N], Xim[N], Xim_dot[N],E[N] ; // DFT of x (real and imaginary parts)
 		    float P[N];           // power spectrum of x

			for (int k=0 ; k<N ; ++k)
			{
				// Real part of X[k]
				Xre[k] = 0;
				for (int j=0 ; j<N ; ++j) Xre[k] += u[j] * cos(j * k * M_PI * 2 / (N));
		
				// Imaginary part of X[k]
				Xim[k] = 0;
				for (int j=0 ; j<N ; ++j) Xim[k] -= u[j] * sin(j * k * M_PI * 2 / (N));
		
				Xim_dot[k] = 0;
				for (int j=0 ; j<N ; ++j) Xim_dot[k] -= du[j] * sin(j * k * M_PI*2/ (N));

				E[k] = 0;
				E[k] = 0.5*pow(Xim_dot[k],2) +2 *  pow(Xim[k],2)*sin(k * M_PI*2/ (2*(N)));

				// Power at kth frequency bin
				P[k] = Xre[k]*Xre[k] + Xim[k]*Xim[k];
			}


			string filename1;		
			filename1 = "data2/" + to_string(t) + ".dump";
			file_data2.open (filename1);
            for(int k=0; k<N/2; k++) file_data2 << k << " " << fixed << P[k] <<"\n";
			file_data2.close();

			string filename2;		
			filename2 = "data3/" + to_string(t) + ".dump";
			file_data3.open (filename2);
            for(int k=0; k<N/2; k++) file_data3 << k << " " << fixed << E[k] <<"\n";
			file_data3.close();

			// Output results to MATLAB / Octave M-file for plotting
			FILE *f = fopen("dftplots.m", "w");
			fprintf(f, "n = [0:%d];\n", N-1);
			fprintf(f, "x = [ ");
			for (j=0 ; j<N ; ++j) fprintf(f, "%f ", u[j]);
			fprintf(f, "];\n");
			fprintf(f, "Xre = [ ");
			for (k=0 ; k<N ; ++k) fprintf(f, "%f ", Xre[k]);
			fprintf(f, "];\n");
			fprintf(f, "Xim = [ ");
			for (k=0 ; k<N ; ++k) fprintf(f, "%f ", Xim[k]);
			fprintf(f, "];\n");
			fprintf(f, "P = [ ");
			for (k=0 ; k<N ; ++k) fprintf(f, "%f ", P[k]);
			fprintf(f, "];\n");
			fprintf(f, "subplot(3,1,1)\nplot(n,x)\n");
			fprintf(f, "xlim([0 %d])\n", N-1);
			fprintf(f, "subplot(3,1,2)\nplot(n,Xre,n,Xim)\n");
			fprintf(f, "xlim([0 %d])\n", N/2);
			fprintf(f, "subplot(3,1,3)\nstem(n,P)\n");
			fprintf(f, "xlim([0 %d])\n", N/2);
			fclose(f);
		}

		//u[0] = 0;
		//u[N-1] = 0;

        // calculate forces and next step over realisations
		ddu[0] = - pair_force(u[1] - u[0], a, C) + pair_force(u[0] - u[N-1], a, C) ;  //periodic left
		//ddu[0] = - pair_force(u[1] - u[0], a, C);
		//forces
		for(int j=1; j < N-1; j++) {
			ddu[j] = - pair_force(u[j+1]- u[j], a, C)  + pair_force(u[j] - u[j-1], a, C) ;
		}
		
		ddu[N-1] = - pair_force(u[0] - u[N-1], a, C)  + pair_force(u[N-1] - u[N-2], a, C); //periodic right
		//ddu[N-1] = - pair_force(u[N-1] - u[N-2], a, C); //periodic right

		//integration
		for(int j=0; j<N; j++){
			du[j] += dt*ddu[j];
			u[j] += dt*du[j];
        }

    }

    
    return 0;
}

