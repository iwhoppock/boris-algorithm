#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <sys/time.h>
//#include <mpi.h>

#define number_of_particles 10
#define time_steps 20000
#define TeTp 0.5   	//~1 (temperature of electrons : temperature of protons)
#define TaTp 4.0    //slow solar wind (temperature of alphas : " )
#define beta 1.0 	//proton beta ~ (v_th / v_A)^2
#define MeMp (1./1836.) //(mass of electron : mass of proton)
#define BOX 10.0   //one side of the volume of the space centred around origin
#define vAc 3e-4   	//Alfvén velocity over speed of light
#define PI 4*atan(1)
#define dt 0.01

#define mass (1.0)
#define charge (1.0)
#define stdev (sqrt(0.5 * beta))
//#define stdev (sqrt(0.5 * beta * TaTp / mass))	//use with alphas
//#define stdev (sqrt(0.5 * beta * TeTp / MeMp))	//use with electrons


static inline double WTime(void){
    //timing the simulation for efficiency
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}

void Boris(double x[3], double v[3], double E[3], double B[3]){
    //Boris Rotation for non-relativistic particle tracing -- "It really is state of the art" - D. Verscharen
    double v_minus[3];
    double v_prime[3];
    double v_plus[3];
    double t[3];
    double s[3];
    double t_mag2;
    
    //t vector
    t[0] = charge / mass * B[0] * 0.5 * dt;
    t[1] = charge / mass * B[1] * 0.5 * dt;
    t[2] = charge / mass * B[2] * 0.5 * dt;
    
    //|t|^2
    t_mag2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
    
    //s vector
    s[0] = 2.0 * t[0] / (1.0 + t_mag2);
    s[1] = 2.0 * t[1] / (1.0 + t_mag2);
    s[2] = 2.0 * t[2] / (1.0 + t_mag2);

    //v minus
    v_minus[0] = v[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
    v_minus[1] = v[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
    v_minus[2] = v[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;
    
    //v prime
    v_prime[0] = v_minus[0] + ( v_minus[1] * t[2] - v_minus[2] * t[1]);
    v_prime[1] = v_minus[1] + (-v_minus[0] * t[2] + v_minus[2] * t[0]);
    v_prime[2] = v_minus[2] + ( v_minus[0] * t[1] - v_minus[1] * t[0]);
    
    //v plus:
    v_plus[0] = v_minus[0] + ( v_prime[1] * s[2] - v_prime[2] * s[1]);
    v_plus[1] = v_minus[1] + (-v_prime[0] * s[2] + v_prime[2] * s[0]);
    v_plus[2] = v_minus[2] + ( v_prime[0] * s[1] - v_prime[1] * s[0]);
    
    //final v_n+1/2
    v[0] = v_plus[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
    v[1] = v_plus[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
    v[2] = v_plus[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;

    //pusher
    x[0] += v[0] * dt;
    x[1] += v[1] * dt;
    x[2] += v[2] * dt;    
}

double gaussian_number(void){
    //obtain a gaussian number from two random numbers on a range of [0,1] inclusive.
    //this procedure CAN produce two gaussian numbers, but we only need one...c'est la vie
    //also, as we multiply by a standard deviation that is beta-dependent, this is really Maxwellian, aber es ist egal...
    double number1 = (double)rand() / (double)RAND_MAX ;
    double number2 = (double)rand() / (double)RAND_MAX ;
    return stdev * sqrt(-2. * log(number1)) * sin(2. * PI * number2);
}

void get_fields(double E[3], double B[3], double x[3]){
    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 1.;
}

void fill_vector(double X[3][number_of_particles], double x[3], int particle_number){
    //used for shared-memory OpenMP caculations
    x[0] = X[0][particle_number];
    x[1] = X[1][particle_number];
    x[2] = X[2][particle_number];
}

void update_matrix(double X[3][number_of_particles], double x[3], int particle_number){
    //used for shared memory OpenMP calculations
	X[0][particle_number] = x[0];
    X[1][particle_number] = x[1];
    X[2][particle_number] = x[2];
}

void fill_matrices(double X[3][number_of_particles], double V[3][number_of_particles]){
    //initialise the matrices for position and velocity 
    //the if statement simply allows one to input specific values for tracing a single particle
    if (number_of_particles > 1){
        for (int j = 0; j< 3; j++){ 
            for (int i = 0; i < number_of_particles; i++) {
                V[j][i] = gaussian_number();
                X[j][i] = BOX * ((double)rand() / (double)RAND_MAX) - (BOX / 2.);
            }
        }
    }
    else{
        V[0][0] = 0.;
        V[1][0] = 1.;
        V[2][0] = 0.;
        X[0][0] = -1.;
        X[1][0] = 0.;
        X[2][0] = 0.;
    }
}

int main(){
    //initialise position and velocity matrices
    double X[3][number_of_particles];
    double V[3][number_of_particles];

    //FILL MATRICES:
    fill_matrices(X,V);
      
    //Openfile
	FILE *fout1 = NULL;
	fout1 = fopen("particle_tracer.csv", "w");

	//TIMEING
    double time1 = WTime();

    //TIME LOOP:
    for (int nt = 0; nt < time_steps; nt++){

#pragma omp parallel for

        //PARTICLE LOOP:
        for (int particle = 0; particle < number_of_particles; particle++){
            double E[3];
            double B[3];
            double x[3]; 
            double v[3];
            fill_vector(X,x,particle);
            fill_vector(V,v,particle);
            get_fields(E, B, x);
            Boris(x, v, E, B);
            update_matrix(X,x,particle); 
            update_matrix(V,v,particle); 
            fprintf(fout1,"%d, %g, %g, %g, %g, %g, %g, %g\n",particle,nt*dt,x[0],x[1],x[2],v[0],v[1],v[2]); 
        }
    }

    fclose(fout1); 

    double time2 = WTime();
    printf("TIME = %g\n",time2-time1);

    return 0;
} 