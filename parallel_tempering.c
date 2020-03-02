#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <omp.h>
#include "util.h"

#define DEBUG 0 /* 0 turns debugging off, 1 turns debugging on */

#define CODE_NAME "Parallel Tempering for 2D Ferromagnetic Ising Model"
#define AUTHOR "Andrew Noble"
#define VERSION "v0.1"
#define PROF "Prof. Helmut Katzgraber"

# define NUM_THREADS 5
# define SPACE_DIM 2
# define LATTICE_COORD_NUM 4
# define SEED 1234

clock_t START, END;

/*----------------------- Get magnetization of spin configuration -----------*/
int get_m(const int N, const int s[]) 
{   
    int i, m = 0;

    for (i=0; i<N; i++) 
        m += s[i];

    return m;
}

/*----------------------- Get internal energy of spin configuration ---------*/
int get_energy(const int N, const int L, const int s[]) 
{   
    int i, nn, delta, energy = 0;

    for (i=0; i<N; i++) {
        
        delta = 0;

        if ((nn = i + 1) >= N) nn -= N;
        delta += s[nn];
        if ((nn = i - 1) < 0) nn += N;
        delta += s[nn];
        if ((nn = i + L) >= N) nn -= N;
        delta += s[nn];
        if ((nn = i - L) < 0) nn += N;
        delta += s[nn];

        energy += s[i] * delta;
    }

    return -energy / 2;
}

/*----------------------- Perform Replica Exchange Update -------------------*/
void exchange(const int N, const int L, const int sizeTemperatureSet, 
    int s[][N], const float deltaBetaSet[], int exchangeAcc[], 
    int tempToReplicaIndex[], const gsl_rng *rng)
{
    int i;
    int energy[sizeTemperatureSet];
    int tmp;
    double arg;

    for (i=0; i<sizeTemperatureSet; i++) 
        energy[i] = get_energy(N, L, s[tempToReplicaIndex[i]]);

    for (i=0; i<(sizeTemperatureSet-1); i++) {

        arg = deltaBetaSet[i] * (energy[i+1] - energy[i]);

        if (arg >= 0 | exp(arg) > gsl_rng_uniform(rng)) {
            tmp = tempToReplicaIndex[i];
            tempToReplicaIndex[i] = tempToReplicaIndex[i+1];
            tempToReplicaIndex[i+1] = tmp;

            exchangeAcc[i] += 1;
        }
    }
}

/*----------------------- Perform Sweep of Spin Updates ---------------------*/
/*--------------- Acknowledgement: ------------------------------------------*/
/*--------------- Routine adapted from Newman and Barkema Textbook-----------*/
/*---------------------------------------------------------------------------*/
void sweep(const int N, const int L, int s[], double acceptanceRatios[], 
    const gsl_rng *rng)
{
    int i, nn, ps, psIndex, delta;

    for (i=0; i<N; i++) {

        psIndex = N * gsl_rng_uniform(rng);
        ps = s[psIndex];
        delta = 0;

        if ((nn = psIndex + 1) >= N) nn -= N;
        delta += s[nn];
        if ((nn = psIndex - 1) < 0) nn += N;
        delta += s[nn];
        if ((nn = psIndex + L) >= N) nn -= N;
        delta += s[nn];
        if ((nn = psIndex - L) < 0) nn += N;
        delta += s[nn];

        delta *= ps;

        if (delta <= 0 | acceptanceRatios[delta] > gsl_rng_uniform(rng))
            s[psIndex] = -ps;
    }
}

int main()
{
    /*------------------------------------------------------------------------*/ 
    /*----------------------- Begin read configure file ----------------------*/
    /*------------------------------------------------------------------------*/ 
    FILE *fp;
    fp = fopen("config.txt", "r");

    char buf[BUFSIZ];

    int L;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d\n", &L);

    int numBurninSweeps;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d\n", &numBurninSweeps);
    
    int numSweeps;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d\n", &numSweeps);

    int numSweepsPerSample;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d\n", &numSweepsPerSample);

    int sizeTemperatureSet;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%d\n", &sizeTemperatureSet);

    float temperatureSet[sizeTemperatureSet];
    int i = 0;
    float a;
    fgets(buf, BUFSIZ, fp);
    fgets(buf, BUFSIZ, fp);
    while (i < sizeTemperatureSet) {
        fgets(buf, BUFSIZ, fp);
        sscanf(buf, "%f\n", &a);
        temperatureSet[i++] = a;
    }

    #if (DEBUG == 1)
        print_float_array(temperatureSet, sizeTemperatureSet);
        printf("\n");
    #endif
    /*-----------------------------------------------------------------------*/ 
    /*----------------------- End read configure file -----------------------*/
    /*-----------------------------------------------------------------------*/ 

    /*-----------------------------------------------------------------------*/ 
    /*----------------------- Begin write output ----------------------------*/ 
    /*-----------------------------------------------------------------------*/ 
    fp = fopen("output.txt","w");

    fprintf(fp, "# code name : %s\n", CODE_NAME);
    fprintf(fp, "# author : %s\n", AUTHOR);
    fprintf(fp, "# version : %s\n", VERSION);
    fprintf(fp, "# inspired by : %s\n\n", PROF);

    int N = L * L;
    fprintf(fp, "# system size : %d\n", L);
    fprintf(fp, "# number of spins : %d\n", N);
    fprintf(fp, "# space dimension : %d\n", SPACE_DIM);
    fprintf(fp, "# lattice coordination number : %d\n", LATTICE_COORD_NUM);

    pid_t pid = getpid();
    fprintf(fp, "# process ID : %d\n", pid);
    gethostname(buf, sizeof(buf));
    fprintf(fp, "# hostname : %s\n", buf);
    time_t t;
    time(&t);
    fprintf(fp, "# job start : %s\n", ctime(&t));

    fprintf(fp, "# temperature set\n");
    for (i = 0; i < sizeTemperatureSet; i++) 
        fprintf(fp, "# | %.2f\n", temperatureSet[i]);
    fprintf(fp, "\n");
    fprintf(fp, "# T MCS <energy> <m>\n");

    /*---------- Initialize deltaBeta for replica exchange updates ----------*/ 
    float deltaBetaSet[sizeTemperatureSet];
    for (i = 0; i < (sizeTemperatureSet-1); i++) 
        deltaBetaSet[i] = 1. / temperatureSet[i+1] - 1. / temperatureSet[i];

    /*---------- Initialize map from temp to replica index ------------------*/ 
    int tempToReplicaIndex[sizeTemperatureSet];
    for (i = 0; i < sizeTemperatureSet; i++) 
        tempToReplicaIndex[i] = i;

    /*---------- Initialize spins -------------------------------------------*/ 
    int s[sizeTemperatureSet][N];
    int j = 0;
    for (i = 0; i < sizeTemperatureSet; i++) {
        for (j = 0; j < N; j++) 
                s[i][j] = 1;
    }

    /*---------- Initialize acceptanceRatios for spin updates ---------------*/ 
    double acceptanceRatios[sizeTemperatureSet][LATTICE_COORD_NUM + 1];
    for (i = 0; i < sizeTemperatureSet; i++) {
        for (j=2; j<(LATTICE_COORD_NUM + 1); j+=2)
           acceptanceRatios[i][j] = exp(-2. * j / temperatureSet[i]);
    }

    /*---------- Initialize arrays to hold samples of observables -----------*/
    int energy[sizeTemperatureSet];
    int m[sizeTemperatureSet];

    /*---------- Initialize array to count accepted exchange updates --------*/
    int exchangeAcc[sizeTemperatureSet-1];
    for (i = 0; i < (sizeTemperatureSet-1); i++)
        exchangeAcc[i] = 0;

    /*---------- Initialize random number generators, one per thread --------*/
    gsl_rng **rng_array;
    rng_array = malloc(NUM_THREADS * sizeof(gsl_rng*));

    for (i = 0; i < NUM_THREADS; i++) {   
        rng_array[i] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng_array[i], (i + 1) * SEED);
    }

    /*---------- Burnin, threading over temperatures ------------------------*/
    for (i = 0; i < numBurninSweeps; i++) {

        #pragma omp parallel
        {
            #pragma omp for
            for (j = 0; j < sizeTemperatureSet; j++) {
                sweep(N, L, s[tempToReplicaIndex[j]], acceptanceRatios[j], 
                    rng_array[omp_get_thread_num()]);

                #if (DEBUG == 1)
                    printf("Burnin: Hello from thread = %d\n", 
                        omp_get_thread_num());
                #endif
            }
        }

        exchange(N, L, sizeTemperatureSet, s, deltaBetaSet, exchangeAcc, 
            tempToReplicaIndex, rng_array[0]);

    } 

    /*---------- Zero out burnin counts of accepted exchange updates --------*/
    for (i = 0; i < (sizeTemperatureSet-1); i++)
        exchangeAcc[i] = 0;

    /*---------- Collect samples following burnin ---------------------------*/
    float energy_float, abs_m_float;

    for (i = 0; i < numSweeps; i++) {

        #pragma omp parallel
        {
            #pragma omp for
            for (j = 0; j < sizeTemperatureSet; j++) {
                sweep(N, L, s[tempToReplicaIndex[j]], acceptanceRatios[j], 
                    rng_array[omp_get_thread_num()]);
                
                #if (DEBUG == 1)
                    printf("Hello from thread = %d\n", omp_get_thread_num());
                #endif
            }
        }

        exchange(N, L, sizeTemperatureSet, s, deltaBetaSet, exchangeAcc, 
            tempToReplicaIndex, rng_array[0]);

        /*---------- Sample statistics written to output here ---------------*/
        if ((i+1) % numSweepsPerSample == 0) {

            for (j=0; j<sizeTemperatureSet; j++) {
                energy[j] = get_energy(N, L, s[tempToReplicaIndex[j]]);
                m[j] = get_m(N, s[tempToReplicaIndex[j]]);

                energy_float = (float) energy[j];
                abs_m_float = (float) abs(m[j]);
                fprintf(fp, "%.4f %d %.4f %.4f\n", temperatureSet[j], i+1, 
                    energy_float / N, abs_m_float / N);
            }

            fprintf(fp, "\n");

            #if (DEBUG == 1) 
                print_array(energy, sizeTemperatureSet);
                print_array(m, sizeTemperatureSet);
            #endif
        }

    } 

    /*---------- Write counts of accepted exchange updates to output --------*/
    float exchangeAcc_float;
    fprintf(fp, "# T A(T)\n");
    for (i=0; i<(sizeTemperatureSet-1); i++) {
        exchangeAcc_float = (float) exchangeAcc[i];
        fprintf(fp, "%.4f %.4f\n", temperatureSet[i], 
            exchangeAcc_float / numSweeps);
    }

    time(&t);
    fprintf(fp, "\n# job end : %s", ctime(&t));
    /*-----------------------------------------------------------------------*/ 
    /*----------------------- End Write Output ------------------------------*/ 
    /*-----------------------------------------------------------------------*/ 

    return 0;
}
