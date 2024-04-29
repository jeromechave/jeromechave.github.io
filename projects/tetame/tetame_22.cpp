/*****************************************************
    Estimation of neutral parameters by maximum likelihood

    J. Chave and F. Jabot
    last update 05-23-2008
******************************************************/

// compilation on cygwin : g++ -O3 -o tetame_v2 tetame_version2.cpp -mno-cygwin


// Option for use in a parallel processor (0 if off, 1 otherwise)
//#define MPI0 0

//#define  INTEL_EXTENDED_IEEE_FORMAT   0 

// Libraries
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <math.h>
//#include <malloc.h>
#include <time.h>    
using namespace std;



#if MPI0
#include <mpi.h>
#endif



// Routines called in the program
// Simplex method for ML estimation
long double simplex(long double (*func)(long double[]), long double start[],int n, long double EPSILON, long double scale, long double LEWENS, long double EWENS_THETA);
long double llik(long double x[]);

// Search of the max of two numbers
inline long double max(long double x, long double y){                            
    if (x > y) return x;
    else return y;
}
// Search of the min of two numbers
inline long double mini(long double x, long double y){                            
    if (x > y) return y;
    else return x;
}
// The Pochhammer symbol as defined by the difference of two Gamma 
//functions
inline long double lpochham(long double x, int n){
    return lgammal(x+n)-lgammal(x);
}

// Function for the quick sort routine
static  int intcompare2(const void *i, const void *j)
{
  int *u,*v;
  u = (int *) i;
  v = (int *) j;
  if (*u > *v)
    return (1);
  if (*u < *v)
    return (-1);
  return (0);
}

// Global variables
long double *K;                                   // K[A]=K(D,A)in Etienne's paper
long double factor;
long double J,SPP;
int **Abund;
long double *xi;
long double SPP2,JM,correc;


#define MAX_IT      1000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */

long double llik(long double x[]){                 // loglikelihood function where x[0]=theta, x[1]=I,K[A]=log(K(D,A)), cf Etienne, 2005 
  long double A;
  long double summand0, summand1, lsummand;
  summand0=0.0;
  summand1=0.0;
  long double divisor=0.0;
  long double newdivisor=0.0;
  for(A=SPP;A<=J;A++) {
  
    lsummand = factor+logl(x[0])*SPP-lpochham(x[1],int(J))+K[int(A)]+A*logl(x[1])-lpochham(x[0],int(A))-divisor;
    if (lsummand>11300){
        newdivisor=(lsummand-11300);
        divisor +=newdivisor;
        summand0=(summand0/expl(newdivisor))+expl(11300);
    }
    else {
    if ((lsummand>-11333.2)&&(summand0<expl(11330))){
        summand0+=expl(lsummand);         
    }
    else {
        if (summand0>expl(11330)){
            divisor+=1;
            summand0=(summand0/expl(1))+expl(lsummand-1);
        }
        else{ //NEW in version 2.1
            if (summand1==0){
                summand1=lsummand;
            }
            else{
                summand1+=logl(1+(expl(lsummand-summand1)));
            }
        }
    }           
  }
  }
  
    if (summand0>0){
      return -logl(summand0)-4500.0*logl(10)-divisor;
    }
    else{//NEW in version 2.1
        return -summand1-4500.0*logl(10);
    }
 

}

long double llikm(long double x[], int samp){                 // loglikelihood function where x[0]=theta, x[1]=I,K[A]=log(K(D,A)), cf Etienne, 2005 
  long double A;
  long double summand0, lsummand;
  summand0=0.0;
  long double divisor=0.0;
  long double newdivisor=0.0;
  long double product=0.0;
  for(int i=0;i<SPP;i++) {
    if (Abund[samp][i]>0){
        product+=lpochham((x[1]*xi[i]),Abund[samp][i]);
        }
  }
  lsummand = factor-lpochham(x[1],int(J))+product;
  return -lsummand;
 

}


long double lliklog(long double x[]){                 // loglikelihood function where x[0]=theta, x[1]=log(m),K[A]=log(K(D,A)), cf Etienne, 2005 (here, I=expl(x[1])*(J-1)/(1-expl(x[1])))
  long double A;
  long double summand0, summand1, lsummand;
  summand0=0.0;
  summand1=0.0;
  long double divisor=0.0;
  long double newdivisor=0.0;
  for(A=SPP;A<=J;A++) {
  
    lsummand = factor+logl(x[0])*SPP-lpochham(expl(x[1])*(J-1)/(1-expl(x[1])),int(J))+K[int(A)]+A*logl(expl(x[1])*(J-1)/(1-expl(x[1])))-lpochham(x[0],int(A))-divisor;
    if (lsummand>11300){
        newdivisor=(lsummand-11300);
        divisor +=newdivisor;
        summand0=(summand0/expl(newdivisor))+expl(11300);
    }
    else {
    if ((lsummand>-11333.2)&&(summand0<expl(11330))){
        summand0+=expl(lsummand);            
    }
    else {
        if (summand0>expl(11330)){
            divisor+=1;
            summand0=(summand0/expl(1))+expl(lsummand-1);
        }
        else{//NEW in version 2.1
            if (summand1==0){
                summand1=lsummand;
            }
            else{
                summand1+=logl(1+(expl(lsummand-summand1)));
            }
        }
    }           
  }
  }
    if (summand0>0){
      return -logl(summand0)-4500.0*logl(10)-divisor;
    }
    else {//NEW in version 2.1
      return -summand1-4500.0*logl(10);
    }
}

//search of the minimum of -loglikelihood using a simplex method
long double simplex(long double (*func)(long double[]), long double start[],int n, long double EPSILON, long double scale, long double LEWENS, long double EWENS_THETA){

  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  int i,j,m,row;
  int k;      /* track the number of function evaluations */
  int itr;    /* track the number of iterations */
  
  long double **v;     /* holds vertices of simplex */
  long double pn,qn;   /* values used to create initial simplex */
  long double *f;      /* value of function at each vertex */
  long double fr;      /* value of function at reflection point */
  long double fe;      /* value of function at expansion point */
  long double fc;      /* value of function at contraction point */
  long double *vr;     /* reflection - coordinates */
  long double *ve;     /* expansion - coordinates */
  long double *vc;     /* contraction - coordinates */
  long double *vm;     /* centroid - coordinates */
  long double min;
  
  long double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  v =  new long double*[n+1];
  f =  new long double[n+1];
  vr = new long double[n];
  ve = new long double[n];  
  vc = new long double[n];  
  vm = new long double[n];  
  
  /* allocate the columns of the arrays */
  for (i=0;i<=n;i++) {
    v[i] = new long double[n];
  }
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));             
  qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));
  
  for (i=0;i<2;i++) {                
    v[0][i] = start[i];
  }
  
  for (i=1;i<=n;i++) {
    for (j=0;j<2;j++) {
      if (i-1 == j) {
    v[i][j] = pn + start[j];
      }
      else {
    v[i][j] = qn + start[j];
      }
    }
  }
  
  /* find the initial function values */
  for (j=0;j<=n;j++) {
    f[j] = func(v[j]);
  }
  
  k = n+1;
  
  /* begin the main loop of the minimization */
  for (itr=1;itr<=MAX_IT;itr++) {     
    /* find the index of the largest value */
    vg=0;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vg]) {
    vg = j;
      }
    }

    /* find the index of the smallest value */
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
    vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
    vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
    if (m!=vg) {
      cent += v[m][j];
    }
      }
      vm[j] = cent/n;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0;j<=n-1;j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = max(vm[j]+ALPHA*(vm[j]-v[vg][j]),0.01);
    }
    fr = func(vr);
    k++;
    
    if (fr < f[vh] && fr >= f[vs]) {  
      for (j=0;j<=n-1;j++) {
      v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0;j<=n-1;j++) {
      ve[j] = max(vm[j]+GAMMA*(vr[j]-vm[j]),0.01);
      }
      fe = func(ve);
      k++;

    /* by making fe < fr as opposed to fe < f[vs],               
     Rosenbrocks function takes 63 iterations as opposed 
     to 64 when using long double variables. */
      
        if (fe < fr) {
            for (j=0;j<=n-1;j++) {
            v[vg][j] = ve[j];
            }
            f[vg] = fe;
        }
        else {
            for (j=0;j<=n-1;j++) {
            v[vg][j] = vr[j];
            }
        f[vg] = fr;
        }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0;j<=n-1;j++) {
            vc[j] = max(vm[j]+BETA*(vr[j]-vm[j]),0.01);
        }
        fc = func(vc);
        k++;
      }
      else {
        /* perform inside contraction */
        for (j=0;j<=n-1;j++) {
            vc[j] = max(vm[j]-BETA*(vm[j]-v[vg][j]),0.01);
        }
        fc = func(vc);
        k++;
      }

      
      if (fc < f[vg]) {
        for (j=0;j<=n-1;j++) {
            v[vg][j] = vc[j];
        }
        f[vg] = fc;
      }
      
      /* at this point the contraction is not successful,
     we must halve the distance from vs to all the 
     vertices of the simplex and then continue.
     10/31/97 - modified to account for ALL vertices.              
      */
      else {
        for (row=0;row<=n;row++) {
            if (row != vs) {
                for (j=0;j<=n-1;j++) {
                    v[row][j] = max(v[vs][j]+(v[row][j]-v[vs][j])/2.0,0.01);
                }
            }
        }
        f[vg] = func(v[vg]);
        k++;
        f[vh] = func(v[vh]);
        k++;    
      }
    }
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0;j<=n;j++) {
      fsum += f[j];
    }
    favg = fsum/(n+1);
    s = 0.0;
    for (j=0;j<=n;j++) {
      s += powl((f[j]-favg),2.0)/(n);
    }
    s = sqrt(s);
    if (s < EPSILON) break;
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  for (j=0;j<2;j++) {//NEW in version 2.1 : n replaced by 2
    start[j] = v[vs][j];
  }
  min=f[vs];                                           
  
  delete[] f;
  delete[] vr;
  delete[] ve;
  delete[] vc;
  delete[] vm;
  for (i=0;i<=n;i++) {
	  delete[] v[i];
  }
  delete[] v;
  return min;
}

long double simplexm(long double (*func)(long double[], int), long double start[],int n, long double EPSILON, long double scale, int samp){

  int vs;         /* vertex with smallest value */
  int vh;         /* vertex with next smallest value */
  int vg;         /* vertex with largest value */
  
  int i,j,m,row;
  int k;      /* track the number of function evaluations */
  int itr;    /* track the number of iterations */
  
  long double **v;     /* holds vertices of simplex */
  long double pn,qn;   /* values used to create initial simplex */
  long double *f;      /* value of function at each vertex */
  long double fr;      /* value of function at reflection point */
  long double fe;      /* value of function at expansion point */
  long double fc;      /* value of function at contraction point */
  long double *vr;     /* reflection - coordinates */
  long double *ve;     /* expansion - coordinates */
  long double *vc;     /* contraction - coordinates */
  long double *vm;     /* centroid - coordinates */
  long double min;
  
  long double fsum,favg,s,cent;
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */
  v =  new long double*[n+1];
  f =  new long double[n+1];
  vr = new long double[n];
  ve = new long double[n];  
  vc = new long double[n];  
  vm = new long double[n];  
  
  /* allocate the columns of the arrays */
  for (i=0;i<=n;i++) {
    v[i] = new long double[n];
  }
  
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));             
  qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));
  
  for (i=0;i<2;i++) {                
    v[0][i] = start[i];
  }
  
  for (i=1;i<=n;i++) {
    for (j=0;j<2;j++) {
      if (i-1 == j) {
    v[i][j] = pn + start[j];
      }
      else {
    v[i][j] = qn + start[j];
      }
    }
  }
  
  /* find the initial function values */
  for (j=0;j<=n;j++) {
    f[j] = func(v[j],samp);
  }
  
  k = n+1;
  
  /* begin the main loop of the minimization */
  for (itr=1;itr<=MAX_IT;itr++) {     
    /* find the index of the largest value */
    vg=0;
    //cerr <<"itr="<< itr<<"\t";
    for (j=0;j<=n;j++) {
    //cerr <<v[j][0]<<"\t"<<v[j][1]<<"\t"<<f[j]<<endl;
      if (f[j] > f[vg]) {
    vg = j;
      }
    }

    /* find the index of the smallest value */
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
    vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
    vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
    if (m!=vg) {
      cent += v[m][j];
    }
      }
      vm[j] = cent/n;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0;j<=n-1;j++) {
      /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
      vr[j] = max(vm[j]+ALPHA*(vm[j]-v[vg][j]),0.01);
    }
    fr = func(vr,samp);
    //cerr <<vr[0]<<"\t"<<vr[1]<< "fr="<<fr;
    k++;
    
    if (fr < f[vh] && fr >= f[vs]) {  
      for (j=0;j<=n-1;j++) {
      v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    if ( fr <  f[vs]) {
      for (j=0;j<=n-1;j++) {
      ve[j] = max(vm[j]+GAMMA*(vr[j]-vm[j]),0.01);
      }
      fe = func(ve,samp);
      k++;

    /* by making fe < fr as opposed to fe < f[vs],               
     Rosenbrocks function takes 63 iterations as opposed 
     to 64 when using long double variables. */
      
        if (fe < fr) {
            for (j=0;j<=n-1;j++) {
            v[vg][j] = ve[j];
            }
            f[vg] = fe;
            //cerr << "case1";
        }
        else {
            for (j=0;j<=n-1;j++) {
            v[vg][j] = vr[j];
            }
        f[vg] = fr;
        //cerr << "case2";
        }
    }
    
    /* check to see if a contraction is necessary */
    if (fr >= f[vh]) {
      if (fr < f[vg] && fr >= f[vh]) {
        /* perform outside contraction */
        for (j=0;j<=n-1;j++) {
            vc[j] = max(vm[j]+BETA*(vr[j]-vm[j]),0.01);
        }
        fc = func(vc,samp);
        k++;
        //cerr << "case3";
      }
      else {
        /* perform inside contraction */
        for (j=0;j<=n-1;j++) {
            vc[j] = max(vm[j]-BETA*(vm[j]-v[vg][j]),0.01);
        }
        fc = func(vc,samp);
        k++;
        //cerr << "case4";
      }

      
      if (fc < f[vg]) {
        for (j=0;j<=n-1;j++) {
            v[vg][j] = vc[j];
        }
        f[vg] = fc;
        //cerr << "case5";
      }
      
      /* at this point the contraction is not successful,
     we must halve the distance from vs to all the 
     vertices of the simplex and then continue.
     10/31/97 - modified to account for ALL vertices.              
      */
      else {
        for (row=0;row<=n;row++) {
            if (row != vs) {
                for (j=0;j<=n-1;j++) {
                    v[row][j] = max(v[vs][j]+(v[row][j]-v[vs][j])/2.0,0.01);
                }
            }
        }
        f[vg] = func(v[vg],samp);
        k++;
        f[vh] = func(v[vh],samp);
        k++;
        //cerr << "case6";

    
      }
    }
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0;j<=n;j++) {
      fsum += f[j];
    }
    favg = fsum/(n+1);
    //cerr <<"favg="<<favg<<endl;
    s = 0.0;
    for (j=0;j<=n;j++) {
      s += powl((f[j]-favg),2.0)/(n);
    }
    s = sqrt(s);
    //cerr <<"s="<<s<<endl;
    if (s < EPSILON) break;
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }
  
  for (j=0;j<2;j++) { //n replaced by 2 in VERSION 2.1
    start[j] = v[vs][j];
  }
  min=f[vs];                                           
  //cout << k << " Function evaluations\n";
  //cout << itr << " Optimizations\n";
  
  delete[] f;
  delete[] vr;
  delete[] ve;
  delete[] vc;
  delete[] vm;
  for (i=0;i<=n;i++) {
	  delete[] v[i];
  }
  delete[] v;
  return min;
}


int main() {

cerr << "ESTIMATING NEUTRAL PARAMETERS BY MAXIMUM LIKELIHOOD"<<endl;
cerr << "This program can be used for:"<<endl;
cerr<<"1-estimating theta and m of Hubbell's 2001 neutral theory using Etienne's 2005 method."<<endl;
cerr<<"2-estimating m in several samples belonging to the same regional pool using Jabot et al. 2008 method."<<endl;
cerr<<"For more details, see the manual."<<endl;
cerr << "Options for answering yes/no questions: 1, y, or yes for 'yes'; 0, n, or no for 'no'.\n \n";
// Parallel stuff
#if MPI0
    int rank,size;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    cerr << "Nb procs " << size << " rank " << rank << endl;
#endif

// Commands to read input files
    int i,n,im,s,sample;
    char nomfi[256], nomfo[256], nomfom[256],nomfobis[256],nomfogrille[256],nomfogrilleR[256];
int nb;
char bufi[128];
cerr << "Please enter the data file name (without'.txt') ";
cin >> bufi;


    sprintf(nomfi,"%s.txt",bufi);                                  
    sprintf(nomfo,"%s_out.txt",bufi); 
    sprintf(nomfom,"%s_outm.txt",bufi);       
    sprintf(nomfobis,"%s_out2.txt",bufi);    
    sprintf(nomfogrille,"%s_outg.txt",bufi);  
    sprintf(nomfogrilleR,"%s_outR.txt",bufi); 

#if MPI0
    sprintf(nomfi,"/users/p0509/chave/PARI/%s%d.txt",bufi,rank);
    sprintf(nomfo,"/users/p0509/chave/PARI/%s_out%d.txt",bufi,rank);
    sprintf(nomfobis,"/users/p0509/chave/PARI/%s_stat%d.txt",bufi,rank);
    sprintf(nomfogrille,"/users/p0509/chave/PARI/%s_statg%d.txt",bufi,rank);
    cerr << "Nom fichier entree " << nomfi << endl;
    cerr << nomfo << endl;
    MPI_Barrier(MPI_COMM_WORLD);
#endif

// First stage: count the species abundances, rank them
    ifstream inf;                                                         
    inf.open(nomfi);
    if(!inf){
        cerr<<"Failed to open data file\n";
        cerr<<"Please close the window and start again\n";
        int toto;
        cin>>toto;
    }
// Number of samples on which to find the MLE 
    int nbsamples;

/* The input file structure is 
        abundance_sample1_1
        abundance_sample1_2
        ...
        abundance_sample1_n
        &
        abundance_sample2_1
        ...
    The file is read twice, first to locate the parameter and 
    dimension the memory files, then to input the data
    Maximum: 10000 samples
*/
    cerr << "Input file: " << nomfi << endl;                           
    cerr << "Reading the file stats ...\n";
    int *Species0 = new int[10000];
    for(i=0;i<10000;i++)
        Species0[i]=0;
    i=0;
    char buffer[5];                                         
    while(inf.eof() == 0){                                     
        inf >> buffer;                                       
        if(buffer[0] == '&') i++;
        else Species0[i]++;
    }
    i++;                                                     
    inf.close();
    nbsamples = i;
    cerr << "Number of samples: " << nbsamples << endl;
    Species0[nbsamples-1]-=1;
    int *Species;
    Species = new int[nbsamples];
    
    ifstream inf2(nomfi);
    int *Speciesvrai;
    Speciesvrai = new int[nbsamples];
    for(sample=0;sample<nbsamples;sample++)
        Species[sample]=Species0[sample];                                 
    delete [] Species0;


    char met='e';
    if (nbsamples>1) {
        cerr<<endl;
        cerr<<"TeTame detected several samples in your data file."<<endl;
        cerr<<"If you want to do the simultaneous estimation of Theta and m, using Etienne (2005)'s likelihood formula, press e and enter"<<endl;
        cerr<<"If you want to estimate m in each sample, assuming that they belong to the same metacommunity, using Jabot et al.(2008)'s method, press j and enter"<<endl;
        cerr<<"In the case you want to do the multi-samples m inference, make sure you entered the pooled-over-the-samples species abundances in the beginning of the data file."<<endl;
        cin >> met;
    }
    
if (met=='e'){
    Abund = new int*[nbsamples];                          
    int **Abund2;
    Abund2 = new int*[nbsamples]; 

    
    for(sample=0;sample<nbsamples;sample++){
        Abund2[sample] = new int[Species[sample]];
        int espece=0;
        for(s=0;s<Species[sample];s++) {
            inf2 >> Abund2[sample][s];
            if (Abund2[sample][s]>0){
                espece++;
            }                          
        }
        Abund[sample] = new int[espece];
        espece=0;
        for(s=0;s<Species[sample];s++) {
            if (Abund2[sample][s]>0){
                Abund[sample][espece]=Abund2[sample][s];
                espece++;
            }                          
        }
        inf2 >> buffer;
        Species[sample]=espece;
        cerr << "In sample " << sample+1<< ", number of species: " << Species[sample] <<endl;
        qsort(Abund[sample],Species[sample],sizeof(int),intcompare2);  
    }
    inf2.close();
    ofstream out(nomfo);    
    out  <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";
    
    
  for(sample=0;sample<nbsamples;sample++){ //CORRECTION OF THE ERROR IN VERSION 2.0: COMPUTE K(D,A) for each sample!!!
    J=0;
    SPP = Species[sample];
    for(s=0;s<SPP;s++)
        J += Abund[sample][s];
     cerr << "Number of individuals: "<< J << endl;
   

    int MaxA = 0;
    MaxA=Abund[sample][(int)SPP-1];

    cerr <<" "<<endl;
    cerr << "Sample " << 1+sample << endl;
    cerr << "Maximal abundance: " << MaxA << endl;

    

    // abundance distribution
    int *Phi = new int[MaxA+1];
    for(s=0;s<=MaxA;s++) Phi[s]=0;
    for(s=0;s<SPP;s++) Phi[Abund[sample][s]]++;

    // Number of distinct abundances
    int NDA=0;
    for(s=0;s<=MaxA;s++) if(Phi[s] > 0) {NDA++;}
    
    
    cerr << "Start computing Stirling numbers ...\n";
    // FIRST STAGE: compute the Stirling numbers
    // I use the relationship S(n,m)=S(n-1,m-1)-(n-1)S(n-1,m)
    // where n >= m >= 1
    // The relation of interest is sum_m S(n,m)/S(n,1)*S(m,1) * x^m
    // defined to be equal to sum_m T(n,m) * x^m
    // The recurrence relation on T(n,m) is 
    // T(n,m)= T(n-1,m) + T(n-1,m-1)*(m-1)/(n-1)
   
    int *f = new int[NDA];
    int *g = new int[NDA];
    i=0;
    for(s=0;s<NDA;s++) {f[s]=0;g[s]=0;}
    for(n=0;n<=MaxA;n++) if(Phi[n] > 0) {             
        f[i] = Phi[n];                                  
        g[i] = n;                                        
        i++;
        }
    long double **T= new long double*[NDA];          // T(n,m) just for the n which are useful
    T[0] = new long double[g[0]+1];
    T[0][0]=0;T[0][1]=1;
    if (g[0]!=1){
        long double *lS2 = new long double[g[0]+1]; 
        lS2[0]=0;lS2[1]=1;
        for (n=2;n<=g[0];n++) {
            long double *lS1 = new long double[n+1];                
            for(im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
            delete[] lS1;            
        }
        for(im=2;im<=g[0];im++) {
            T[0][im]=lS2[im];
        }
        delete[] lS2;
    }
    for (int in=1;in<i;in++){
        T[in]= new long double[g[in]+1];
        T[in][0]=0;T[in][1]=1;
        long double *lS2 = new long double[g[in]+1];         
        for(im=0;im<=g[in-1];im++) {
                lS2[im] = T[in-1][im];
            }
        for (n=g[in-1]+1;n<=g[in];n++) {
            long double *lS1 = new long double[n+1];                
            for(im=0;im<=n-1;im++) {
                lS1[im] = lS2[im];
            }
            lS1[n]=0;
            for(im=2;im<=n;im++) {
                lS2[im] = lS1[im]+lS1[im-1]*(im-1)/(n-1); 
            }
            delete[] lS1;            
        }
        for(im=2;im<=g[in];im++) {
            T[in][im]=lS2[im];
        }
        delete[] lS2;
    }
    // After this stage we have stored in T[i][m] T(g[i],m)
    // with T(n,m) = S(n,m)*S(m,1)/S(n,1) for i>0

    cerr << "Start computing ln(K(D,A)) ...\n";
    // SECOND STAGE: compute the K(D,A)
    // I follow Etienne's route. Compute the product of polynomials 
    // of length J
    int j,nn,mm;
    K = new long double[int(J)+1];
    long double *poly2 = new long double[int(J)+1];
    for(i=0;i<=J;i++){
        K[i] = poly2[i] = 0.0;
    }
    K[0]=1;
    int degree = 0;
    int spe=0;
    for(i=0;i<NDA;i++) // loop over number of distinct abundances
        for(j=0;j<f[i];j++){ // loop over abundances per class             
            for(nn=0;nn<=degree;nn++)
                for(mm=1;mm<=g[i];mm++){
                    if (K[nn]>0){                   
                       poly2[nn+mm] += T[i][mm]*K[nn];
                    }
                    
                }              
            degree += g[i];
            for(nn=0;nn<=degree;nn++){            
                K[nn] = (poly2[nn]/powl(10,(4500.0/SPP)));
                poly2[nn] = 0.0;
            }
            spe++;
        }
    
    for(i=int(SPP);i<=J;i++){
        K[i] = logl(K[i]);                                    // now K[A]=ln(K(D,A)/10^4500) in Etienne's paper
    }
    for(i=0;i<NDA;i++) delete[] T[i];
	delete[] T;
    delete[] poly2;
    delete[] f;
    delete[] g;

// search of "infinite" values in K[A]
    int borneinf=int(SPP-1);
    int bornesup=int(J+1);
    long double maxlog=11333.2;
    int infinity=0;
    for(i=int(SPP);i<=J;i++){
        if ((K[i]>maxlog)||(K[i]<-maxlog)) {
            infinity=1;
            break;
        }
        borneinf++;
    }  //after that, borneinf=indice next to infinity but before
    for(int i=0;i<=J-SPP;i++){
        if ((K[(int)J-i]>maxlog)||(K[(int)J-i]<-maxlog)) {
            infinity=1;
            break;
        }
        bornesup--;
    }    //after that, bornesup=indice next to infinity but after
if (infinity==1){
cerr << "WARNING : the sample is too large to compute an exact likelihood, the program is thus doing approximations. The following results are to be taken with caution"<<endl;
cerr << "Value of A above which K(D,A) is computed approximately ="<<borneinf<<endl;
cerr << "Value of A below which K(D,A) is computed approximately ="<<bornesup<<endl;

//fitting of the infinite values of K[A] by a polynom of degree 3
    //computing of the derivatives at the critic points
    long double Kprimeinf = K[borneinf]-K[borneinf-1];
    long double Kprimesup = K[bornesup+1]-K[bornesup];
    // definition of the parameters of the fitted polynom aX^3+bX^2+cX+d
    long double a,b,c,d;
    //inversion of the linear system of equations (with the Gauss method)
    long double borneinf2=(long double)borneinf*(long double)borneinf;
    long double borneinf3=(long double)borneinf2*(long double)borneinf;
    long double bornesup2=(long double)bornesup*(long double)bornesup;
    long double bornesup3=(long double)bornesup2*(long double)bornesup;
    d=(Kprimesup-3*bornesup2*K[borneinf]/borneinf3+(2*bornesup/(long double)borneinf-3*bornesup2/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf)))/((6*bornesup2/borneinf3)-(6*bornesup/borneinf2)-((1+3*bornesup2/borneinf2-4*bornesup/(long double)borneinf)/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2))*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3));
    c=((K[bornesup]-bornesup3*K[borneinf]/borneinf3+(bornesup2/(long double)borneinf-bornesup3/borneinf2)*(Kprimeinf-3*K[borneinf]/(long double)borneinf))-d*(1-3*bornesup2/borneinf2+2*bornesup3/borneinf3))/(bornesup-2*bornesup2/(long double)borneinf+bornesup3/borneinf2);
    b=(Kprimeinf-3*K[borneinf]/(long double)borneinf+2*c+3*d/(long double)borneinf)/(0.0-(long double)borneinf);
    a=(K[borneinf]-b*borneinf2-c*(long double)borneinf-d)/borneinf3;
    
    //reconstruction of K[A] with the fitted polynom
    for (int i=borneinf+1;i<bornesup;i++) {
     K[i]=(a*i*i*i+b*i*i+c*i+d);
  }
}

    // THIRD STEP: define the log-Likelihood
    // L(theta,I) = theta^S/(I)_J * sum_A K(D,A) I^A/(theta)_A
    // logL(theta,I) = S*log(theta)-log((I)_J)  
    //                  + log(sum_A K(D,A) I^A/(theta)_A)
    // where (x)_N = x(x+1)...(x+N-1)
    // I = m(J-1)/(1-m)
    //
    // theta is between 1 and S
    // m is between 0 and 1

    cerr << "Compute the Ewens theta and log-likelihood ...\n";
    long double sume,EWENS_THETA=1.0;

    for(int ii=0;ii<1000;ii++) {                                           //  initial value to lauch the simplex algorithm 
        sume=0.0;
        for(j=0;j<J;j++) sume +=1.0/(EWENS_THETA+j);
        EWENS_THETA = double(SPP)/sume;
    }
    
    long double LEWENS;
    factor = lgammal(J+1);                                      
    for(s=0;s<SPP;s++) factor -= logl(max(1,Abund[sample][s]));
    for(s=0;s<=MaxA;s++) {
      factor -= lgammal(Phi[s]+1);
    }
    delete[] Phi;

  LEWENS = lpochham(EWENS_THETA,int(J))-logl(EWENS_THETA)*SPP-factor;               //  LEWENS=-loglikelihood(EWENS_THETA) 
  cerr << "Ewens' -log-likelihood: " <<  LEWENS << endl;
  long double x[2];

#if MPI0
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  cerr << "Maximizing the likelihood ...\n";
  long double start[2] = {EWENS_THETA,J/10}; 
  long double startb[2] = {J/10,EWENS_THETA}; 
  
  char ini[128];
  int endquestion=0;
  while (endquestion==0){
  cerr << "Would you like to provide initial values for the optimization procedure? ";
  cin >> ini;
  if (!strcmp(ini,"1")||!strcmp(ini,"y")||!strcmp(ini,"yes")) {
    endquestion=1;
    cerr << "initial theta value: "<<endl;
    cin >> start[0];
    cerr << "initial m value: "<<endl;
    cin >> start[1];
    start[1]=start[1]*(J-1)/(1-start[1]);
    startb[0]=start[1];
    startb[1]=start[0];
  }
  else{
    if (!strcmp(ini,"0")||!strcmp(ini,"n")||!strcmp(ini,"no")) {
        endquestion=1;
    }
  }
  }
  int end = 0;
  while (end==0) {
  long double min,minb;
  //cerr<<"start "<<start[0]<<" "<<start[1]<<endl;
  min=simplex(&llik,start,5,1.0e-13,1,LEWENS,EWENS_THETA);  //launching of the simplex algorithm   
  end=1;         
  
  //cerr<<"startb "<<startb[0]<<" "<<startb[1]<<endl;
  minb=simplex(&llik,startb,5,1.0e-13,1,LEWENS,EWENS_THETA);  //launching of the simplex algorithm   
  
    if (minb>min){
        //out  <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";
        // computation of "standard deviations" of the parameters estimated by examination of the local curvature of the loglikelihood.
        long double fsp[2] = {start[0]*1.01,start[1]};
        long double fsm[2] = {start[0]*0.99,start[1]};
        long double D2theta_ll = (llik(fsp)+llik(fsm)-2.0*min)/(0.0001*start[0]*start[0]);
        fsp[0] = start[0];fsp[1]=start[1]*1.01;
        fsm[0] = start[0];fsm[1]=start[1]*0.99;
        long double D2I_ll = (llik(fsp)+llik(fsm)-2.0*min)/(0.0001*start[1]*start[1]);
        
        cout <<"\n RESULTS (also output in file named: "<<bufi<<"_out.txt):\n";
        cout <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";           
        cout << SPP << "\t" << J << "\t" << start[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llik(start) <<"\t"  <<EWENS_THETA << "\t" << LEWENS;
        out  << SPP << "\t" << J << "\t" << start[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llik(start) <<"\t"  <<EWENS_THETA << "\t" << LEWENS;
        
        fsp[0] = startb[0]*1.01;fsp[1]=startb[1];
        fsm[0] = startb[0]*0.99;fsm[1]=startb[1];
        D2theta_ll = (llik(fsp)+llik(fsm)-2.0*minb)/(0.0001*startb[0]*startb[0]);
        fsp[0] = startb[0];fsp[1]=startb[1]*1.01;
        fsm[0] = startb[0];fsm[1]=startb[1]*0.99;
        D2I_ll = (llik(fsp)+llik(fsm)-2.0*minb)/(0.0001*startb[1]*startb[1]);
        
        cout <<"\t" << startb[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< startb[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << startb[1]/(startb[1]+J-1)<<"\t"<<((J-1)/((startb[1]+J-1)*(startb[1]+J-1)))*1.0/sqrt(D2I_ll)<<"\t"<<minb<<endl;
        out <<"\t" << startb[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< startb[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << startb[1]/(startb[1]+J-1)<<"\t"<<((J-1)/((startb[1]+J-1)*(startb[1]+J-1)))*1.0/sqrt(D2I_ll)<<"\t"<<minb<<endl;
        cout <<""<<endl;
        
        out.flush();
        
    }
    else {
        //out  <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";
        // computation of "standard deviations" of the parameters estimated by examination of the local curvature of the loglikelihood.
        long double fsp[2] = {startb[0]*1.01,startb[1]};
        long double fsm[2] = {startb[0]*0.99,startb[1]};
        long double D2theta_ll = (llik(fsp)+llik(fsm)-2.0*minb)/(0.0001*startb[0]*startb[0]);
        fsp[0] = startb[0];fsp[1]=startb[1]*1.01;
        fsm[0] = startb[0];fsm[1]=startb[1]*0.99;
        long double D2I_ll = (llik(fsp)+llik(fsm)-2.0*minb)/(0.0001*startb[1]*startb[1]);
        
        cout <<"\n RESULTS (also output in file named: "<<bufi<<"_out.txt):\n";
        cout <<"S\t J\t Theta\t Std_Theta\t I\t Std_I\t m\t Std_m\t loglike_min\t Theta_Ewens\t loglike_Ewens\t Theta2\t Std_Theta2\t I2\t Std_I2\t m2\t Std_m2\t loglike_min2\n";           
        cout << SPP << "\t" << J << "\t" << startb[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< startb[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << startb[1]/(startb[1]+J-1)<<"\t"<<((J-1)/((startb[1]+J-1)*(startb[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llik(startb) <<"\t"  <<EWENS_THETA << "\t" << LEWENS;
        out  << SPP << "\t" << J << "\t" << startb[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< startb[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << startb[1]/(startb[1]+J-1)<<"\t"<<((J-1)/((startb[1]+J-1)*(startb[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llik(startb) <<"\t"  <<EWENS_THETA << "\t" << LEWENS;
        
        fsp[0] = start[0]*1.01;fsp[1]=start[1];
        fsm[0] = start[0]*0.99;fsm[1]=start[1];
        D2theta_ll = (llik(fsp)+llik(fsm)-2.0*min)/(0.0001*start[0]*start[0]);
        fsp[0] = start[0];fsp[1]=start[1]*1.01;
        fsm[0] = start[0];fsm[1]=start[1]*0.99;
        D2I_ll = (llik(fsp)+llik(fsm)-2.0*min)/(0.0001*start[1]*start[1]);
        
        cout <<"\t" << start[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<<"\t"<<min<<endl;
        out <<"\t" << start[0] <<"\t" << 1.0/sqrt(D2theta_ll) <<"\t"<< start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<<"\t"<<min<<endl;
        cout <<""<<endl;
        
        out.flush();
    }
  
  
  
  
    
    int yes;
    int endquestion=0;
    while (endquestion==0){
        cerr << "Would you like to plot the likelihood surface? ";
        cin >> ini;
        if (!strcmp(ini,"1")||!strcmp(ini,"y")||!strcmp(ini,"yes")) {
            endquestion=1;
            yes=1;
        }
        else{
            if (!strcmp(ini,"0")||!strcmp(ini,"n")||!strcmp(ini,"no")) {
                endquestion=1;
                yes=0;
            }
        }
    }
    if (yes==1) {
    ofstream outb(nomfobis);  
    ofstream outg(nomfogrille); 
    ofstream outgR(nomfogrilleR);
    int points2;
    cerr << "How many grid points would you like to have? ";
    cin >> points2;
    int testlog2;
    int testlog3;
    int endquestion=0;
    while (endquestion==0){
        cerr << "Would you like to use a log scale for m? ";
        cin >> ini;
        if (!strcmp(ini,"1")||!strcmp(ini,"y")||!strcmp(ini,"yes")) {
            endquestion=1;
            testlog2=1;
        }
        else{
            if (!strcmp(ini,"0")||!strcmp(ini,"n")||!strcmp(ini,"no")) {
                endquestion=1;
                testlog2=0;
            }
        }
    }
    endquestion=0;
    while (endquestion==0){
        cerr << "Would you like to use a log scale for Theta? ";
        cin >> ini;
        if (!strcmp(ini,"1")||!strcmp(ini,"y")||!strcmp(ini,"yes")) {
            endquestion=1;
            testlog3=1;
        }
        else{
            if (!strcmp(ini,"0")||!strcmp(ini,"n")||!strcmp(ini,"no")) {
                endquestion=1;
                testlog3=0;
            }
        }
    }
     

        // construction of the likelihood plot by computing the likelihood of points on a grid
        long double thetamax;
        long double mmax;
        long double thetamin;
        long double mmin;
        cerr << "Value of theta_min (must be a positive number smaller than "<<start[0]<<"): ";
        cin >> thetamin;
        cerr << "Value of theta_max (must be a positive number greater than "<< start[0]<<"): ";
        cin >> thetamax;
        cerr << "Value of m_min (must be a number between 0 and "<<start[1]/(start[1]+J-1)<<": ";
        cin >> mmin;
        cerr << "Value of m_max (must be a number between "<<start[1]/(start[1]+J-1)<<" and 1): ";
        cin >> mmax;
        cerr << "Computing the likelihood of the grid points ..."<<endl;
   
        long double deep=10.0;
        int points = int(ceil(sqrt((double)points2)));
        long double **table = new long double *[(points+1)*(points+1)+1];
        
        // Computing of the position of the maximum_likelihood in the grid
        int im,im2,it,it2;
        
        
    if (testlog2==0) { 
      if (testlog3==0){
        it= (int) floor(points*(start[0]-thetamin)/(thetamax-thetamin));
        it2=it;
        im= (int) floor(points*((start[1]/(J-1+start[1]))-mmin)/(mmax-mmin));
        im2=im; 
        table[0]=new long double[3];
        // initialization of the grid
        table[0][0]=start[0];table[0][1]=start[1]/(J-1+start[1]);table[0][2]=llik(start);
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                table[(t-1)*(points+1)+t2]=new long double[3];
                table[(t-1)*(points+1)+t2][0]=(long double) thetamin+ (thetamax-thetamin)*((double)(t-1)/(double)points);
                table[(t-1)*(points+1)+t2][1]=(long double) mmin+ (mmax-mmin)*((double)(t2-1)/(double)points);
                table[(t-1)*(points+1)+t2][2]=min+deep;
            }
        }  
        int y=0;
        long double start2[2];
        outb << "theta"<<"\t"<<"m"<<"\t"<<"llik"<<endl;
        outb << table[0][0]<<"\t"<<table[0][1]<<"\t"<<table[0][2]<<endl;
        // Creation of the R file
        outgR << "x<-c(";
        
        for(int t=0;t<points;t++){
            outgR <<(long double) thetamin+ (thetamax-thetamin)*((double)(t)/(double)points)<<",";
        }
        outgR <<(long double) thetamax;
        outgR << ")"<<endl;
        outgR << "y<-c(";
        for(int t=0;t<points;t++){
            outgR <<(long double) mmin+ (mmax-mmin)*((double)(t)/(double)points)<<",";
        }
        outgR << mmax;
        outgR << ")"<<endl;
        outgR << "x=sort(x)"<<endl;
        outgR << "y=sort(y)"<<endl;
        outgR << "data = read.table('"<<bufi<<"_outg.txt',header=FALSE)"<<endl;
        outgR << "datamatrix=as.matrix(data)"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]+deep<<","<<table[0][2]+500<<"),col=rainbow(100,start=0.65,end=0.7),main='"<<bufi<<"',xlab='Theta',ylab='m')"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]<<","<<table[0][2]+deep<<"),col=rainbow(100,start=0,end=0.65),add=TRUE)"<<endl;
        outgR << "contour(x,y,datamatrix,nlevels=1,levels=c("<<table[0][2]+3<<"),add=TRUE)"<<endl;

        outgR.flush();
        
        // Computing of the grid 
        
        for (it=0;it<points+1;it++){
            for (im=0;im<points+1;im++){
                start2[0]=table[it*(points+1)+im][0];
                start2[1]=(J-1)*table[it*(points+1)+im][1]/(1-table[it*(points+1)+im][1]);
                table[it*(points+1)+im][2]=llik(start2);
            }
        }
        
        //Printing of the grid
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                if (table[(t-1)*(points+1)+t2][2]<min+deep){
                    outg << table[(t-1)*(points+1)+t2][2]<<"\t";
                }
                else {
                    outg << min+deep<<"\t";
                }
            outb << table[(t-1)*(points+1)+t2][0]<<"\t"<<table[(t-1)*(points+1)+t2][1]<<"\t"<<table[(t-1)*(points+1)+t2][2]<<endl;
            }
            outg<<""<<endl;
        } 
      }
      else {
        it= (int) floor(points*(logl(start[0])-logl(thetamin))/logl(thetamax/thetamin));
        it2=it;
        im= (int) floor(points*((start[1]/(J-1+start[1]))-mmin)/(mmax-mmin));
        im2=im; 
        table[0]=new long double[3];
        // initialization of the grid
        table[0][0]=start[0];table[0][1]=start[1]/(J-1+start[1]);table[0][2]=llik(start);
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                table[(t-1)*(points+1)+t2]=new long double[3];
                table[(t-1)*(points+1)+t2][0]=(long double) expl(logl(thetamin)+ ((double)(t-1)/(double)points)*logl(thetamax/thetamin));
                table[(t-1)*(points+1)+t2][1]=(long double) mmin+ (mmax-mmin)*((double)(t2-1)/(double)points);
                table[(t-1)*(points+1)+t2][2]=min+deep;
            }
        }  
        int y=0;
        long double start2[2];
        outb << "theta"<<"\t"<<"m"<<"\t"<<"llik"<<endl;
        outb << table[0][0]<<"\t"<<table[0][1]<<"\t"<<table[0][2]<<endl;
        // Creation of the R file
        outgR << "x<-c(";
        
        for(int t=0;t<points;t++){
            outgR <<(long double) expl(logl(thetamin)+ ((double)(t)/(double)points)*logl(thetamax/thetamin))<<",";
        }
        outgR <<(long double) thetamax;
        outgR << ")"<<endl;
        outgR << "y<-c(";
        for(int t=0;t<points;t++){
            outgR <<(long double) mmin+ (mmax-mmin)*((double)(t)/(double)points)<<",";
        }
        outgR << mmax;
        outgR << ")"<<endl;
        outgR << "x=sort(x)"<<endl;
        outgR << "y=sort(y)"<<endl;
        outgR << "data = read.table('"<<bufi<<"_outg.txt',header=FALSE)"<<endl;
        outgR << "datamatrix=as.matrix(data)"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]+deep<<","<<table[0][2]+500<<"),col=rainbow(100,start=0.65,end=0.7),main='"<<bufi<<"',xlab='ln(Theta)',ylab='m')"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]<<","<<table[0][2]+deep<<"),col=rainbow(100,start=0,end=0.65),add=TRUE)"<<endl;
        outgR << "contour(x,y,datamatrix,nlevels=1,levels=c("<<table[0][2]+3<<"),add=TRUE)"<<endl;

        outgR.flush();
        
        // Computing of the grid 
        
        for (it=0;it<points+1;it++){
            for (im=0;im<points+1;im++){
                start2[0]=table[it*(points+1)+im][0];
                start2[1]=(J-1)*table[it*(points+1)+im][1]/(1-table[it*(points+1)+im][1]);
                table[it*(points+1)+im][2]=llik(start2);
            }
        }
        
        //Printing of the grid
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                if (table[(t-1)*(points+1)+t2][2]<min+deep){
                    outg << table[(t-1)*(points+1)+t2][2]<<"\t";
                }
                else {
                    outg << min+deep<<"\t";
                }
            outb << table[(t-1)*(points+1)+t2][0]<<"\t"<<table[(t-1)*(points+1)+t2][1]<<"\t"<<table[(t-1)*(points+1)+t2][2]<<endl;
            }
            outg<<""<<endl;
        } 
      } 
    }
    else {
        if (testlog3==0){
        it= (int) floor(points*(start[0]-thetamin)/(thetamax-thetamin));
        it2=it;
        im= (int) floor(points*(logl(start[1]/(J-1+start[1]))-logl(mmin))/logl(mmax/mmin));
        im2=im; 
        table[0]=new long double[3];
        // initialization of the grid
        table[0][0]=start[0];table[0][1]=logl(start[1]/(J-1+start[1]));
        start[1]=table[0][1]; table[0][2]=lliklog(start);
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                table[(t-1)*(points+1)+t2]=new long double[3];
                table[(t-1)*(points+1)+t2][0]=(long double) thetamin+ (thetamax-thetamin)*((double)(t-1)/(double)points);
                table[(t-1)*(points+1)+t2][1]=(long double) logl(mmin)+ logl(mmax/mmin)*((double)(t2-1)/(double)points);
                table[(t-1)*(points+1)+t2][2]=min+deep;
            }
        }  
        int y=0;
        long double start2[2];
        outb << "theta"<<"\t"<<"m"<<"\t"<<"llik"<<endl;
        outb << table[0][0]<<"\t"<<table[0][1]<<"\t"<<table[0][2]<<endl;
        
        // Creation of the R file
        outgR << "x<-c(";
        
        for(int t=0;t<points;t++){
            outgR <<(long double) thetamin+ (thetamax-thetamin)*((double)(t)/(double)points)<<",";
        }
        outgR << thetamax;
        outgR << ")"<<endl;
        outgR << "y<-c(";
        for(int t=0;t<points;t++){
            outgR <<(long double) logl(mmin)+ logl(mmax/mmin)*((double)(t)/(double)points)<<",";
        }
        outgR << logl(mmax);
        outgR << ")"<<endl;
        outgR << "x=sort(x)"<<endl;
        outgR << "y=sort(y)"<<endl;
        outgR << "data = read.table('"<<bufi<<"_outg.txt',header=FALSE)"<<endl;
        outgR << "datamatrix=as.matrix(data)"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]+deep<<","<<table[0][2]+500<<"),col=rainbow(100,start=0.65,end=0.7),main='"<<bufi<<"',xlab='Theta',ylab='ln(m)')"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]<<","<<table[0][2]+deep<<"),col=rainbow(100,start=0,end=0.65),add=TRUE)"<<endl;
        outgR << "contour(x,y,datamatrix,nlevels=1,levels=c("<<table[0][2]+3<<"),add=TRUE)"<<endl;

        outgR.flush();
        
        // Computing of the grid 
        
        for (it=0;it<points+1;it++){
            for (im=0;im<points+1;im++){
                start2[0]=table[it*(points+1)+im][0];
                start2[1]=table[it*(points+1)+im][1];
                table[it*(points+1)+im][2]=lliklog(start2);
            }
        }
        
        //Printing of the grid
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                if (table[(t-1)*(points+1)+t2][2]<min+deep){
                    outg << table[(t-1)*(points+1)+t2][2]<<"\t";
                }
                else {
                    outg << min+deep<<"\t";
                }
            outb << table[(t-1)*(points+1)+t2][0]<<"\t"<<table[(t-1)*(points+1)+t2][1]<<"\t"<<table[(t-1)*(points+1)+t2][2]<<endl;
            }
            outg<<""<<endl;
        }  
    }
    else {
    it= (int) floor(points*(logl(start[0]/thetamin))/logl(thetamax/thetamin));
        it2=it;
        im= (int) floor(points*(logl(start[1]/(J-1+start[1]))-logl(mmin))/logl(mmax/mmin));
        im2=im; 
        table[0]=new long double[3];
        // initialization of the grid
        table[0][0]=logl(start[0]);table[0][1]=logl(start[1]/(J-1+start[1]));
        start[1]=table[0][1]; table[0][2]=lliklog(start);
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                table[(t-1)*(points+1)+t2]=new long double[3];
                table[(t-1)*(points+1)+t2][0]=(long double) logl(thetamin)+ logl(thetamax/thetamin)*((double)(t-1)/(double)points);
                table[(t-1)*(points+1)+t2][1]=(long double) logl(mmin)+ logl(mmax/mmin)*((double)(t2-1)/(double)points);
                table[(t-1)*(points+1)+t2][2]=min+deep;
            }
        }  
        int y=0;
        long double start2[2];
        outb << "theta"<<"\t"<<"m"<<"\t"<<"llik"<<endl;
        outb << table[0][0]<<"\t"<<table[0][1]<<"\t"<<table[0][2]<<endl;
        // Creation of the R file
        outgR << "x<-c(";
        
        for(int t=0;t<points;t++){
            outgR <<(long double) logl(thetamin)+ logl(thetamax/thetamin)*((double)(t)/(double)points)<<",";
        }
        outgR << logl(thetamax);
        outgR << ")"<<endl;
        outgR << "y<-c(";
        for(int t=0;t<points;t++){
            outgR <<(long double) logl(mmin)+ logl(mmax/mmin)*((double)(t)/(double)points)<<",";
        }
        outgR << logl(mmax);
        outgR << ")"<<endl;
        outgR << "x=sort(x)"<<endl;
        outgR << "y=sort(y)"<<endl;
        outgR << "data = read.table('"<<bufi<<"_outg.txt',header=FALSE)"<<endl;
        outgR << "datamatrix=as.matrix(data)"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]+deep<<","<<table[0][2]+500<<"),col=rainbow(100,start=0.65,end=0.7),main='"<<bufi<<"',xlab='ln(Theta)',ylab='ln(m)')"<<endl;
        outgR << "image(x,y,datamatrix,zlim=c("<<table[0][2]<<","<<table[0][2]+deep<<"),col=rainbow(100,start=0,end=0.65),add=TRUE)"<<endl;
        outgR << "contour(x,y,datamatrix,nlevels=1,levels=c("<<table[0][2]+3<<"),add=TRUE)"<<endl;

        outgR.flush();
        
        // Computing of the grid
        
        for (it=0;it<points+1;it++){
            for (im=0;im<points+1;im++){
                start2[0]=expl(table[it*(points+1)+im][0]);
                start2[1]=table[it*(points+1)+im][1];
                table[it*(points+1)+im][2]=lliklog(start2);
            }
        }
        
        //Printing of the grid
        for (int t=1; t<points+2 ; t++) {
            for (int t2=1; t2<points+2 ; t2++) {
                if (table[(t-1)*(points+1)+t2][2]<min+deep){
                    outg << table[(t-1)*(points+1)+t2][2]<<"\t";
                }
                else {
                    outg << min+deep<<"\t";
                }
            outb << table[(t-1)*(points+1)+t2][0]<<"\t"<<table[(t-1)*(points+1)+t2][1]<<"\t"<<table[(t-1)*(points+1)+t2][2]<<endl;
            }
            outg<<""<<endl;
        }
    }
    }
    
    long double minimum=table[0][2];        // procedure of control that the maximum-likelihood estimator is not a local maximum
    int tt=0;
    for (int t=1;t<(points+1)*(points+1)+1;t++){
        if (table[t][2]<minimum) {
            minimum=table[t][2];
            tt=t;
        }
    }
    int warning=0;
    if (minimum<table[0][2]) {
        cerr <<"WARNING !!! : the grid computing procedure found a point with a larger likelihood that our estimate -> the maximum-likelihood estimator is in fact a local maximum" <<endl;
        cerr <<"We advise you to restart the optimization procedure, if you agree push 1 and enter"<<endl;
        cin >> warning;
    }
    if (warning==1) {
        end=0;
        start[0]=table[tt][0];
        if (testlog2==0){
            start[1]=(J-1)*table[tt][1]/(1-table[tt][1]);;
        }
        else {
            start[1]=expl(x[1])*(J-1)/(1-expl(x[1]));
        }
    }
    outb.flush();
    outg.flush();
//}
}
}
}
}
else {
    if (met=='j'){
    Abund = new int*[(nbsamples+1)];                          
    Abund[0] = new int[Species[0]];
    for(s=0;s<Species[0];s++) {
        Abund[0][s]=0;
    }
    for(sample=0;sample<nbsamples;sample++){
        Abund[(sample+1)] = new int[Species[sample]];
        Speciesvrai[sample]=0;
        for(s=0;s<Species[sample];s++) {
            inf2 >> Abund[(sample+1)][s];
            Abund[0][s]+=Abund[(sample+1)][s];
            if (Abund[(sample+1)][s]>0){
                Speciesvrai[sample]++;
            }                                
        }
        inf2 >> buffer; 
    }
    inf2.close();
    for(sample=0;sample<nbsamples;sample++){
        cerr << "In sample " << (sample+1)<< ", number of species: " << Speciesvrai[sample] <<endl;  
    }
    inf2.close();
    ofstream out(nomfom);                                                            
    out  <<"S\t J\t I\t Std_I\t m\t Std_m\t loglike_min\n";
    
    // Computing of the relative abundances in the metacommunity
    JM=0;
    SPP = Species[0];
    for(s=0;s<SPP;s++){
        JM += Abund[0][s];
    }
    cerr << "Number of individuals in the pooled samples: "<< JM << endl;
    xi=new long double [(int)SPP];
    for(s=0;s<SPP;s++){
        xi[s]= (long double) Abund[0][s]/ (long double) JM;
        //cerr << s << "\t"<<xi[s]<<endl;
    }
    
for(sample=1;sample<(nbsamples+1);sample++){
    // Total number of individuals
    J=0;
    SPP = Species[(sample-1)];
    int MaxA = 0;
    for(s=0;s<SPP;s++){
            J += Abund[sample][s];
            if (Abund[sample][s]>MaxA){
                MaxA=Abund[sample][s];
            }
        }
    cerr<<"In sample "<<sample<<":"<<endl;
    cerr << "Number of individuals: "<< J << endl;
    cerr << "Maximal abundance: " << MaxA << endl;

    // abundance distribution
    int *Phi = new int[MaxA+1];
    for(s=0;s<=MaxA;s++) Phi[s]=0;
    for(s=0;s<SPP;s++) Phi[Abund[sample][s]]++;
    for(s=0;s<=MaxA;s++){
        //cerr <<"phi"<<Phi[s]<<endl;
    }
    SPP2=SPP-Phi[0];
    cerr << "Number of species present in the sample: "<<SPP2<<endl;
    // Number of distinct abundances
    int NDA=0;
    for(s=0;s<=MaxA;s++) if(Phi[s] > 0) {NDA++;}
    //cerr << "NbDistinctAbund: " << NDA << endl;

    // THIRD STEP: define the log-Likelihood
    // L(theta,I) = theta^S/(I)_J * sum_A K(D,A) I^A/(theta)_A
    // logL(theta,I) = S*log(theta)-log((I)_J)  
    //                  + log(sum_A K(D,A) I^A/(theta)_A)
    // where (x)_N = x(x+1)...(x+N-1)
    // I = m(J-1)/(1-m)
    //
    // theta is between 1 and S
    // m is between 0 and 1

    factor = lgammal(J+1);  
    for(s=0;s<SPP;s++) {
        if (Abund[sample][s]>0){
            factor -= lgammal(Abund[sample][s]+1);
        }
    }                                    
    delete[] Phi;
    long double x[2];
 
  cerr << "Maximizing the likelihood ...\n";
  long double start[2] = {50,J};
  char ini[128];
  int end = 0;
  while (end==0) {
  long double min,minb;
  //cerr <<llikm(start,sample)<<endl;
  min=simplexm(&llikm,start,5,1.0e-13,1,sample);  //launching of the simplex algorithm   
  end=1;         
  // computation of "standard deviations" of the parameters estimated by examination of the local curvature of the loglikelihood.
  long double fsp[2] = {start[0]*1.01,start[1]};
  long double fsm[2] = {start[0]*0.99,start[1]};
  long double D2theta_ll = (llikm(fsp,sample)+llikm(fsm,sample)-2.0*min)/(0.01*start[0]);
  fsp[0] = start[0];fsp[1]=start[1]*1.01;
  fsm[0] = start[0];fsm[1]=start[1]*0.99;

  long double D2I_ll = (llikm(fsp,sample)+llikm(fsm,sample)-2.0*min)/(0.01*start[1]);
  cout <<"\n RESULTS (also output in file named: "<<bufi<<"_outm.txt):\n";
  cout  <<"S\t J\t I\t Std_I\t m\t Std_m\t loglike_min\n";           
  cout  << SPP2 << "\t" << J  <<"\t"<< start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llikm(start,sample) << endl;
  cout <<""<<endl;
  out  << SPP2 << "\t" << J << "\t" << start[1] << "\t" << 1.0/sqrt(D2I_ll) << "\t" << start[1]/(start[1]+J-1)<<"\t"<<((J-1)/((start[1]+J-1)*(start[1]+J-1)))*1.0/sqrt(D2I_ll)<< "\t" << llikm(start,sample) <<endl;
  out.flush();
    }
    delete[] K;

}

}
else{
    cerr<<"Unrecognized answer, please close the program and start again.";
    int toto;
    cin >> toto;
}
}
#if MPI0
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 0;
  
  cerr << "End of program";
}
