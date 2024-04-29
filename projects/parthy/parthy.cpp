/*************************************************************** 
   Simulation of density-dependent communities
   
   Franck Jabot
   Last Modified 21/09/2009
   Pour compiler :
   g++ -O3 -o ABC_densite-dep ABC_densite-dep.cpp  -mno-cygwin
   *************************************************************/

/* LIBRARIES */
/*************/
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <cmath>
#include <string.h>
using namespace std;

//Definitions of random numbers generators until line 1633.

/************************************************
MULTINOMIAL RANDOM GENERATION
*************************************************/
#define RANDOMC_H


// Define 32 bit signed and unsigned integers.
// Change these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else    
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seed
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};    


void CRandomMersenne::Init0(uint32 seed) {
   // Detect computer architecture
   union {double f; uint32 i[2];} convert;
   convert.f = 1.0;
   if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
   else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
   else Architecture = NONIEEE;

   // Seed generator
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(uint32 seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}


void CRandomMersenne::RandomInitByArray(uint32 seeds[], int length) {
   // Seed by more than 32 bits
   int i, j, k;

   // Initialize
   Init0(19650218);

   if (length <= 0) return;

   // Randomize mt[] using whole seeds[] array
   i = 1;  j = 0;
   k = (MERS_N > length ? MERS_N : length);
   for (; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
      i++; j++;
      if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
      if (j >= length) j=0;}
   for (k = MERS_N-1; k; k--) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
      if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
   mt[0] = 0x80000000UL;  // MSB is 1; assuring non-zero initial array

   // Randomize some more
   mti = 0;
   for (int i = 0; i <= MERS_N; i++) BRandom();
}


uint32 CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32 y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32 UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {    
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }

   y = mt[mti++];

#if 1
   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;
#endif

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   union {double f; uint32 i[2];} convert;
   uint32 r = BRandom();               // Get 32 random bits
   // The fastest way to convert random bits to floating point is as follows:
   // Set the binary explonent of a floating point number to 1+bias and set
   // the mantissa to random bits. This will give a random number in the 
   // interval [1,2). Then subtract 1.0 to get a random number in the interval
   // [0,1). This procedure requires that we know how floating point numbers
   // are stored. The storing method is tested in function RandomInit and saved 
   // in the variable Architecture.

   // This shortcut allows the compiler to optimize away the following switch
   // statement for the most common architectures:
#if defined(_M_IX86) || defined(_M_X64) || defined(__LITTLE_ENDIAN__)
   Architecture = LITTLE_ENDIAN1;
#elif defined(__BIG_ENDIAN__)
   Architecture = BIG_ENDIAN1;
#endif

   switch (Architecture) {
   case LITTLE_ENDIAN1:
      convert.i[0] =  r << 20;
      convert.i[1] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case BIG_ENDIAN1:
      convert.i[1] =  r << 20;
      convert.i[0] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case NONIEEE: default: ;
   } 
   // This somewhat slower method works for all architectures, including 
   // non-IEEE floating point representation:
   return (double)r * (1./((double)(uint32)(-1L)+1.));
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((max - min + 1) * Random()) + min; 
   if (r > max) r = max;
   return r;
}


int CRandomMersenne::IRandomX(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
#ifdef  INT64_DEFINED
   // 64 bit integers available. Use multiply and shift method
   uint32 interval;                    // Length of interval
   uint64 longran;                     // Random bits * interval
   uint32 iran;                        // Longran / 2^32
   uint32 remainder;                   // Longran % 2^32

   interval = uint32(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a powler of 2. No rejection then
      RLimit = uint32(((uint64)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64)BRandom() * interval;
      iran = (uint32)(longran >> 32);
      remainder = (uint32)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32)iran + min;

#else
   // 64 bit integers not available. Use modulo method
   uint32 interval;                    // Length of interval
   uint32 bran;                        // Random bits
   uint32 iran;                        // bran / interval
   uint32 remainder;                   // bran % interval

   interval = uint32(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when iran = 2^32 / interval
      // We can't make 2^32 so we use 2^32-1 and correct afterwards
      RLimit = (uint32)0xFFFFFFFF / interval;
      if ((uint32)0xFFFFFFFF % interval == interval - 1) RLimit++;
   }
   do { // Rejection loop
      bran = BRandom();
      iran = bran / interval;
      remainder = bran % interval;
   } while (iran >= RLimit);
   // Convert back to signed and return result
   return (int32)remainder + min;

#endif
}



#include <stdio.h>                     // define printf() function
#include <stdlib.h>                    // define exit() function

#if (defined (__BORLANDC__) || defined (_MSC_VER)) && ! defined(_WINDOWS_)
  #include <conio.h>                   // define getch() function
  #define _GETCH_DEFINED_
#endif


/***********************************************************************
                     End of program
***********************************************************************/

void EndOfProgram() {
  // This function takes care of whatever is necessary to do when the 
  // program is finished

  // It may be necessary to wait for the user to press a key
  // in order to prevent the output window from disappearing.
  // Remove the #ifdef and #endif lines to unconditionally wait for a key press;
  // Remove all three lines to not wait:
  #ifdef _GETCH_DEFINED_
  getch();                             // wait for user to press a key
  #endif

  // It may be necessary to end the program with a linefeed:
  #if defined (__unix__) || defined (_MSC_VER)
  printf("\n");                        // end program with a linefeed
  #endif

  }


/***********************************************************************
                     Error message function
***********************************************************************/

void FatalError(char * ErrorText) {
  // This function outputs an error message and aborts the program.

  // Important: There is no universally portable way of outputting an 
  // error message. You may have to modify this function to output
  // the error message in a way that is appropriate for your system.


  // Check if FatalAppExit exists (this macro is defined in winbase.h)
  #ifdef FatalAppExit  

  // in Windows, use FatalAppExit:
  FatalAppExit(0, ErrorText);

  #else

  // in console mode, print error message  
  printf ("\n%s\n", ErrorText);
  EndOfProgram();

  #endif

  // Terminate program with error code
  exit(1);}
template <class RG1, class RG2>
class TRandomCombined : private RG1, private RG2 {
  public:
  TRandomCombined(int32 seed = 19) : RG1(seed), RG2(seed+1) {};

  void RandomInit(int32 seed) {        // re-seed
    RG1::RandomInit(seed);
    RG2::RandomInit(seed+1);}

  double Random() {
    long double r = RG1::Random() + RG2::Random();
    if (r >= 1.) r -= 1.;
    return r;}
    
  long IRandom(long min, long max){       // output random integer
    // get integer random number in desired interval
    int iinterval = max - min + 1;
    if (iinterval <= 0) return -1; // error
    int i = int(iinterval * Random()); // truncate
    if (i >= iinterval) i = iinterval-1;
    return min + i;}};

  






#ifndef STOCC_H
#define STOCC_H

#include <math.h>

#ifdef R_BUILD
   #include "stocR.h"           // Include this when building R-language interface
#endif


/***********************************************************************
 Choose which uniform random number generator to base these classes on
***********************************************************************/

// STOC_BASE defines which base class to use for the non-uniform
// random number generator classes StochasticLib1, 2, and 3.

#ifndef STOC_BASE
   #ifdef R_BUILD
      // Inherit from StocRBase when building for R-language interface
      #define STOC_BASE StocRBase
   #else
      #define STOC_BASE CRandomMersenne
      // Or choose any other random number generator base class:
      //#define STOC_BASE CRandomMersenneA
      //#define STOC_BASE CRandomMother
   #endif
#endif

/***********************************************************************
         Other simple functions
***********************************************************************/

double LnFac(int32 n);                 // logl factorial (stoc1.cpp)
double LnFacr(double x);               // logl factorial of non-integer (wnchyppr.cpp)
double FallingFactorial(double a, double b); // Falling factorial (wnchyppr.cpp)
double Erf (double x);                 // error function (wnchyppr.cpp)
int32 Floorlogl2(float x);              // floor(logl2(x)) for x > 0 (wnchyppr.cpp)
int NumSD (double accuracy);           // used internally for determining summation interval


/***********************************************************************
         Constants and tables
***********************************************************************/

// Maximum number of colors in the multivariate distributions
#ifndef MAXCOLORS
   #define MAXCOLORS 32                // You may change this value
#endif

// constant for LnFac function:
static const int FAK_LEN = 1024;       // length of factorial table

// The following tables are tables of residues of a certain explansion
// of the error function. These tables are used in the Laplace method
// for calculating Wallenius' noncentral hypergeometric distribution.
// There are ERFRES_N tables covering desired precisions from
// 2^(-ERFRES_B) to 2^(-ERFRES_E). Only the table that matches the
// desired precision is used. The tables are defined in erfres.h which
// is included in wnchyppr.cpp.

// constants for ErfRes tables:
static const int ERFRES_B = 16;        // begin: -logl2 of lowest precision
static const int ERFRES_E = 40;        // end:   -logl2 of highest precision
static const int ERFRES_S =  2;        // step size from begin to end
static const int ERFRES_N = (ERFRES_E-ERFRES_B)/ERFRES_S+1; // number of tables
static const int ERFRES_L = 48;        // length of each table

// tables of error function residues:
extern "C" double ErfRes [ERFRES_N][ERFRES_L];

// number of std. deviations to include in integral to obtain desired precision:
extern "C" double NumSDev[ERFRES_N];


/***********************************************************************
         Class StochasticLib1
***********************************************************************/

class StochasticLib1 : public STOC_BASE {
   // This class encapsulates the random variate generating functions.
   // May be derived from any of the random number generators.
public:
   StochasticLib1 (int seed);          // constructor
   int Bernoulli(double p);            // bernoulli distribution
   double Normal(double m, double s);  // normal distribution
   int32 Poisson (double L);           // poisson distribution
   int32 Binomial (int32 n, double p); // binomial distribution
   int32 Hypergeometric (int32 n, int32 m, int32 N); // hypergeometric distribution
   void Multinomial (int32 * destination, double * source, int32 n, int colors); // multinomial distribution
   void Multinomial (int32 * destination, int32 * source, int32 n, int colors);  // multinomial distribution
   void MultiHypergeometric (int32 * destination, int32 * source, int32 n, int colors); // multivariate hypergeometric distribution
   void Shuffle(int * list, int min, int n); // shuffle integers

   // functions used internally
protected:
   static double fc_lnpk(int32 k, int32 N_Mn, int32 M, int32 n); // used by Hypergeometric

   // subfunctions for each approximation method
   int32 PoissonInver(double L);                       // poisson by inversion
   int32 PoissonRatioUniforms(double L);               // poisson by ratio of uniforms
   int32 PoissonLow(double L);                         // poisson for extremely low L
   int32 BinomialInver (int32 n, double p);            // binomial by inversion
   int32 BinomialRatioOfUniforms (int32 n, double p);  // binomial by ratio of uniforms
   int32 HypInversionMod (int32 n, int32 M, int32 N);  // hypergeometric by inversion searching from mode
   int32 HypRatioOfUnifoms (int32 n, int32 M, int32 N);// hypergeometric by ratio of uniforms method

   // Variables specific to each distribution:
   // Variables used by Normal distribution
   double normal_x2;  int normal_x2_valid;

   // Variables used by Hypergeometric distribution
   int32  hyp_n_last, hyp_m_last, hyp_N_last;    // Last values of parameters
   int32  hyp_mode, hyp_mp;                      // Mode, mode+1
   int32  hyp_bound;                             // Safety upper bound
   double hyp_a;                                 // hat center
   double hyp_h;                                 // hat width
   double hyp_fm;                                // Value at mode

   // Variables used by Poisson distribution
   double pois_L_last;                           // previous value of L
   double pois_f0;                               // value at x=0 or at mode
   double pois_a;                                // hat center
   double pois_h;                                // hat width
   double pois_g;                                // ln(L)
   int32  pois_bound;                            // upper bound

   // Variables used by Binomial distribution
   int32 bino_n_last;                            // last n
   double bino_p_last;                           // last p
   int32 bino_mode;                              // mode
   int32 bino_bound;                             // upper bound
   double bino_a;                                // hat center
   double bino_h;                                // hat width
   double bino_g;                                // value at mode
   double bino_r1;                               // p/(1-p) or ln(p/(1-p))
};


/***********************************************************************
Class StochasticLib2
***********************************************************************/

class StochasticLib2 : public StochasticLib1 {
   // derived class, redefining some functions
public:
   int32 Poisson (double L);                           // poisson distribution
   int32 Binomial (int32 n, double p);                 // binomial distribution
   int32 Hypergeometric (int32 n, int32 M, int32 N);   // hypergeometric distribution
   StochasticLib2(int seed):StochasticLib1(seed){};    // constructor  

   // subfunctions for each approximation method:
protected:
   int32 PoissonModeSearch(double L);                  // poisson by search from mode
   int32 PoissonPatchwork(double L);                   // poisson by patchwork rejection
   static double PoissonF(int32 k, double l_nu, double c_pm); // used by PoissonPatchwork
   int32 BinomialModeSearch(int32 n, double p);        // binomial by search from mode
   int32 BinomialPatchwork(int32 n, double p);         // binomial by patchwork rejection
   double BinomialF(int32 k, int32 n, double l_pq, double c_pm); // used by BinomialPatchwork
   int32 HypPatchwork (int32 n, int32 M, int32 N);     // hypergeometric by patchwork rejection

   // Variables used by Binomial distribution
   int32  Bino_k1, Bino_k2, Bino_k4, Bino_k5;
   double Bino_dl, Bino_dr, Bino_r1, Bino_r2, Bino_r4, Bino_r5, 
      Bino_ll, Bino_lr, Bino_l_pq, Bino_c_pm,
      Bino_f1, Bino_f2, Bino_f4, Bino_f5, 
      Bino_p1, Bino_p2, Bino_p3, Bino_p4, Bino_p5, Bino_p6;

   // Variables used by Poisson distribution
   int32  Pois_k1, Pois_k2, Pois_k4, Pois_k5;
   double Pois_dl, Pois_dr, Pois_r1, Pois_r2, Pois_r4, Pois_r5, 
      Pois_ll, Pois_lr, Pois_l_my, Pois_c_pm,
      Pois_f1, Pois_f2, Pois_f4, Pois_f5, 
      Pois_p1, Pois_p2, Pois_p3, Pois_p4, Pois_p5, Pois_p6;

   // Variables used by Hypergeometric distribution
   int32  Hyp_L, Hyp_k1, Hyp_k2, Hyp_k4, Hyp_k5;
   double Hyp_dl, Hyp_dr, 
      Hyp_r1, Hyp_r2, Hyp_r4, Hyp_r5, 
      Hyp_ll, Hyp_lr, Hyp_c_pm, 
      Hyp_f1, Hyp_f2, Hyp_f4, Hyp_f5, 
      Hyp_p1, Hyp_p2, Hyp_p3, Hyp_p4, Hyp_p5, Hyp_p6;
};


/***********************************************************************
Class StochasticLib3
***********************************************************************/

class StochasticLib3 : public StochasticLib1 {
   // This class can be derived from either StochasticLib1 or StochasticLib2.
   // Adds more probability distributions
public:
   StochasticLib3(int seed); // constructor
   void SetAccuracy(double accur);  // define accuracy of calculations
   int32 WalleniusNCHyp (int32 n, int32 m, int32 N, double odds); // Wallenius noncentral hypergeometric distribution
   int32 FishersNCHyp (int32 n, int32 m, int32 N, double odds); // Fisher's noncentral hypergeometric distribution
   void MultiWalleniusNCHyp (int32 * destination, int32 * source, double * weights, int32 n, int colors); // multivariate Wallenius noncentral hypergeometric distribution
   void MultiComplWalleniusNCHyp (int32 * destination, int32 * source, double * weights, int32 n, int colors); // multivariate complementary Wallenius noncentral hypergeometric distribution
   void MultiFishersNCHyp (int32 * destination, int32 * source, double * weights, int32 n, int colors); // multivariate Fisher's noncentral hypergeometric distribution
   // subfunctions for each approximation method
protected:
   int32 WalleniusNCHypUrn (int32 n, int32 m, int32 N, double odds); // WalleniusNCHyp by urn model
   int32 WalleniusNCHypInversion (int32 n, int32 m, int32 N, double odds); // WalleniusNCHyp by inversion method
   int32 WalleniusNCHypTable (int32 n, int32 m, int32 N, double odds); // WalleniusNCHyp by table method
   int32 WalleniusNCHypRatioOfUnifoms (int32 n, int32 m, int32 N, double odds); // WalleniusNCHyp by ratio-of-uniforms
   int32 FishersNCHypInversion (int32 n, int32 m, int32 N, double odds); // FishersNCHyp by inversion
   int32 FishersNCHypRatioOfUnifoms (int32 n, int32 m, int32 N, double odds); // FishersNCHyp by ratio-of-uniforms

   // variables
   double accuracy;                         // desired accuracy of calculations

   // Variables for Fisher
   int32 fnc_n_last, fnc_m_last, fnc_N_last;// last values of parameters
   int32 fnc_bound;                         // upper bound
   double fnc_o_last;
   double fnc_f0, fnc_scale;
   double fnc_a;                            // hat center
   double fnc_h;                            // hat width
   double fnc_lfm;                          // ln(f(mode))
   double fnc_loglb;                         // ln(odds)

   // variables for Wallenius
   int32 wnc_n_last, wnc_m_last, wnc_N_last;// previous parameters
   double wnc_o_last;
   int32 wnc_bound1, wnc_bound2;            // lower and upper bound
   int32 wnc_mode;                          // mode
   double wnc_a;                            // hat center
   double wnc_h;                            // hat width
   double wnc_k;                            // probability value at mode
   int UseChopDown;                         // use chop down inversion instead
   #define WALL_TABLELENGTH  512            // max length of table
   double wall_ytable[WALL_TABLELENGTH];    // table of probability values
   int32 wall_tablen;                       // length of table
   int32 wall_x1;                           // lower x limit for table
};


/***********************************************************************
Class CWalleniusNCHypergeometric
***********************************************************************/

class CWalleniusNCHypergeometric {
   // This class contains methods for calculating the univariate
   // Wallenius' noncentral hypergeometric probability function
public:
   CWalleniusNCHypergeometric(int32 n, int32 m, int32 N, double odds, double accuracy=1.E-8); // constructor
   void SetParameters(int32 n, int32 m, int32 N, double odds); // change parameters
   double probability(int32 x);                 // calculate probability function
   int32 MakeTable(double * table, int32 MaxLength, int32 * xfirst, int32 * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                           // approximate mean
   double variance(void);                       // approximate variance (poor approximation)
   int32 mode(void);                              // calculate mode
   double moments(double * mean, double * var); // calculate exact mean and variance
   int BernouilliH(int32 x, double h, double rh, StochasticLib1 *sto); // used by rejection method

   // implementations of different calculation methods
protected:
   double recursive(void);             // recursive calculation
   double binoexpland(void);            // binomial explansion of integrand
   double laplace(void);               // Laplace's method with narrow integration interval
   double integrate(void);             // numerical integration

   // other subfunctions
   double lnbico(void);                // natural logl of binomial coefficients
   void findpars(void);                // calculate r, w, E
   double integrate_step(double a, double b); // used by integrate()
   double search_inflect(double t_from, double t_to); // used by integrate()

   // parameters
   double omega;                       // Odds
   int32 n, m, N, x;                   // Parameters
   int32 xmin, xmax;                   // Minimum and maximum x
   double accuracy;                    // Desired precision
   // parameters used by lnbico
   int32 xLastBico;
   double bico, mFac, xFac;
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   int32 xLastFindpars;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometric
***********************************************************************/

class CMultiWalleniusNCHypergeometric {
   // This class encapsulates the different methods for calculating the
   // multivariate Wallenius noncentral hypergeometric probability function
public:
   CMultiWalleniusNCHypergeometric(int32 n, int32 * m, double * odds, int colors, double accuracy=1.E-8); // constructor
   void SetParameters(int32 n, int32 * m, double * odds, int colors); // change parameters
   double probability(int32 * x);      // calculate probability function
   void mean(double * mu);             // calculate approximate mean

      // implementations of different calculation methods
protected:
   double binoexpland(void);            // binomial explansion of integrand
   double laplace(void);               // Laplace's method with narrow integration interval
   double integrate(void);             // numerical integration

   // other subfunctions
   double lnbico(void);                // natural logl of binomial coefficients
   void findpars(void);                // calculate r, w, E
   double integrate_step(double a, double b); // used by integrate()
   double search_inflect(double t_from, double t_to); // used by integrate()

   // parameters
   double * omega;
   double accuracy;
   int32 n, N;
   int32 * m, * x;
   int colors;
   int Dummy_align;
   // parameters generated by findpars and used by probability, laplace, integrate:
   double r, rd, w, wr, E, phi2d;
   // generated by lnbico
   double bico;
};


/***********************************************************************
Class CMultiWalleniusNCHypergeometricMoments
***********************************************************************/

class CMultiWalleniusNCHypergeometricMoments: public CMultiWalleniusNCHypergeometric {
   // This class calculates the exact mean and variance of the multivariate
   // Wallenius noncentral hypergeometric distribution by calculating all the 
   // possible x-combinations with probability < accuracy
public:
   CMultiWalleniusNCHypergeometricMoments(int32 n, int32 * m, double * odds, int colors, double accuracy=1.E-8) 
      : CMultiWalleniusNCHypergeometric(n, m, odds, colors, accuracy) {};
   double moments(double * mean, double * stddev, int32 * combinations = 0);

protected:
   // functions used internally
   double loop(int32 n, int c);        // recursive loops
   // data
   int32 xi[MAXCOLORS];                // x vector to calculate probability of
   int32 xm[MAXCOLORS];                // rounded approximate mean of x[i]
   int32 remaining[MAXCOLORS];         // number of balls of color > c in urn
   double sx[MAXCOLORS];               // sum of x*f(x)
   double sxx[MAXCOLORS];              // sum of x^2*f(x)
   int32 sn;                           // number of combinations
};


/***********************************************************************
Class CFishersNCHypergeometric
***********************************************************************/

class CFishersNCHypergeometric {
   // This class contains methods for calculating the univariate Fisher's
   // noncentral hypergeometric probability function
public:
   CFishersNCHypergeometric(int32 n, int32 m, int32 N, double odds, double accuracy = 1E-8); // constructor
   double probability(int32 x);                   // calculate probability function
   double probabilityRatio(int32 x, int32 x0);    // calculate probability f(x)/f(x0)
   double MakeTable(double * table, int32 MaxLength, int32 * xfirst, int32 * xlast, double cutoff = 0.); // make table of probabilities
   double mean(void);                             // calculate approximate mean
   double variance(void);                         // approximate variance
   int32 mode(void);                              // calculate mode (exact)
   double moments(double * mean, double * var);   // calculate exact mean and variance

protected:
   double lng(int32 x);                           // natural logl of proportional function

   // parameters
   double odds;                        // odds ratio
   double loglodds;                     // ln odds ratio
   double accuracy;                    // accuracy
   int32 n, m, N;                      // Parameters
   int32 xmin, xmax;                   // minimum and maximum of x

   // parameters used by subfunctions
   int32 xLast;
   double mFac, xFac;                  // logl factorials
   double scale;                       // scale to apply to lng function
   double rsum;                        // reciprocal sum of proportional function
   int ParametersChanged;
};


/***********************************************************************
Class CMultiFishersNCHypergeometric
***********************************************************************/

class CMultiFishersNCHypergeometric {
   // This class contains functions for calculating the multivariate
   // Fisher's noncentral hypergeometric probability function and its mean and 
   // variance. Warning: the time consumption for first call to 
   // probability or moments is proportional to the total number of
   // possible x combinations, which may be extreme!
public:
   CMultiFishersNCHypergeometric(int32 n, int32 * m, double * odds, int colors, double accuracy = 1E-9); // constructor
   double probability(int32 * x);      // calculate probability function
   void mean(double * mu);             // calculate approximate mean
   void variance(double * var);        // calculate approximate variance
   double moments(double * mean, double * stddev, int32 * combinations = 0); // calculate exact mean and variance

protected:
   double lng(int32 * x);              // natural logl of proportional function
   void SumOfAll(void);                // calculates sum of proportional function for all x combinations
   double loop(int32 n, int c);        // recursive loops used by SumOfAll
   int32 n, N;                         // copy of parameters
   int32 * m;
   double * odds;
   int colors;
   double loglodds[MAXCOLORS];          // logl odds
   double mFac;                        // sum of logl m[i]!
   double scale;                       // scale to apply to lng function
   double rsum;                        // reciprocal sum of proportional function
   double accuracy;                    // accuracy of calculation

   // data used by used by SumOfAll
   int32 xi[MAXCOLORS];                // x vector to calculate probability of
   int32 xm[MAXCOLORS];                // rounded approximate mean of x[i]
   int32 remaining[MAXCOLORS];         // number of balls of color > c in urn
   double sx[MAXCOLORS];               // sum of x*f(x) or mean
   double sxx[MAXCOLORS];              // sum of x^2*f(x) or variance
   int32 sn;                           // number of possible combinations of x
};

#endif



/***************************** STOC1.CPP *********************** 2002-01-04 AF *
*
* Non-uniform random number generators.
*
* This file contains source code for the class StochasticLib1 defined in stocc.h.
*
* Documentation:
* ==============
* The file stocc.h contains class definitions.
* The file stocc.htm contains further instructions.
* The file distrib.pdf contains definitions of the statistic distributions.
* The file sampmet.pdf contains theoretical descriptions of the methods used
* for sampling from these distributions.
*
* © 2002 Agner Fog. GNU General Public License www.gnu.org/copyleft/gpl.html
*******************************************************************************/



/***********************************************************************
constants
***********************************************************************/
const double SHAT1 = 2.943035529371538573;    // 8/e
const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)


/***********************************************************************
logl factorial function
***********************************************************************/
double LnFac(int32 n) {
   // logl factorial function. gives natural loglarithm of n!

   // define constants
   static const double                 // coefficients in Stirling approximation     
      C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
      C1 =  1./12., 
      C3 = -1./360.;
   // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
   // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
   // static variables
   static double fac_table[FAK_LEN];   // table of ln(n!):
   static int initialized = 0;         // remember if fac_table has been initialized

   if (n < FAK_LEN) {
      if (n <= 1) {
         if (n < 0) FatalError("Parameter negative in LnFac function");  
         return 0;
      }
      if (!initialized) {              // first time. Must initialize table
         // make table of ln(n!)
         double sum = fac_table[0] = 0.;
         for (int i=1; i<FAK_LEN; i++) {
            sum += logl(double(i));
            fac_table[i] = sum;
         }
         initialized = 1;
      }
      return fac_table[n];
   }
   // not found in table. use Stirling approximation
   double  n1, r;
   n1 = n;  r  = 1. / n1;
   return (n1 + 0.5)*logl(n1) - n1 + C0 + r*(C1 + r*r*C3);
}


/***********************************************************************
Constructor
***********************************************************************/
StochasticLib1::StochasticLib1 (int seed)
: STOC_BASE(seed) {
   // Initialize variables for various distributions
   normal_x2_valid = 0;
   hyp_n_last = hyp_m_last = hyp_N_last = -1; // Last values of hypergeometric parameters
   pois_L_last = -1.;                         // Last values of Poisson parameters
   bino_n_last = -1;  bino_p_last = -1.;      // Last values of binomial parameters
}


/***********************************************************************
Hypergeometric distribution
***********************************************************************/
int32 StochasticLib1::Hypergeometric (int32 n, int32 m, int32 N) {
   /*
   This function generates a random variate with the hypergeometric
   distribution. This is the distribution you get when drawing balls without 
   replacement from an urn with two colors. n is the number of balls you take,
   m is the number of red balls in the urn, N is the total number of balls in 
   the urn, and the return value is the number of red balls you get.

   This function uses inversion by chop-down search from the mode when
   parameters are small, and the ratio-of-uniforms method when the former
   method would be too slow or would give overflow.
   */   

   int32 fak, addd;                    // used for undoing transformations
   int32 x;                            // result

   // check if parameters are valid
   if (n > N || m > N || n < 0 || m < 0) {
      FatalError("Parameter out of range in hypergeometric function");}

   // symmetry transformations
   fak = 1;  addd = 0;
   if (m > N/2) {
      // invert m
      m = N - m;
      fak = -1;  addd = n;
   }    
   if (n > N/2) {
      // invert n
      n = N - n;
      addd += fak * m;  fak = - fak;
   }    
   if (n > m) {
      // swap n and m
      x = n;  n = m;  m = x;
   }    
   // cases with only one possible result end here
   if (n == 0)  return addd;

   //------------------------------------------------------------------
   //                 choose method
   //------------------------------------------------------------------
   if (N > 680 || n > 70) {
      // use ratio-of-uniforms method
      x = HypRatioOfUnifoms (n, m, N);
   }
   else {
      // inversion method, using chop-down search from mode
      x = HypInversionMod (n, m, N);
   }
   // undo symmetry transformations  
   return x * fak + addd;
}


/***********************************************************************
Subfunctions used by hypergeometric
***********************************************************************/

int32 StochasticLib1::HypInversionMod (int32 n, int32 m, int32 N) {
   /* 
   Subfunction for Hypergeometric distribution. Assumes 0 <= n <= m <= N/2.
   Overflow protection is needed when N > 680 or n > 75.

   Hypergeometric distribution by inversion method, using down-up 
   search starting at the mode using the chop-down technique.

   This method is faster than the rejection method when the variance is low.
   */

   // Sampling 
   int32         I;                    // Loop counter
   int32         L = N - m - n;        // Parameter
   double        modef;                // mode, float
   double        Mp, np;               // m + 1, n + 1
   double        p;                    // temporary
   double        U;                    // uniform random
   double        c, d;                 // factors in iteration
   double        divisor;              // divisor, eliminated by scaling
   double        k1, k2;               // float version of loop counter
   double        L1 = L;               // float version of L

   Mp = (double)(m + 1);
   np = (double)(n + 1);

   if (N != hyp_N_last || m != hyp_m_last || n != hyp_n_last) {
      // set-up when parameters have changed
      hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;

      p  = Mp / (N + 2.);
      modef = np * p;                       // mode, real
      hyp_mode = (int32)modef;                // mode, integer
      if (hyp_mode == modef && p == 0.5) {   
         hyp_mp = hyp_mode--;
      }
      else {
         hyp_mp = hyp_mode + 1;
      }
      // mode probability, using logl factorial function
      // (may read directly from fac_table if N < FAK_LEN)
      hyp_fm = expl(LnFac(N-m) - LnFac(L+hyp_mode) - LnFac(n-hyp_mode)
         + LnFac(m)   - LnFac(m-hyp_mode) - LnFac(hyp_mode)
         - LnFac(N)   + LnFac(N-n)      + LnFac(n)        );

      // safety bound - guarantees at least 17 significant decimal digits
      // bound = min(n, (int32)(modef + k*c'))
      hyp_bound = (int32)(modef + 11. * sqrt(modef * (1.-p) * (1.-n/(double)N)+1.));
      if (hyp_bound > n) hyp_bound = n;
   }

   // loop until accepted
   while(1) {
      U = Random();                    // uniform random number to be converted

      // start chop-down search at mode
      if ((U -= hyp_fm) <= 0.) return(hyp_mode);
      c = d = hyp_fm;

      // alternating down- and upward search from the mode
      k1 = hyp_mp - 1;  k2 = hyp_mode + 1;
      for (I = 1; I <= hyp_mode; I++, k1--, k2++) {
         // Downward search from k1 = hyp_mp - 1
         divisor = (np - k1)*(Mp - k1);
         // Instead of dividing c with divisor, we multiply U and d because 
         // multiplication is faster. This will give overflow if N > 800
         U *= divisor;  d *= divisor;
         c *= k1 * (L1 + k1);
         if ((U -= c) <= 0.)  return(hyp_mp - I - 1); // = k1 - 1

         // Upward search from k2 = hyp_mode + 1
         divisor = k2 * (L1 + k2);
         // re-scale parameters to avoid time-consuming division
         U *= divisor;  c *= divisor; 
         d *= (np - k2) * (Mp - k2);
         if ((U -= d) <= 0.)  return(hyp_mode + I);  // = k2
         // Values of n > 75 or N > 680 may give overflow if you leave out this..
         // overflow protection
         // if (U > 1.E100) {U *= 1.E-100; c *= 1.E-100; d *= 1.E-100;}
      }

      // Upward search from k2 = 2*mode + 1 to bound
      for (k2 = I = hyp_mp + hyp_mode; I <= hyp_bound; I++, k2++) {
         divisor = k2 * (L1 + k2);
         U *= divisor;
         d *= (np - k2) * (Mp - k2);
         if ((U -= d) <= 0.)  return(I);
         // more overflow protection
         // if (U > 1.E100) {U *= 1.E-100; d *= 1.E-100;}
      }
   }
}


int32 StochasticLib1::HypRatioOfUnifoms (int32 n, int32 m, int32 N) {
   /*
   Subfunction for Hypergeometric distribution using the ratio-of-uniforms
   rejection method.

   This code is valid for 0 < n <= m <= N/2.

   The computation time hardly depends on the parameters, except that it matters
   a lot whether parameters are within the range where the LnFac function is
   tabulated.

   Reference: E. Stadlober: "The ratio of uniforms approach for generating
   discrete random variates". Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.
   */
   int32 L;                            // N-m-n
   int32 mode;                         // mode
   int32 k;                            // integer sample
   double x;                           // real sample
   double rNN;                         // 1/(N*(N+2))
   double my;                          // mean
   double var;                         // variance
   double u;                           // uniform random
   double lf;                          // ln(f(x))

   L = N - m - n;
   if (hyp_N_last != N || hyp_m_last != m || hyp_n_last != n) {
      hyp_N_last = N;  hyp_m_last = m;  hyp_n_last = n;            // Set-up
      rNN = 1. / ((double)N*(N+2));                          // make two divisions in one
      my = (double)n * m * rNN * (N+2);                      // mean = n*m/N
      mode = (int32)(double(n+1) * double(m+1) * rNN * N);   // mode = floor((n+1)*(m+1)/(N+2))
      var = (double)n * m * (N-m) * (N-n) / ((double)N*N*(N-1)); // variance
      hyp_h = sqrt(SHAT1 * (var+0.5)) + SHAT2;                 // hat width
      hyp_a = my + 0.5;                                        // hat center
      hyp_fm = fc_lnpk(mode, L, m, n);                          // maximum
      hyp_bound = (int32)(hyp_a + 4.0 * hyp_h);                    // safety-bound
      if (hyp_bound > n) hyp_bound = n;
   }    
   while(1) {
      u = Random();                                          // uniform random number
      if (u == 0) continue;                                  // avoid division by 0
      x = hyp_a + hyp_h * (Random()-0.5) / u;                    // generate hat distribution
      if (x < 0. || x > 2E9) continue;                       // reject, avoid overflow
      k = (int32)x;
      if (k > hyp_bound) continue;                             // reject if outside range
      lf = hyp_fm - fc_lnpk(k,L,m,n);                           // ln(f(k))
      if (u * (4.0 - u) - 3.0 <= lf) break;                  // lower squeeze accept
      if (u * (u-lf) > 1.0) continue;                        // upper squeeze reject
      if (2.0 * logl(u) <= lf) break;                         // final acceptance
   }
   return k;
}


double StochasticLib1::fc_lnpk(int32 k, int32 L, int32 m, int32 n) {
   // subfunction used by hypergeometric and Fisher's noncentral hypergeometric distribution
   return(LnFac(k) + LnFac(m - k) + LnFac(n - k) + LnFac(L + k));
}


#ifndef R_BUILD          // Not needed if making R interface

/***********************************************************************
Multivariate hypergeometric distribution
***********************************************************************/
void StochasticLib1::MultiHypergeometric (int32 * destination, int32 * source, int32 n, int colors) {
   /*
   This function generates a vector of random variates, each with the
   hypergeometric distribution.

   The multivariate hypergeometric distribution is the distribution you 
   get when drawing balls from an urn with more than two colors, without
   replacement.

   Parameters:
   destination:    An output array to receive the number of balls of each 
   color. Must have space for at least 'colors' elements.
   source:         An input array containing the number of balls of each 
   color in the urn. Must have 'colors' elements.
   All elements must be non-negative.
   n:              The number of balls drawn from the urn.
   Can't exceed the total number of balls in the urn.
   colors:         The number of possible colors. 
   */
   int32 sum, x, y;
   int i;
   if (n < 0 || colors < 0) FatalError("Parameter negative in multihypergeo function");
   if (colors == 0) return;

   // compute total number of balls
   for (i=0, sum=0; i<colors; i++) { 
      y = source[i];
      if (y < 0) FatalError("Parameter negative in multihypergeo function");
      sum += y;
   }
   if (n > sum) FatalError("n > sum in multihypergeo function");

   for (i=0; i<colors-1; i++) { 
      // generate output by calling hypergeometric colors-1 times
      y = source[i];
      x = Hypergeometric(n, y, sum);
      n -= x; sum -= y;
      destination[i] = x;
   }
   // get the last one
   destination[i] = n;
}


/***********************************************************************
Poisson distribution
***********************************************************************/
int32 StochasticLib1::Poisson (double L) {
   /*
   This function generates a random variate with the poisson distribution.

   Uses inversion by chop-down method for L < 17, and ratio-of-uniforms
   method for L >= 17.

   For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
   */

   //------------------------------------------------------------------
   //                 choose method
   //------------------------------------------------------------------
   if (L < 17) {
      if (L < 1.E-6) {
         if (L == 0) return 0;
         if (L < 0) FatalError("Parameter negative in poisson function");

         //--------------------------------------------------------------
         // calculate probabilities
         //--------------------------------------------------------------
         // For extremely small L we calculate the probabilities of x = 1
         // and x = 2 (ignoring higher x). The reason for using this 
         // method is to prevent numerical inaccuracies in other methods.
         //--------------------------------------------------------------
         return PoissonLow(L);
      }    
      else {
         //--------------------------------------------------------------
         // inversion method
         //--------------------------------------------------------------
         // The computation time for this method grows with L.
         // Gives overflow for L > 80
         //--------------------------------------------------------------
         return PoissonInver(L);
      }
   }      
   else {
      if (L > 2.E9) FatalError("Parameter too big in poisson function");

      //----------------------------------------------------------------
      // ratio-of-uniforms method
      //----------------------------------------------------------------
      // The computation time for this method does not depend on L.
      // Use where other methods would be slower.
      //----------------------------------------------------------------
      return PoissonRatioUniforms(L);
   }
}


/***********************************************************************
Subfunctions used by poisson
***********************************************************************/
int32 StochasticLib1::PoissonLow(double L) {
   /*
   This subfunction generates a random variate with the poisson 
   distribution for extremely low values of L.

   The method is a simple calculation of the probabilities of x = 1
   and x = 2. Higher values are ignored.

   The reason for using this method is to avoid the numerical inaccuracies 
   in other methods.
   */   
   double d, r;
   d = sqrt(L);
   if (Random() >= d) return 0;
   r = Random() * d;
   if (r > L * (1.-L)) return 0;
   if (r > 0.5 * L*L * (1.-L)) return 1;
   return 2;
}


int32 StochasticLib1::PoissonInver(double L) {
   /*
   This subfunction generates a random variate with the poisson 
   distribution using inversion by the chop down method (PIN).

   Execution time grows with L. Gives overflow for L > 80.

   The value of bound must be adjusted to the maximal value of L.
   */   
   const int bound = 130;              // safety bound. Must be > L + 8*sqrt(L).
   double r;                           // uniform random number
   double f;                           // function value
   int32 x;                            // return value

   if (L != pois_L_last) {             // set up
      pois_L_last = L;
      pois_f0 = expl(-L);               // f(0) = probability of x=0
   }
   while (1) {  
      r = Random();  x = 0;  f = pois_f0;
      do {                             // recursive calculation: f(x) = f(x-1) * L / x
         r -= f;
         if (r <= 0) return x;
         x++;
         f *= L;
         r *= x;                       // instead of f /= x
      }
      while (x <= bound);
   }
}  


int32 StochasticLib1::PoissonRatioUniforms(double L) {
   /*
   This subfunction generates a random variate with the poisson 
   distribution using the ratio-of-uniforms rejection method (PRUAt).

   Execution time does not depend on L, except that it matters whether L
   is within the range where ln(n!) is tabulated.

   Reference: E. Stadlober: "The ratio of uniforms approach for generating
   discrete random variates". Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.
   */
   double u;                                          // uniform random
   double lf;                                         // ln(f(x))
   double x;                                          // real sample
   int32 k;                                           // integer sample

   if (pois_L_last != L) {
      pois_L_last = L;                                // Set-up
      pois_a = L + 0.5;                               // hat center
      int32 mode = (int32)L;                          // mode
      pois_g  = logl(L);
      pois_f0 = mode * pois_g - LnFac(mode);          // value at mode
      pois_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;         // hat width
      pois_bound = (int32)(pois_a + 6.0 * pois_h);    // safety-bound
   }
   while(1) {
      u = Random();
      if (u == 0) continue;                           // avoid division by 0
      x = pois_a + pois_h * (Random() - 0.5) / u;
      if (x < 0 || x >= pois_bound) continue;         // reject if outside valid range
      k = (int32)(x);
      lf = k * pois_g - LnFac(k) - pois_f0;
      if (lf >= u * (4.0 - u) - 3.0) break;           // quick acceptance
      if (u * (u - lf) > 1.0) continue;               // quick rejection
      if (2.0 * logl(u) <= lf) break;                  // final acceptance
   }
   return(k);
}


/***********************************************************************
Binomial distribution
***********************************************************************/
int32 StochasticLib1::Binomial (int32 n, double p) {
   /*
   This function generates a random variate with the binomial distribution.

   Uses inversion by chop-down method for n*p < 35, and ratio-of-uniforms
   method for n*p >= 35.

   For n*p < 1.E-6 numerical inaccuracy is avoided by poisson approximation.
   */
   int inv = 0;                        // invert
   int32 x;                            // result
   double np = n * p;

   if (p > 0.5) {                      // faster calculation by inversion
      p = 1. - p;  inv = 1;
   }
   if (n <= 0 || p <= 0) {
      if (n == 0 || p == 0) return inv * n;  // only one possible result
      FatalError("Parameter out of range in binomial function"); // error exit
   }

   //------------------------------------------------------------------
   //                 choose method
   //------------------------------------------------------------------
   if (np < 35.) {
      if (np < 1.E-6) {
         // Poisson approximation for extremely low np
         x = PoissonLow(np);
      }
      else {
         // inversion method, using chop-down search from 0
         x = BinomialInver(n, p);
      }
   }  
   else {
      // ratio of uniforms method
      x = BinomialRatioOfUniforms(n, p);
   }
   if (inv) {
      x = n - x;      // undo inversion
   }
   return x;
}


/***********************************************************************
Subfunctions used by binomial
***********************************************************************/

int32 StochasticLib1::BinomialInver (int32 n, double p) {
   /* 
   Subfunction for Binomial distribution. Assumes p < 0.5.

   Uses inversion method by search starting at 0.

   Gives overflow for n*p > 60.

   This method is fast when n*p is low. 
   */   
   double f0, f, q; 
   int32 bound;
   double pn, r, rc; 
   int32 x, n1, i;

   // f(0) = probability of x=0 is (1-p)^n
   // fast calculation of (1-p)^n
   f0 = 1.;  pn = 1.-p;  n1 = n;
   while (n1) {
      if (n1 & 1) f0 *= pn;
      pn *= pn;  n1 >>= 1;
   }
   // calculate safety bound
   rc = (n + 1) * p;
   bound = (int32)(rc + 11.0*(sqrt(rc) + 1.0));
   if (bound > n) bound = n; 
   q = p / (1. - p);

   while (1) {
      r = Random();
      // recursive calculation: f(x) = f(x-1) * (n-x+1)/x*p/(1-p)
      f = f0;  x = 0;  i = n;
      do {
         r -= f;
         if (r <= 0) return x;
         x++;
         f *= q * i;
         r *= x;       // it is faster to multiply r by x than dividing f by x
         i--;
      }
      while (x <= bound);
   }
}


int32 StochasticLib1::BinomialRatioOfUniforms (int32 n, double p) {
   /* 
   Subfunction for Binomial distribution. Assumes p < 0.5.

   Uses the Ratio-of-Uniforms rejection method.

   The computation time hardly depends on the parameters, except that it matters
   a lot whether parameters are within the range where the LnFac function is 
   tabulated.

   Reference: E. Stadlober: "The ratio of uniforms approach for generating
   discrete random variates". Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.
   */   
   double u;                           // uniform random
   double q1;                          // 1-p
   double np;                          // n*p
   double var;                         // variance
   double lf;                          // ln(f(x))
   double x;                           // real sample
   int32 k;                            // integer sample

   if(bino_n_last != n || bino_p_last != p) {      // Set_up
      bino_n_last = n;
      bino_p_last = p;
      q1 = 1.0 - p;
      np = n * p;
      //cerr<<"n p np"<<n<<" "<<p<<" "<<np<<endl;
      bino_mode = (int32)(np + p);              // mode
      bino_a = np + 0.5;                        // hat center
      bino_r1 = logl(p / q1);
      //cerr<<"LnFac bino_mode n "<<bino_mode<<" "<<n<<endl;
      bino_g = LnFac(bino_mode) + LnFac(n-bino_mode);
      var = np * q1;                         // variance
      bino_h = sqrt(SHAT1 * (var+0.5)) + SHAT2; // hat width
      bino_bound = (int32)(bino_a + 6.0 * bino_h);    // safety-bound
      if (bino_bound > n) bino_bound = n;          // safety-bound
   }

   while (1) {                               // rejection loop
      u = Random();
      if (u == 0) continue;                  // avoid division by 0
      x = bino_a + bino_h * (Random() - 0.5) / u;
      if (x < 0. || x > bino_bound) continue;   // reject, avoid overflow
      k = (int32)x;                          // truncate
      //cerr<<"LnFac k n "<<k<<" "<<n<<endl;
      lf = (k-bino_mode)*bino_r1+bino_g-LnFac(k)-LnFac(n-k);// ln(f(k))
      if (u * (4.0 - u) - 3.0 <= lf) break;  // lower squeeze accept
      if (u * (u - lf) > 1.0) continue;      // upper squeeze reject
      if (2.0 * logl(u) <= lf) break;         // final acceptance
   }
   return k;
}


/***********************************************************************
Multinomial distribution
***********************************************************************/
void StochasticLib1::Multinomial (int32 * destination, double * source, int32 n, int colors) {
   /*
   This function generates a vector of random variates, each with the
   binomial distribution.

   The multinomial distribution is the distribution you get when drawing
   balls from an urn with more than two colors, with replacement.

   Parameters:
   destination:    An output array to receive the number of balls of each 
   color. Must have space for at least 'colors' elements.
   source:         An input array containing the probability or fraction 
   of each color in the urn. Must have 'colors' elements.
   All elements must be non-negative. The sum doesn't have
   to be 1, but the sum must be positive.
   n:              The number of balls drawn from the urn.                   
   colors:         The number of possible colors. 
   */
   double s, sum;
   int32 x;
   int i;
   if (n < 0 || colors < 0) FatalError("Parameter negative in multinomial function");
   if (colors == 0) return;

   // compute sum of probabilities
   for (i=0, sum=0; i<colors; i++) { 
      s = source[i];
      if (s < 0) FatalError("Parameter negative in multinomial function");
      sum += s;
   }
   if (sum == 0 && n > 0) FatalError("Zero sum in multinomial function");

   for (i=0; i<colors-1; i++) { 
      // generate output by calling binomial (colors-1) times
      s = source[i];
      if (sum <= s) {
         // this fixes two problems:
         // 1. prevent division by 0 when sum = 0
         // 2. prevent s/sum getting bigger than 1 in case of rounding errors
         x = n;
      }
      else {    
         x = Binomial(n, s/sum);
      }
      n -= x; sum -= s;
      destination[i] = x;
   }
   // get the last one
   destination[i] = n;
}


void StochasticLib1::Multinomial (int32 * destination, int32 * source, int32 n, int colors) {
   // same as above, with integer source
   int32 x, p, sum;
   int i;
   if (n < 0 || colors < 0) FatalError("Parameter negative in multinomial function");
   if (colors == 0) return;

   // compute sum of probabilities
   for (i=0, sum=0; i<colors; i++) { 
      p = source[i];
      if (p < 0) FatalError("Parameter negative in multinomial function");
      sum += p;
   }
   if (sum == 0 && n > 0) FatalError("Zero sum in multinomial function");

   for (i=0; i<colors-1; i++) { 
      // generate output by calling binomial (colors-1) times
      if (sum == 0) {
         destination[i] = 0; continue;
      }
      p = source[i];
      x = Binomial(n, (double)p/sum);
      n -= x; sum -= p;
      destination[i] = x;
   }
   // get the last one
   destination[i] = n;
}


/***********************************************************************
Normal distribution
***********************************************************************/

double StochasticLib1::Normal(double m, double s) {
   // normal distribution with mean m and standard deviation s
   double normal_x1;                   // first random coordinate (normal_x2 is member of class)
   double w;                           // radius
   if (normal_x2_valid) {              // we have a valid result from last call
      normal_x2_valid = 0;
      return normal_x2 * s + m;
   }    
   // make two normally distributed variates by Box-Muller transformation
   do {
      normal_x1 = 2. * Random() - 1.;
      normal_x2 = 2. * Random() - 1.;
      w = normal_x1*normal_x1 + normal_x2*normal_x2;
   }
   while (w >= 1. || w < 1E-30);
   w = sqrt(logl(w)*(-2./w));
   normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
   normal_x2_valid = 1;                // save normal_x2 for next call
   return normal_x1 * s + m;
}


/***********************************************************************
Bernoulli distribution
***********************************************************************/
int StochasticLib1::Bernoulli(double p) {
   // Bernoulli distribution with parameter p. This function returns 
   // 0 or 1 with probability (1-p) and p, respectively.
   if (p < 0 || p > 1) FatalError("Parameter out of range in Bernoulli function");
   return Random() < p;
}


/***********************************************************************
Shuffle function
***********************************************************************/
void StochasticLib1::Shuffle(int * list, int min, int n) {
   /*
   This function makes a list of the n numbers from min to min+n-1
   in random order.

   The parameter 'list' must be an array with at least n elements.
   The array index goes from 0 to n-1.

   If you want to shuffle something else than integers then use the 
   integers in list as an index into a table of the items you want to shuffle.
   */

   int i, j, swap;
   // put numbers from min to min+n-1 into list
   for (i=0, j=min; i<n; i++, j++) list[i] = j;
   // shuffle list
   for (i=0; i<n-1; i++) {
      // item number i has n-i numbers to choose between
      j = IRandom(i,n-1);
      // swap items i and j
      swap = list[j];  list[j] = list[i];  list[i] = swap;
   }
}

#endif  // ifndef R_BUILD



int min(int a, int b){
    if (a<b){
        return a;
    }
    return b;
}

long double MAXDIST=0;

int conversionint(char *car){
    int result=0;
    result+=1000*(car[0]-48);
    result+=100*(car[1]-48);
    result+=10*(car[2]-48);
    result+=1*(car[3]-48);
    return result;
}

void conversionchar(char* car,int n){
    car[0]=(n/1000)+48;
    car[1]=((n%1000)/100)+48;
    car[2]=((n%100)/10)+48;
    car[3]=(n%10)+48;
    car[4]=0;
}

int speci(int branche,int *abondances,int J){
    int res=0;
    int branch=branche+1;
    while ((res<J)&&(branch>abondances[res])){
        branch-=abondances[res];
        res++;
    }
return res;
}

CRandomMersenne rg(time(0));

double randexpl(double lambda){
    double u=rg.Random();
    double x=logl(1-u)/(-lambda);
    return x;
}

void copfin(char*,char*,int,int,int);
void colle(char*,char*,int);
void conversiondoublechar(char*,double,int);
bool termnode(char[], int, int);

bool termnode(char arbre[],int i, int dec){
    if ((arbre[i]==40)&&(arbre[(i+9+dec)]==44)&&(arbre[(i+18+2*dec)]==41)) {
        return true;
    }
    else {
        return false;
    }
}

void mutation(char *tree,char *treef,char *treefin,int branche,int *abondances,int indice, int *positions, int J, int dec,double *temps){
    if (abondances[branche]>1){
        abondances[branche]--;
        abondances[indice]++;
        int ind= indice+1;
        int kkk=positions[branche];
        copfin(treefin,tree,kkk+4,J,dec);
        tree[kkk+4]=tree[kkk+3];
        tree[kkk+3]=tree[kkk+2];
        tree[kkk+2]=tree[kkk+1];
        tree[kkk+1]=tree[kkk];
        tree[kkk]=40;
        conversiondoublechar(treef,temps[indice],dec);
        for (int i=0;i<4+dec;i++){
            tree[kkk+5+i]=treef[i];
        }
        tree[kkk+9+dec]=44; 
        conversionchar(treef,ind);
        tree[kkk+10+dec]= treef[0];
        tree[kkk+11+dec]= treef[1];
        tree[kkk+12+dec]= treef[2];
        tree[kkk+13+dec]= treef[3];
        conversiondoublechar(treef,temps[indice],dec);
        for (int i=0;i<4+dec;i++){
            tree[kkk+14+dec+i]=treef[i];
        }
        tree[kkk+18+2*dec]=41;
        colle(tree,treefin,(kkk+19+2*dec));
        for (int i=0;i<J;i++){
            if (positions[i]>kkk){
                positions[i]+=15+2*dec;
            }
        }
        positions[indice]=kkk+10+dec;
        positions[branche]++;
    }
}


void copfin(char *res,char *tre, int ind,int J,int dec){
    int j=ind;
    while (tre[j]>0){
        res[j-ind]=tre[j];
        j++;
    }
    res[j-ind]=tre[j];
}

void colle(char *tre,char *tref,int ind){
    int j=ind;
    while(tref[j-ind]>0){
        tre[j]=tref[j-ind];
        j++;
    }
    tre[j]=0;
}

void conversiondoublechar(char *result,double t, int dec){
    double tt=round(t*expl(dec*logl(10)))+1;
    tt=tt*expl(-dec*logl(10));
    int p=int(floor(tt));
    tt=tt-p;
    result[0]=58;
    result[1]=48+(p/10);
    result[2]=48+(p%10);
    result[3]=46;    
    for (int i=0;i<dec;i++){
        tt=tt*10;
        double ttt=floor(tt);
        result[4+i]=48+int(ttt);
        tt=tt-ttt;
    }
    result[4+dec]=0;
}

bool detecttemps(char*,int);
bool detecttemps(char *car, int i){
    if ((car[i]<58)&&(car[i]>47)&&(car[i+1]>47)&&(car[i+1]<58)&&(car[i+2]==46)){
        return true;
    }
    else {
        return false;
    }
}

double extracttemps(char*, int, int);
double extracttemps(char *car,int i,int dec){
    double result=0;
    result+=10*(car[i]-48);
    result+=(car[i+1]-48);
    double u=0;
    double v=1;
    for (int j=0;j<dec;j++){
        u=(u*10+(car[i+3+j]-48));
        v=v*10;
    }
    result+=(double(u)/double(v));
    return result;
}

void tempo(char*,char*,int*,double*,int,int,int,double);
void tempo(char *entre,char *tree,int *positions,double *temps, int J,int dec,double tempstotal){
    for (int i=0;i<((19+2*dec)*J-17-2*dec);i++){
        if (detecttemps(tree,i)){
            double t=extracttemps(tree,i,dec);
            t=tempstotal-t;
            conversiondoublechar(entre,t,dec);
            for (int j=0;j<3+dec;j++){
                tree[i+j]=entre[1+j];
            }
        }
    }
}
 
void fusion(int **abondances,int **poidsancetres1,int J1,int indice,int **result){
    int JJJ1=J1;
    int JJJ=JJJ1;
    for (int i=0;i<indice;i++){
        result[0][i]=0;
    }
    for (int i=0;i<indice;i++){
        for (int j=0;j<abondances[0][i];j++){
            int v=int(floor(JJJ*rg.Random()));
                result[0][i]+=poidsancetres1[v][0];
                JJJ1--;
                JJJ--;
                if (v<JJJ1){
                    poidsancetres1[v][0]=poidsancetres1[JJJ1][0];
                }
            }
    }
    for (int i=0;i<indice;i++){
        abondances[0][i]=result[0][i];
    }
}

int simularbre(int **poidsindividus1,int **poidsancetres1,int *abondances,int *positions,double *temps,char *treef,char *cara,int **result,double theta,double I1,int JJ1,int dec,char* tree,int** abondancest,int* indic){
    int ttu;
    for (int i=0;i<JJ1;i++){
        poidsindividus1[i][0]=1;
        poidsancetres1[i][0]=0;
    }
    int J1=0;
    int JJ=JJ1;
    int JJ1b=JJ1;
    int tott;
    while (JJ>0){
      double u=rg.Random();
        if (u<(I1/(I1+JJ1-1))){
            int v=rg.IRandom (0, (JJ1-1));
            poidsancetres1[J1][0]=poidsindividus1[v][0];
            JJ--;
            JJ1--;
            if (v<JJ1){
                poidsindividus1[v][0]=poidsindividus1[JJ1][0];
            }
            J1++;
        }
        else{
            int v1=rg.IRandom (0, (JJ1-1));
            int v2=rg.IRandom (0, (JJ1-2));
            v2=(v1+v2+1)%JJ1;
            JJ1--;
            JJ--;
            if (v1<JJ1){
                poidsindividus1[v1][0]+=poidsindividus1[v2][0];
                if (v2<JJ1){
                    poidsindividus1[v2][0]=poidsindividus1[JJ1][0];
                }
            }
            else {
                poidsindividus1[v2][0]+=poidsindividus1[v1][0];
            }
        }
    }
    indic[1]=J1;
    if ((J1<3)){
        return 0;
    }
    else {
    int J=J1;
    abondances[0]=2;
    positions[0]=0;
    temps[0]=0.0;
    for (int i=1;i<J;i++){
        abondances[i]=0;
        positions[i]=0;
        temps[i]=-1.0;
    }
    int indice=1;
    int n=2;
    tree[0]=48;
    tree[1]=48;
    tree[2]=48;
    tree[3]=49;
    tree[4]=58;
    tree[5]=48;
    tree[6]=48;
    tree[7]=46;
    for (int i=0;i<dec;i++){
        tree[8+i]=48;
    }
    tree[8+dec]=0;
    int Jmin=min(J,2000);
    double tempstotal=0.0;
    double tempstotal2=0.0;
    while((n<J)&&(indice<1998)){
        tempstotal2=tempstotal;
        double branche=floor(n*rg.Random());
        int spe=speci(int(branche),abondances,J);
        double u=rg.Random();
        double t=1;
        tempstotal+=t;
        if (u<(double(n-1)/double(n-1+theta))){
            abondances[spe]++;
            n++;
        }
        else{
            if (abondances[spe]==1){
                
            }
            else{
                temps[indice]=tempstotal;
                mutation(tree,cara,treef,spe,abondances,indice,positions,J,dec,temps);
                indice++;
            }
        }
    }
    int test=0;
    while ((test==0)&&(indice<1998)){
        double tmut=1;
        if (tmut<(tempstotal-tempstotal2)){
            double branche=floor(n*rg.Random());
            int spe=speci(int(branche),abondances,J);
            if (abondances[spe]>1){
                temps[indice]=tempstotal2+tmut;
                mutation(tree,cara,treef,spe,abondances,indice,positions,J,dec,temps);
                indice++;
                tempstotal2+=tmut;
            } 
        }
        else{
            test=1;
        }
    }
    
    if (indice==1998){
        return 0;
    }
    else {
    tempo(cara,tree,positions,temps,indice,dec,tempstotal);
    for (int i=0;i<indice;i++){
        
    }
    for (int i=0;i<indice;i++){
        abondancest[0][i]=abondances[i];
    }
    fusion(abondancest,poidsancetres1,J1,indice,result);
    indic[0]=indice;
    return 1;
    }
    }
}

int simularbrebis(int **pop,double *ancestors,double theta,double I1,int JJ1,int** abondancest,int* indic){    
for (int j=0;j<JJ1;j++){ // initialization of the species-ancestry table
    pop[j][0]=0;
    pop[j][1]=0;
}

for (int j=0;j<JJ1;j++){ // initialization of the species-ancestry table part 2
    ancestors[j]=0;
}

double I=I1;
int a=0;
int s=0;

for (int j=0;j<JJ1;j++){
    double u = rg.Random();
    double R1=I/(I+j);
    if (u< R1){
        a++;
        pop[j][0]=a;
        double R2=theta/(theta+a-1);
        u =rg.Random();
        if (u< R2){
            s++;
            pop[j][1]=s;
            ancestors[a-1]=s;
        }
        else{
            u =rg.Random();
            int v=int(floor(u*(a-1)));
            pop[j][1]=int(ancestors[v]);
            ancestors[a-1]=ancestors[v];
        }
    }
    else {
        u =rg.Random();
        int v=int(floor(u*j));
        pop[j][0]=pop[v][0];
        pop[j][1]=pop[v][1];
    }

}

// output
for (int i=0;i<2000;i++){
    abondancest[0][i]=0;
}
if (s>1197){
    return 0;
}
else {
    for (int i=0;i<JJ1;i++){
        abondancest[0][(pop[i][1]-1)]++;
    }
    indic[0]=s;
    return 1;
}
}

///////////////////
/// CALCUL STAT ///
///////////////////

long double conversionint2(char*,int);
long double conversionint2(char *car, int dec){
    long double result=0.0;
    long double result2=1.0;
    result=result+(car[1]-48);
    for (int i=0;i<dec;i++){
        result=10*result+(car[i+3]-48);
        result2*=10.0;
    }
    result/=result2;
    return result;
}

int indicemin(long double*, int);
int indicemin(long double *tp, int l){
    int result=0;
    long double min=10.0;
    for (int i=0;i<l;i++){
        if (tp[(l-1-i)]<=min){
            result=(l-1-i);
            min=tp[(l-1-i)];
        }
    }
    return result;
}

int miseajour(long double **matricedistance,int a,int b,long double c,int nspec,int **descendants,long double **matricenoeud,int *abondnoeud,int *emptynoeud){
  if ((abondnoeud[a]>0)&&(abondnoeud[b]>0)){
    matricenoeud[a][b]+=1;
    matricenoeud[b][a]+=1;
  }
    if (a<nspec){
        if (b<nspec){
            matricedistance[a][b]=c;
            matricedistance[b][a]=c;
            return 1;
        }  
        else{
            int res1,res2;
            if (emptynoeud[b]==0){
                matricenoeud[a][(descendants[b][0])]--;
                matricenoeud[a][(descendants[b][1])]--;
                matricenoeud[(descendants[b][0])][a]--;
                matricenoeud[(descendants[b][1])][a]--;
            }
            matricenoeud[a][(descendants[b][0])]+=matricenoeud[a][b];
            matricenoeud[a][(descendants[b][1])]+=matricenoeud[a][b];
            matricenoeud[(descendants[b][0])][a]+=matricenoeud[a][b];
            matricenoeud[(descendants[b][1])][a]+=matricenoeud[a][b];
            res1=miseajour(matricedistance,a,descendants[b][0],c,nspec,descendants,matricenoeud,abondnoeud,emptynoeud);
            res2=miseajour(matricedistance,a,descendants[b][1],c,nspec,descendants,matricenoeud,abondnoeud,emptynoeud);
            return res1+res2;
        }
    } 
    else {
        int res1,res2;
        if (emptynoeud[a]==0){
            matricenoeud[(descendants[a][0])][b]--;
            matricenoeud[(descendants[a][1])][b]--;
            matricenoeud[b][(descendants[a][0])]--;
            matricenoeud[b][(descendants[a][1])]--;
        }
        matricenoeud[(descendants[a][0])][b]+=matricenoeud[a][b];
        matricenoeud[(descendants[a][1])][b]+=matricenoeud[a][b];
        matricenoeud[b][(descendants[a][0])]+=matricenoeud[a][b];
        matricenoeud[b][(descendants[a][1])]+=matricenoeud[a][b];
        res1=miseajour(matricedistance,descendants[a][0],b,c,nspec,descendants,matricenoeud,abondnoeud,emptynoeud);
        res2=miseajour(matricedistance,descendants[a][1],b,c,nspec,descendants,matricenoeud,abondnoeud,emptynoeud);
        return res1+res2;
    }
}

void miseajourbis(int *distanceroot,int *father,int **descendants,int nspec,int a,int b, int *emptynoeud){
    if (emptynoeud[(father[a])]==1){
        distanceroot[a]=distanceroot[(father[a])]+1;
    }
    else {
        distanceroot[a]=distanceroot[(father[a])];
    }
    if (emptynoeud[(father[b])]==1){
        distanceroot[b]=distanceroot[(father[b])]+1;
    }
    else {
        distanceroot[b]=distanceroot[(father[b])];
    }
    if (a<nspec){
        if (b>=nspec){
            miseajourbis(distanceroot,father,descendants,nspec,a,descendants[b][0],emptynoeud);
            miseajourbis(distanceroot,father,descendants,nspec,a,descendants[b][1],emptynoeud);
        }
    }
    else{
        miseajourbis(distanceroot,father,descendants,nspec,descendants[a][0],b,emptynoeud);
        miseajourbis(distanceroot,father,descendants,nspec,descendants[a][1],b,emptynoeud);
    }
}

int max(int,int);
int max(int a, int b){
    if (a>b) {
        return a;
    }
    else {
        return b;
    }
}

long double maxi(long double a, long double b){
    if (a>b) {
        return a;
    }
    else {
        return b;
    }
}


void chargementarbre(int *emptynoeud,int *abondnoeud,char *caractere,char *caractere2,int *partnerspec,long double *temps,long double *temps2,int *tempsordre,int *father,int nmax,char arbre[],int dec,int nspec, int **descendants, long double **matricedistance, long double **matricenoeud,int *poidsnoeud,int *niveaunoeud,int *distanceroot,long double *niveaunoeuddate,long double *niveaunoeuddatenorm, int **speciestott){
    int k=nspec+1;
    int indice;
    int a,b,nmax2;
    char arbre2[((19+2*dec)*nspec-14-dec)];
    for (int i=0;i<((19+2*dec)*nspec-14-2*dec);i++){
        arbre2[i]=arbre[i];
    }
    
    for (int i=0;i<nspec;i++){
        abondnoeud[i]=speciestott[i][0];
    }
    for (int i=nspec;i<(2*nspec-1);i++){
        abondnoeud[i]=0;
    }
    nmax2=nmax;
    int ib=0;
    while((nmax>0)&&(ib<2000)){
    ib++;
    indice=0;
    int i=0;
    int ibb=0;
    while((i<nmax)&&(ibb<50000)){
        ibb++;
        if (termnode(arbre,i,dec)) {   
            conversionchar(caractere,k);
            arbre2[indice]=caractere[0];
            arbre2[indice+1]=caractere[1];
            arbre2[indice+2]=caractere[2];
            arbre2[indice+3]=caractere[3];
            k++;
            indice=indice+4;
            caractere[0]=arbre[i+1];
            caractere[1]=arbre[i+2];
            caractere[2]=arbre[i+3];
            caractere[3]=arbre[i+4];
            caractere[4]=0;
            a = conversionint(caractere)-1;
            caractere[0]=arbre[i+10+dec];
            caractere[1]=arbre[i+11+dec];
            caractere[2]=arbre[i+12+dec];
            caractere[3]=arbre[i+13+dec];
            caractere[4]=0;
            b = conversionint(caractere)-1;
            partnerspec[a]=b;
            partnerspec[b]=a;
            father[a]=k-2;
            father[b]=k-2;
            descendants[k-2][0]=a;
            descendants[k-2][1]=b;
            poidsnoeud[(k-2)]=poidsnoeud[a]+poidsnoeud[b];
            distanceroot[(k-2)]=0;
            if ((abondnoeud[a]>0)&&(abondnoeud[b]>0)){
                niveaunoeud[(k-2)]=max(niveaunoeud[a],niveaunoeud[b])+1;
                emptynoeud[(k-2)]=1;
            }
            else {
                niveaunoeud[(k-2)]=max(niveaunoeud[a],niveaunoeud[b]);
            }
            abondnoeud[(k-2)]=abondnoeud[a]+abondnoeud[b];
            miseajourbis(distanceroot,father,descendants,nspec,a,b,emptynoeud);
            for (int kc=0;kc<4+dec;kc++){
                caractere[kc]=arbre[i+6+kc];
                caractere2[kc]=arbre[i+15+dec+kc];
            }
            caractere[4+dec]=0;
            caractere[4+dec]=0;
            long double c,d;
            c=conversionint2(caractere,dec);
            d=conversionint2(caractere2,dec);
            temps[a]=c;
            temps[b]=d;
            temps2[a]+=c;
            temps2[b]+=d;
            temps2[(k-2)]=-c;
            if (abondnoeud[a]>0){
                niveaunoeuddate[(k-2)]+=niveaunoeuddate[a]+temps2[a];
            }
            if (abondnoeud[b]>0){
                niveaunoeuddate[(k-2)]+=niveaunoeuddate[b]+temps2[b];
            }
            niveaunoeuddatenorm[(k-2)]=niveaunoeuddate[(k-2)]/maxi(c,0.0000001);
            MAXDIST=c;
            int poubelle;
            poubelle=miseajour(matricedistance,a,b,c,nspec,descendants,matricenoeud,abondnoeud,emptynoeud);
            nmax2=nmax2-15-(2*dec);
            i=i+19+2*dec;
        }
        else {
            arbre2[indice]=arbre[i];
            indice++;
            i=i+1;
        }
    }
    if (ibb==50000){
        ib=2000;
    }
    while (arbre[i]>0){
        arbre2[indice]=arbre[i];
        i++;
        indice++;
    }
    arbre2[indice]=0;
    indice++;
    for (int j=indice;j<((19+2*dec)*nspec-14-dec);j++){
        arbre2[j]=0;
    }
    for (int j=0;j<((13+2*dec)*nspec-9-2*dec);j++){
        arbre[j]=arbre2[j];
    }
    nmax=nmax2;
    }
    if (ib==2000){
        emptynoeud[0]=-10;
    }
    partnerspec[(2*nspec-2)]=(2*nspec-2);
    father[(2*nspec-2)]=(2*nspec-2);
    int indtps=0;
}
    

void calculstat(int *emptynoeud,char *arbre,int **speciestott,int *speciestot,int *partnerspec,int *father,long double *temps,long double *temps2,int *tempsordre,int **descendants,long double **matricedistance,long double **matricenoeud,int *poidsnoeud,int *niveaunoeud,long double *niveaunoeuddate,long double *niveaunoeuddatenorm,int *distanceroot,int *abondnoeud,char *caractere,char *caractere2,char* tree,int** abondancest,int dec,int* indic,double* stat){
    int ntot=0;
    int nspec=indic[0];
        // tableau d'abondance des espèces et des noeuds internes
    
    for (int i=0;i<nspec;i++){
        speciestott[i][0]=abondancest[0][i];
        ntot+=speciestott[i][0];
        speciestot[i]=0;
    }
    ntot=0;
    int s1=0;
    for (int i=0;i<nspec;i++){
        speciestot[i]=speciestott[i][0];
        ntot+=speciestot[i];
        if (speciestot[i]>0){
            s1+=1;
        }
    }
    stat[4]=s1;
    double shan1=ntot*logl(ntot);
    for (int i=0;i<nspec;i++){
      if (speciestot[i]>0){
        shan1-=speciestot[i]*logl(speciestot[i]);
      }
    }
    shan1/=ntot;
    stat[5]=shan1;
    
/*
    int nnn=3*nspec+2*ntot;
    for (int i=nspec;i<(2*nspec-1);i++){
        speciestott[i][0]=0;
        speciestot[i]=0;
    }
    
    ///Lecture de l'arbre
    ///1-stockage de l'arbre en chaîne de caractères
    int iii=0;
    while (tree[iii]>0){
        arbre[iii]=tree[iii];
        iii++;
    }
    arbre[iii]=0;
    //cerr<<"lecture faite";
    
    
    for (int u=0;u<(2*nspec-1);u++){
        temps[u]=9.999;
        temps2[u]=0;
    }
    
    for (int u=0;u<(2*nspec-1);u++){
        tempsordre[u]=-1;
    }
    
    
    int nmax=((15+2*dec)*nspec-6-dec)-18-2*dec;
    //cerr<<"chargement arbre"<<endl;
    for (int i=0;i<nspec;i++){
        descendants[i][0]=i;
        descendants[i][1]=i;
    }
    for (int i=nspec;i<(2*nspec-1);i++){
        descendants[i][0]=-1;
        descendants[i][1]=-1;
    }
    for (int i=0;i<nspec;i++){
        for (int j=0;j<nspec;j++){
            matricedistance[i][j]=0;
        }
    }
    for (int i=0;i<(2*nspec-1);i++){
        for (int j=0;j<(2*nspec-1);j++){
            matricenoeud[i][j]=0;
        }
    }
    
    for (int h=0;h<nspec;h++){
        poidsnoeud[h]=0;
        emptynoeud[h]=0;
        if (speciestott[h][0]>0){
            poidsnoeud[h]=1;
        }
        niveaunoeud[h]=0;
        niveaunoeuddate[h]=0;
        niveaunoeuddatenorm[h]=0;
        distanceroot[h]=0;
    }
    for (int h=nspec;h<(2*nspec-1);h++){
        emptynoeud[h]=0;
        poidsnoeud[h]=0;
        niveaunoeud[h]=0;
        niveaunoeuddate[h]=0;
        niveaunoeuddatenorm[h]=0;
        distanceroot[h]=0;
    }
    //cerr<<"c1";
    chargementarbre(emptynoeud,abondnoeud,caractere,caractere2,partnerspec,temps,temps2,tempsordre,father,nmax,arbre,dec,nspec,descendants,matricedistance,matricenoeud,poidsnoeud,niveaunoeud,distanceroot,niveaunoeuddate,niveaunoeuddatenorm,speciestott);
    //cerr<<"c1";
    if (emptynoeud[0]==-10){
        stat[5]=0;
        stat[6]=0;
        stat[7]=0;
        stat[8]=0;
        stat[9]=0;
        stat[10]=0;
        stat[11]=0;
        stat[12]=0;
        stat[13]=0;
        stat[14]=0;
        stat[15]=0;
        stat[16]=0;
        stat[17]=0;
        stat[18]=0;
        stat[19]=0;
        stat[20]=0;
        stat[21]=0;
        stat[22]=0;
    }
    else {
    
    //CALCUL DES SUMSTATS
     
    //SAMPLE 1
    ntot=0;
    int s1=0;
    for (int i=0;i<nspec;i++){
        speciestot[i]=speciestott[i][0];
        ntot+=speciestot[i];
        if (speciestot[i]>0){
            s1+=1;
        }
    }
    
    
    
    
    
    long double distancephylomoy=0;
    long double distancephylomoyspec=0;
    long double distancenoeudmoy1=0;
    long double distancenoeudmoyspec1=0;
    
    long double distprochvoisin=0;
    long double distprochvoisinspec=0;
    long double minimatdist=0;
    long double minimatdistnoeud=0;
    long double hetero1=ntot*ntot;
    long double shan1=ntot*logl(ntot);

    
    for (int i=0;i<(nspec-1);i++){
      if (speciestot[i]>0){
        hetero1-=speciestot[i]*speciestot[i];
        shan1-=speciestot[i]*logl(speciestot[i]);
        for (int j=(i+1);j<nspec;j++){
          if (speciestot[j]>0){
            distancephylomoy+=speciestot[i]*speciestot[j]*matricedistance[i][j];
            distancephylomoyspec+=matricedistance[i][j];
            distancenoeudmoy1+=speciestot[i]*speciestot[j]*matricenoeud[i][j];
            distancenoeudmoyspec1+=matricenoeud[i][j];
          }
        }
      }
    }
    if (speciestot[(nspec-1)]>0){
        hetero1-=speciestot[(nspec-1)]*speciestot[(nspec-1)];
        shan1-=speciestot[(nspec-1)]*logl(speciestot[(nspec-1)]);
    }
    hetero1/=(ntot*ntot);
    shan1/=ntot;
    distancephylomoy/=(ntot*(ntot-1));
    distancephylomoy*=4;
    distancephylomoyspec/=(s1*(s1-1));
    distancephylomoyspec*=4;
    distancenoeudmoy1/=(ntot*(ntot-1));
    distancenoeudmoy1*=2;
    distancenoeudmoyspec1/=(s1*(s1-1));
    distancenoeudmoyspec1*=2;
    for (int i=0;i<nspec;i++){
      if (speciestot[i]>0){
        minimatdist=100;
        minimatdistnoeud=100;
        for (int j=0;j<nspec;j++){
            if ((matricedistance[i][j]<minimatdist)&&(i!=j)&&(speciestot[j]>0)){
                minimatdist=matricedistance[i][j];
            }
        }
        distprochvoisin+=speciestot[i]*minimatdist;
        distprochvoisinspec+=minimatdist;
      }
    }
    distprochvoisin/=ntot;
    distprochvoisinspec/=s1;
    distprochvoisin*=2;
    distprochvoisinspec*=2;
    long double distancephylomoynorm1=distancephylomoy*2.0/MAXDIST;
    long double distancephylomoyspecnorm1=distancephylomoyspec*2.0/MAXDIST;
    long double distprochvoisinnorm1=distprochvoisin*2.0/MAXDIST;
    long double distprochvoisinspecnorm1=distprochvoisinspec*2.0/MAXDIST;
    
    long double vardistphylo1=0;
    long double vardistphylospec1=0;
    long double vardistnoeud1=0;
    long double vardistnoeudspec1=0;
    for (int i=0;i<(nspec-1);i++){
      if (speciestot[i]>0){
        for (int j=(i+1);j<nspec;j++){
          if (speciestot[j]>0){
            vardistphylo1+=speciestot[i]*speciestot[j]*(matricedistance[i][j]-distancephylomoy)*(matricedistance[i][j]-distancephylomoy);
            vardistphylospec1+=(matricedistance[i][j]-distancephylomoyspec)*(matricedistance[i][j]-distancephylomoyspec);
            vardistnoeud1+=speciestot[i]*speciestot[j]*(matricenoeud[i][j]-distancenoeudmoy1)*(matricenoeud[i][j]-distancenoeudmoy1);
            vardistnoeudspec1+=(matricenoeud[i][j]-distancenoeudmoyspec1)*(matricenoeud[i][j]-distancenoeudmoyspec1);
          }
        }
      }
    }
    
    vardistphylo1/=(ntot*(ntot-1));
    vardistphylo1*=4;
    vardistphylospec1/=(s1*(s1-1));
    vardistphylospec1*=4;
    vardistnoeud1/=(ntot*(ntot-1));
    vardistnoeud1*=2;
    vardistnoeudspec1/=(s1*(s1-1));
    vardistnoeudspec1*=2;
    
    long double vardistphylonorm1=vardistphylo1*4.0/(MAXDIST*MAXDIST);
    long double vardistphylospecnorm1=vardistphylospec1*4.0/(MAXDIST*MAXDIST);
    
    long double vardistvoisin1=0;
    long double vardistvoisinspec1=0;
    for (int i=0;i<nspec;i++){
      if (speciestot[i]>0){
        minimatdist=100;
        minimatdistnoeud=100;
        for (int j=0;j<nspec;j++){
            if ((matricedistance[i][j]<minimatdist)&&(i!=j)&&(speciestot[j]>0)){
                minimatdist=matricedistance[i][j];
            }
        }
        vardistvoisin1+=speciestot[i]*(minimatdist-distprochvoisin)*(minimatdist-distprochvoisin);
        vardistvoisinspec1+=(minimatdist-distprochvoisinspec)*(minimatdist-distprochvoisinspec);
       }
    }
    
    vardistvoisin1/=ntot;
    vardistvoisinspec1/=s1;
    vardistvoisin1*=2;
    vardistvoisinspec1*=2;
    
    long double vardistvoisinnorm1=vardistvoisin1*4.0/(MAXDIST*MAXDIST);
    long double vardistvoisinspecnorm1=vardistvoisinspec1*4.0/(MAXDIST*MAXDIST);
    
    double varianceni1=0;
    if (nspec>1){
        double Nmean=double(ntot)/double(s1);
        for (int i=0;i<nspec;i++){
         if (speciestot[i]>0){
            varianceni1+=(speciestot[i]-Nmean)*(speciestot[i]-Nmean);
         }
        }
        varianceni1/=(s1-1);
    }
    
    int ntot1=ntot;

    
    //stats d'ensemble

    long double bun=0;
    int h=nspec;
    int hh=0;
    while ((h<(2*nspec-2))&&(hh<(s1-2))){
        if ((niveaunoeud[h]>0)&&(emptynoeud[h]>0)){
            bun+=1/double(niveaunoeud[h]);
            hh++;
        }
        h++;
    }
       
    stat[5]=bun;
    stat[6]=niveaunoeuddatenorm[(2*nspec-2)];
    stat[7]=s1;
    stat[8]=shan1;
    stat[9]=hetero1;
    stat[10]=varianceni1;
    stat[11]=distprochvoisinnorm1;
    stat[12]=distprochvoisinspecnorm1;
    stat[13]=distancephylomoynorm1;
    stat[14]=distancephylomoyspecnorm1;
    stat[15]=distancenoeudmoy1;
    stat[16]=distancenoeudmoyspec1;
    stat[17]=vardistvoisinnorm1;
    stat[18]=vardistvoisinspecnorm1;
    stat[19]=vardistphylonorm1;
    stat[20]=vardistphylospecnorm1;
    stat[21]=vardistnoeud1;
    stat[22]=vardistnoeudspec1;
} 
*/
}

int position(int32 *comsample,int *abondancest,int u,int indic){
    int res=-1;
    int p=0;
    int uu=u;
    while (res=-1){
        if (uu<(abondancest[p]-comsample[p])){
            res=p;
        }
        else{
            uu-=(abondancest[p]-comsample[p]);
            p++;
        }
    }
return res;
}

void initialisation(int32 *comsample,int *abondancest,int JJsample,int indic){
    for (int i=0;i<JJsample;i++){
        int u=rg.IRandom (0, (indic-1));
        int p=position(comsample,abondancest,u,indic);
        comsample[p]++;
    }   
}

StochasticLib1 rmulti(time(0)+527);

void forwarddensitedep(double *probamultinom,int32 *comsample,int32 *comsample2,int *abondancest,int JJsample,int indic,double alpha,double m){
    //initialisation
    double somme=0.0;
    for (int i=0;i<indic;i++){
        probamultinom[i]=comsample[i]+0.0;
        somme+=comsample[i];
    }
    if (somme==0.0){
        for (int i=0;i<indic;i++){
            probamultinom[i]+=abondancest[i]+0.0;
        }
        int32 n=JJsample;
        rmulti.Multinomial (comsample, probamultinom, n, indic);
        for (int i=0;i<indic;i++){
            probamultinom[i]=comsample[i]+0.0;
        }
    }
    
  for (int j=0;j<100;j++){
   
    //mortalité densité dependante de 1/10ème de la communauté locale
    int32 n=JJsample/10;
    for (int i=0;i<indic;i++){
        probamultinom[i]=double(comsample[i])/(JJsample+0.0);
        if (probamultinom[i]>0){
            probamultinom[i]=expl((-alpha+1)*logl(probamultinom[i]));
        }
    }
    rmulti.Multinomial (comsample2, probamultinom, n, indic);
    
    somme=0.0;
    n=0;
    for (int i=0;i<indic;i++){
        if (comsample[i]>comsample2[i]){
            comsample[i]-=comsample2[i];
            n+=comsample2[i];
        }
        else { 
            n+=comsample[i];
            comsample[i]=0;
        }
        somme+=comsample[i]+0.0;
    }
    //recrutement neutre (non densite dep) de n individus
    
    for (int i=0;i<indic;i++){
        probamultinom[i]=double(abondancest[i]);
    }
    rmulti.Multinomial (comsample2, probamultinom, n, indic);
    for (int i=0;i<indic;i++){
        comsample[i]+=comsample2[i];
    }
  }
}

void densitedependance(double *probamultinom,int *abondmetacom,int32 *comsample,int32 *comsample2,int *abondancest,double alpha,int *indic,int JJsample,double m){
    for (int i=0;i<2000;i++){
        abondmetacom[i]=abondancest[i];
    }
    for (int i=0;i<2000;i++){
        comsample[i]=0;
        probamultinom[i]=0.0;
    }
    //initialisation(comsample,abondancest,JJsample,indic[0]);
    forwarddensitedep(probamultinom,comsample,comsample2,abondancest,JJsample,indic[0],alpha,m);
    indic[1]=0;
    for (int i=0;i<indic[0];i++){
        abondancest[i]=int(comsample[i]);
        if (abondancest[i]>0){
            indic[1]+=1;
        }
    }
    
}

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

#define MAX_IT      1000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* explansion coefficient */

int tototo;

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
  long double fe;      /* value of function at explansion point */
  long double fc;      /* value of function at contraction point */
  long double *vr;     /* reflection - coordinates */
  long double *ve;     /* explansion - coordinates */
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
      s += powl((f[j]-favg),(long double) 2.0)/(n);
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
  
  for (j=0;j<2;j++) {
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

int rounding(double a){
	double b=floor(a);
	if (a-b>=0.5){
		b=b+1;
	}
return int(b);
}

/* MAIN ROUTINE */
/****************/
int main(int argc,char **argv) {
cerr << "ESTIMATING PARAMETER DELTA FOR THE DISPERSAL-LIMITED NON-NEUTRAL MODEL BY ABC"<<endl;
cerr << "This program estimates the parameter delta of deviation from neutrality using the approach of Jabot and Chave (2011). For more details, see the manual."<<endl;
cerr << "Options for answering yes/no questions: 1, y, or yes for 'yes'; 0, n, or no for 'no'.\n \n";

//Lecture data format tetame + Calcul sumstat data
  
char bufi[128];
cerr << "Please enter the data file name (without'.txt') ";
cin >> bufi;
char nomfi[256];
sprintf(nomfi,"%s.txt",bufi);
char nomfo[256];
char nomfou[256];
sprintf(nomfou,"%s_out.txt",bufi); 
char nomfoun[256];
sprintf(nomfoun,"%s_test-neutrality_out.txt",bufi); 
char nomfoR[256];

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
    int i;
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
    int **Abund;
    Abund = new int*[nbsamples];
    double *Shannondata;
    Shannondata= new double[nbsamples];
    int *Jdata;
    Jdata = new int[nbsamples];  
    int Jmaxdata=0;                       
    ifstream inf2(nomfi);
    int sample;
    int s;
    for(sample=0;sample<nbsamples;sample++)
        Species[sample]=Species0[sample];                                 
    for(sample=0;sample<nbsamples;sample++){
        Abund[sample] = new int[Species[sample]];
        int nspecpres=0;
        Jdata[sample]=0;
        Shannondata[sample]=0;
        for(s=0;s<Species[sample];s++) {
            inf2 >> Abund[sample][s];
            if (Abund[sample][s]>0){
                nspecpres++;
                Jdata[sample]+=Abund[sample][s];
                Shannondata[sample]-=Abund[sample][s]*logl(Abund[sample][s]);
            }
        }
        if (Jdata[sample]>0){
            Shannondata[sample]+=Jdata[sample]*logl(Jdata[sample]);
            Shannondata[sample]/=Jdata[sample];
        }
        Species[sample]=nspecpres;
        cerr << "In sample " << sample+1<< ", number of individuals: " << Jdata[sample]<< ", number of species: " << Species[sample]<< ", Shannon's index: " <<Shannondata[sample] <<endl;
        inf2 >> buffer;
        if (Jdata[sample]>Jmaxdata){
            Jmaxdata=Jdata[sample];
        }
        qsort(Abund[sample],Species[sample],sizeof(int),intcompare2);  
    }
    inf2.close();
    
  
  int neutraltestquestion=0;
  char inineutr[128];
  int endquestion=0;
  while (endquestion==0){
  cerr << "Do you want to perform the neutrality test in your sample(s)? ";
  cin >> inineutr;
  if (!strcmp(inineutr,"1")||!strcmp(inineutr,"y")||!strcmp(inineutr,"yes")) {
    endquestion=1;
    neutraltestquestion=1;
  }
  else{
    if (!strcmp(inineutr,"0")||!strcmp(inineutr,"n")||!strcmp(inineutr,"no")) {
        endquestion=1;
    }
  }
  }

  int deltaquestion=0;
  char inidelta[128];
  endquestion=0;
  while (endquestion==0){
  cerr << "Do you want to infer Delta value(s) in your sample(s)? ";
  cin >> inidelta;
  if (!strcmp(inidelta,"1")||!strcmp(inidelta,"y")||!strcmp(inidelta,"yes")) {
    endquestion=1;
    deltaquestion=1;
  }
  else{
    if (!strcmp(inidelta,"0")||!strcmp(inidelta,"n")||!strcmp(inidelta,"no")) {
        endquestion=1;
    }
  }
  }
  

    //Remplir tableaux J_a_faire, S_a_faire
    int nJ=nbsamples;
    int *Jafaire=new int[nJ];
    int *Safaire=new int[nJ];
    for (int u=0;u<nJ;u++){
        Jafaire[u]=Jdata[u];
        Safaire[u]=Species[u];
    }
    
    double theta,I;
    int dec,nbsamp,nbsamptot,JJ1,JJ2;
    dec=1;
    int nbsamp2=1;
    int speciestarget;
    int speciestarget2;
    
    double thetamin;
    double thetamax;
    double mmin;
    double mmax;
    double deltamin;
    double deltamax;
    int JJsample;
    double mforward;
    int speciestargetmax;
    double percentsimul;
    int nsimulneutraltest;
        
    ifstream ini("prior-distributions_ABC.txt");
    ini>>nomfo;
    ini>>thetamin;
    ini>>nomfo;
    ini>>thetamax;
    ini>>nomfo;
    ini>>mmin;
    ini>>nomfo;
    ini>>mmax;
    ini>>nomfo;
    ini>>deltamin;
    ini>>nomfo;
    ini>>deltamax;
    ini>>nomfo;
    ini>>JJ1;
    if (JJ1<3*Jmaxdata){
        JJ1=3*Jmaxdata;
    }
    ini>>nomfo;
    ini>>nbsamptot;
    ini>>nomfo;
    ini>>percentsimul;
    ini>>nomfo;
    ini>>nsimulneutraltest;
    
    speciestarget=3;
    speciestargetmax=1197;
    mforward=1;
    ini.close();
    
    //Stuff for random number generator
    int32 seed = (int32)time(0);  
    int seedmulti = (int)time(0)-100; 
    CRandomMersenne rg(seed);
    StochasticLib1 rmulti(seedmulti);
    
    int tttt=int(floor(thetamin));
    int iiii=int(floor(mmin));
    
    int tttt2=int(floor(thetamax));
    int iiii2=int(floor(mmax));
    int mf=int(floor(10*mforward));    
    
    char *tree;
    tree = new char[((19+2*10)*(2000)-14-2*10)];
    int **abondancest;
    abondancest=new int*[1];
    abondancest[0]=new int[2000];
    int *indic; // S, A1, A2
    indic= new int[3];
    double *stat;
    stat = new double[23];
    for (int ii=0;ii<23;ii++){
        stat[ii]=0;
    }
    int testsimul=0;
    
        int **spoidsindividus1;
        spoidsindividus1= new int*[JJ1];
        for (int ii=0;ii<JJ1;ii++){
            spoidsindividus1[ii]=new int[2];
            spoidsindividus1[ii][0]=1;
            spoidsindividus1[ii][1]=1;
        }
        int **spoidsancetres1;
        spoidsancetres1= new int*[(JJ1)];
        for (int ii=0;ii<(JJ1);ii++){
            spoidsancetres1[ii]= new int[1];
            spoidsancetres1[ii][0]=0;
        }
        int *sabondances;
        sabondances= new int[JJ1];
        int *spositions;
        spositions= new int[JJ1];
        double *stemps;
        stemps= new double[JJ1];
        char *streef= new char[((19+2*dec)*2000-14-2*dec)];
        char *scara= new char[(5+dec)];
        int **sresult;
        sresult = new int*[1];
        sresult[0]= new int[2000];
    
    int **cspeciestott= new int*[2*2000-1];
    int *cspeciestot= new int[2*2000-1];
    for (int i=0;i<(2*2000-1);i++){
        cspeciestott[i]=new int[2];
    }
    char *carbre=new char[((19+2*dec)*2000-14-dec)];
    int *cpartnerspec;                     //   tableau des partenaires de spéciation
    cpartnerspec= new int[2*2000-1];
    int *cfather;                     //   tableau des noeuds parents
    cfather= new int[2*2000-1];
    long double *ctemps;
    ctemps= new long double[2*2000-1];    // tableau des temps de spéciation
    long double *ctemps2;
    ctemps2= new long double[2*2000-1];    // tableau des temps de spéciation
    int *ctempsordre;
    ctempsordre= new int[2*2000-1];   // tableau des ordres de spéciation
    int **cdescendants;
    cdescendants= new int*[(2*2000-1)];
    for (int i=0;i<2000;i++){
        cdescendants[i]=new int[2];
        cdescendants[i][0]=i;
        cdescendants[i][1]=i;
    }
    for (int i=2000;i<(2*2000-1);i++){
        cdescendants[i]=new int[2];
        cdescendants[i][0]=-1;
        cdescendants[i][1]=-1;
    }
    long double **cmatricedistance;
    cmatricedistance= new long double *[2000];
    for (int i=0;i<2000;i++){
        cmatricedistance[i]=new long double[2000];
        for (int j=0;j<2000;j++){
            cmatricedistance[i][j]=0;
        }
    }
    long double **cmatricenoeud;
    cmatricenoeud= new long double *[(2*2000-1)];
    for (int i=0;i<(2*2000-1);i++){
        cmatricenoeud[i]=new long double[(2*2000-1)];
        for (int j=0;j<(2*2000-1);j++){
            cmatricenoeud[i][j]=0;
        }
    }
    
    int *cpoidsnoeud;
    cpoidsnoeud= new int[(2*2000-1)];
    int *cniveaunoeud;
    cniveaunoeud= new int[(2*2000-1)];
    long double *cniveaunoeuddate;
    cniveaunoeuddate= new long double[(2*2000-1)];
    long double *cniveaunoeuddatenorm;
    cniveaunoeuddatenorm= new long double[(2*2000-1)];
    int *cdistanceroot;
    cdistanceroot= new int[(2*2000-1)];
    for (int h=0;h<2000;h++){
        cpoidsnoeud[h]=0;
        cniveaunoeud[h]=0;
        cniveaunoeuddate[h]=0;
        cniveaunoeuddatenorm[h]=0;
        cdistanceroot[h]=0;
    }
    for (int h=2000;h<(2*2000-1);h++){
        cpoidsnoeud[h]=0;
        cniveaunoeud[h]=0;
        cniveaunoeuddate[h]=0;
        cniveaunoeuddatenorm[h]=0;
        cdistanceroot[h]=0;
    }
    
    int *cabondnoeud= new int[(2*2000-1)];
    int *cemptynoeud= new int[(2*2000-1)];
    char *ccaractere; 
    ccaractere=new char[(5+dec)];
    char *ccaractere2; 
    ccaractere2=new char[(5+dec)];
    
    
    
    int *dabondmetacom=new int [2000];
    for (int i=0;i<2000;i++){
        dabondmetacom[i]=0;
    }
    int32 *dcomsample= new int32[2000];
    int32 *dcomsample2= new int32[2000];
    double *dprobamultinom = new double[2000];
    for (int i=0;i<2000;i++){
        dcomsample[i]=0;
        dcomsample2[i]=0;
        dprobamultinom[i]=0.0;
    }
    
    
    int i2;
    int *nbJS;
    nbJS=new int[nJ];
    
    
   //
    
    
    //NEUTRALITY TEST
    if (neutraltestquestion==1){

        double **shanref;
        shanref = new double *[nJ];
        for (int u=0;u<nJ;u++){
            shanref[u]=new double[nsimulneutraltest];
            for (int v=0;v<nsimulneutraltest;v++){
                shanref[u][v]=0;
            }
        }
    
        //NEUTRAL TEST AVEC OUTPUT ABC SIMUL
            cerr<<endl;
            cerr<<"NEUTRALITY TEST ..."<<endl;
            
            ofstream outdatan(nomfoun);  
            outdatan<<"Sample\t J\t S\t H\t pvalue\n";
            outdatan.flush();
            
            for (int u=0;u<nJ;u++){
                sprintf(nomfo,"ABC_neutral_J%d_S%d.txt",Jafaire[u],Safaire[u]);
                ifstream inf3;                                                         
                inf3.open(nomfo);
                if(!inf3){
                    ofstream outt(nomfo,ios::app);
                    outt<<"J\t Theta\t m\t S\t H\n";
                    outt.flush();
                    outt.close();
                    nbJS[u]=0;  
                }
                else {
                    inf3.close();
                    int testJS=0;
                    nbJS[u]=2;
                    while (testJS==0){
                        sprintf(nomfo,"ABC_neutral_J%d_S%d_%d.txt",Jafaire[u],Safaire[u],nbJS[u]);
                        inf3.open(nomfo);
                        if (!inf3){
                            ofstream outt(nomfo,ios::app);
                            outt<<"J\t Theta\t m\t S\t H\n";
                            outt.flush();
                            outt.close();
                            testJS=1;
                        }
                        else{
                            nbJS[u]++;
                        }
                    }
                }
            }
            
            
        for (int u=0;u<nJ;u++){
            if (nbJS[u]==0){
                sprintf(nomfo,"ABC_neutral_J%d_S%d.txt",Jafaire[u],Safaire[u]);
            }
            else {
                sprintf(nomfo,"ABC_neutral_J%d_S%d_%d.txt",Jafaire[u],Safaire[u],nbJS[u]);
            }
            ofstream outt(nomfo,ios::app); 
            //TETAME
            sample=u;
                // Total number of individuals
    J=0;
    SPP = Species[sample];
    for(s=0;s<SPP;s++)
        J += Abund[sample][s];
    cerr << "SAMPLE "<<(u+1)<<": Number of individuals: "<< J << endl;

    int MaxA = 0;
    MaxA=Abund[sample][(int)SPP-1];
    //cerr <<" "<<endl;
    //cerr << "Sample " << 1+sample << endl;
    //cerr << "Maximal abundance: " << MaxA << endl;

    // abundance distribution
    int *Phi = new int[MaxA+1];
    for(s=0;s<=MaxA;s++) Phi[s]=0;
    for(s=0;s<SPP;s++) Phi[Abund[sample][s]]++;

    // Number of distinct abundances
    int NDA=0;
    for(s=0;s<=MaxA;s++) if(Phi[s] > 0) {NDA++;}
    //cerr << "NbDistinctAbund\t" << NDA << endl;
    
    //cerr << "Start computing Stirling numbers ...\n";
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
    int n,im;
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

    //cerr << "Start computing ln(K(D,A)) ...\n";
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
	    K[i] = logl(K[i]);                                    // now K[A]=ln(K(D,A)) in Etienne's paper
    }
    for(i=0;i<NDA;i++) delete[] T[i];
	delete[] T;
    delete[] poly2;
    delete[] f;
    delete[] g;

// search of "infinite" values in K[A]
    int borneinf=int(SPP-1);
    int bornesup=int(J+1);
    long double maxlogl=11333.2;
    int infinity=0;
    for(i=int(SPP);i<=J;i++){
        if ((K[i]>maxlogl)||(K[i]<-maxlogl)) {
            infinity=1;
	        break;
        }
	    borneinf++;
    }  //after that, borneinf=indice next to infinity but before
    for(i=0;i<=J-SPP;i++){
        if ((K[(int)J-i]>maxlogl)||(K[(int)J-i]<-maxlogl)) {
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

    // THIRD STEP: define the logl-Likelihood
    // L(theta,I) = theta^S/(I)_J * sum_A K(D,A) I^A/(theta)_A
    // logl(theta,I) = S*logl(theta)-logl((I)_J)  
    //                  + logl(sum_A K(D,A) I^A/(theta)_A)
    // where (x)_N = x(x+1)...(x+N-1)
    // I = m(J-1)/(1-m)
    //
    // theta is between 1 and S
    // m is between 0 and 1

    //cerr << "Compute the Ewens theta and logl-likelihood ...\n";
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
  //cerr << "Ewens' -logl-likelihood: " <<  LEWENS << endl;
  long double x[2];

#if MPI0
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  //cerr << "Maximizing the likelihood ...\n";
  long double start[2] = {EWENS_THETA,J/10};
  long double min,minb;
  long double startb[2] = {J/10,EWENS_THETA};
  startb[0]=start[1];
  startb[1]=start[0];
  //cerr<<"factor="<<factor<<endl;
  //cerr<<"expl(5)="<<expl(5)<<endl;
  //cerr<<"logl(5)="<<logl(5)<<endl;
  //cerr<<"lpochham(5,2)="<<lpochham(5,2)<<endl;
  //cerr<<"llik="<<llik(start)<<endl;
  min=simplex(&llik,start,5,1.0e-13,1,LEWENS,EWENS_THETA);  //launching of the simplex algorithm           
  minb=simplex(&llik,startb,5,1.0e-13,1,LEWENS,EWENS_THETA);  //launching of the simplex algorithm
  double prob1;
  double THETA_ETIENNE1;
  double I_ETIENNE1;   
  double THETA_ETIENNE2;
  double I_ETIENNE2; 
  //cerr<<logl(start[0])<<" "<<logl(start[1])<<" "<<logl(startb[0])<<" "<<logl(startb[1])<<endl;  

if (min<0){
    //cerr<<"case e";
    THETA_ETIENNE1=startb[0];
    I_ETIENNE1=startb[1];
    prob1=1;
}
else {
 if (minb<0){
    //cerr<<"case f";
    THETA_ETIENNE1=start[0];
    I_ETIENNE1=start[1];
    prob1=1;
 }
 else {

if (min>minb){
    if (min<(minb+10)){
        //cerr<<"case a";
        THETA_ETIENNE1=startb[0];
        I_ETIENNE1=startb[1];
        THETA_ETIENNE2=start[0];
        I_ETIENNE2=start[1];
        prob1=expl(-minb)/(expl(-minb)+expl(-min));
    }
    else{
        //cerr<<"case b";
        THETA_ETIENNE1=startb[0];
        I_ETIENNE1=startb[1];
        prob1=1;
    }
}
else{
    if (minb<(min+10)){
        //cerr<<"case c";
        THETA_ETIENNE1=start[0];
        I_ETIENNE1=start[1];
        THETA_ETIENNE2=startb[0];
        I_ETIENNE2=startb[1];
        prob1=expl(-min)/(expl(-min)+expl(-minb));
    }
    else{
        //cerr<<"case d";
        THETA_ETIENNE1=start[0];
        I_ETIENNE1=start[1];
        prob1=1;
    }
}
}
}
            
//END OF TETAME
            
            i=0;
            i2=0;
            
            JJsample=Jafaire[u];
            double theta=THETA_ETIENNE1;
            double I1=I_ETIENNE1;
            double m1=I1/(I1+JJsample-1);
            int nsimul1=rounding(nsimulneutraltest*prob1);
            //cerr<<"-logl-likelihood first peak="<<min<<" -logl-likelihood second peak="<<minb<<endl;
            cerr<<"Composition of neutral simulations:"<<endl;
            cerr<<nsimul1<<" simulations with parameter Theta="<<THETA_ETIENNE1<<" and parameter I="<<I_ETIENNE1<<endl;
            if ((nsimulneutraltest-nsimul1)>0){
                cerr<<(nsimulneutraltest-nsimul1)<<" simulations with parameter Theta="<<THETA_ETIENNE2<<" and parameter I="<<I_ETIENNE2<<endl;
            }
            cerr<< "Neutral Simulations running..."<<endl;
            while(i<(nsimul1)){
                    
                    testsimul=0;
                    while (testsimul==0){
                        i2++;
                        indic[0]=0;
                        indic[1]=0;
                        indic[2]=0;
                        //cerr<<"b"<<endl;
                        testsimul=simularbrebis(spoidsindividus1,stemps,theta,I1,JJsample,abondancest,indic);
                        //cerr<< indic[0] << " ";
                        if ((testsimul==1)&&(indic[0]==Safaire[u])){
                                stat[0]=JJsample;
                                stat[1]=theta;
                                stat[2]=m1;
                                stat[3]=indic[0];
                                for (int j=0;j<4;j++){
                                    outt<<stat[j]<<"\t";
                                }
                                stat[5]=JJsample*logl(JJsample);
                                for (int j=0;j<indic[0];j++){
                                    if (abondancest[0][j]>0){
                                        stat[5]-=abondancest[0][j]*logl(abondancest[0][j]);
                                    }
                                }
                                stat[5]/=JJsample;
                                shanref[u][i]=stat[5];
                                outt<<stat[5]<<endl;
                        }
                        else{
                            testsimul=0;
                        }
                    }
                    outt.flush();
                i++;
                if (i%100==0){
                    cerr<<i<<" ";
                }
            }
            
            theta=THETA_ETIENNE2;
            I1=I_ETIENNE2;
            m1=I1/(I1+JJsample-1);
            
            while(i<nsimulneutraltest){
                    //JJsample=Jafaire[u];
                    testsimul=0;
                    while (testsimul==0){
                        i2++;
                        indic[0]=0;
                        indic[1]=0;
                        indic[2]=0;
                        testsimul=simularbrebis(spoidsindividus1,stemps,theta,I1,JJsample,abondancest,indic);
                        if ((testsimul==1)&&(indic[0]==Safaire[u])){
                                stat[0]=JJsample;
                                stat[1]=theta;
                                stat[2]=m1;
                                stat[3]=indic[0];
                                for (int j=0;j<4;j++){
                                    outt<<stat[j]<<"\t";
                                }
                                stat[5]=JJsample*logl(JJsample);
                                for (int j=0;j<indic[0];j++){
                                    if (abondancest[0][j]>0){
                                        stat[5]-=abondancest[0][j]*logl(abondancest[0][j]);
                                    }
                                }
                                stat[5]/=JJsample;
                                shanref[u][i]=stat[5];
                                outt<<stat[5]<<endl;
                        }
                        else{
                            testsimul=0;
                        }
                    }
                    outt.flush();
                i++;
                if (i%100==0){
                    cerr<<i<<" ";
                }
            }
            outt.close();
            cerr<<endl;
            cerr<<endl;
            double pval=0;
            for (int v=0;v<nsimulneutraltest;v++){
                if (shanref[u][v]<Shannondata[u]){
                    pval++;
                }
            }
            pval/=double(nsimulneutraltest);
            outdatan<<(u+1)<<"\t"<<Jdata[u]<<"\t"<<Species[u]<<"\t"<<Shannondata[u]<<"\t"<<pval<<endl;
            outdatan.flush();
            
        }
        outdatan.close();
    }
    
    
        //DELTA INFERENCE
    if (deltaquestion==1){
    
    sprintf(nomfoR,"%s_outR.txt",bufi);
    ofstream outR(nomfoR);
    
outR<<"euclideandistance<-function(statdata,statsimul,nstat,nparam,varsimul){"<<endl;
outR<<"result=0"<<endl;
outR<<"for (i in 1:nstat){"<<endl;
outR<<"result=result+(statsimul[nparam+i]-statdata[i])*(statsimul[nparam+i]-statdata[i])/varsimul[i]"<<endl;
outR<<"}"<<endl;
outR<<"sqrt(result)"<<endl;
outR<<"}"<<endl;

outR<<"calculpoids<-function(delta,distance,nbsimul){"<<endl;
outR<<"result=array(1,nbsimul)"<<endl;
outR<<"for (i in 1:nbsimul){"<<endl;
outR<<"result[i]=max(0,result[i]-((distance[i]/delta)*(distance[i]/delta)))"<<endl;
outR<<"}"<<endl;
outR<<"result"<<endl;
outR<<"}"<<endl;

outR<<"selectionsimulation<-function(statdata,statsimul,nstat,nparam,nbsimulutiles){"<<endl;
outR<<"d=dim(statsimul)[1]"<<endl;
outR<<"distance=array(0,d)"<<endl;
outR<<"varsimul=array(0,nstat)"<<endl;
outR<<"for (i in 1:nstat){"<<endl;
outR<<"varsimul[i]=var(statsimul[,nparam+i])"<<endl;
outR<<"}"<<endl;
outR<<"for (i in 1:d){"<<endl;
outR<<"distance[i]=as.numeric(euclideandistance(statdata,statsimul[i,],nstat,nparam,varsimul))"<<endl;
outR<<"}"<<endl;
outR<<"ordredistance<-order(distance)"<<endl;
outR<<"ordredistanceutile=ordredistance[1:nbsimulutiles]"<<endl;
outR<<"statsimulutile=distance[ordredistanceutile]"<<endl;
outR<<"for (i in 1:(nparam+nstat)){"<<endl;
outR<<"statsimulutile=cbind(statsimulutile,statsimul[,i][ordredistanceutile])"<<endl;
outR<<"}"<<endl;
outR<<"statsimulutile"<<endl;
outR<<"}"<<endl;

outR<<"statdata=read.table('"<<nomfou<<"',h=T)"<<endl;
outR<<"deltas=array(0,"<<nJ<<")"<<endl;
outR<<"for (i in 1:"<<nJ<<"){"<<endl;
outR<<"tableau=read.table(as.character(statdata[i,5]),h=T)"<<endl;
outR<<"statsimul=cbind(tableau$delta,tableau$S,tableau$H)"<<endl;
outR<<"npoints=ceiling(dim(statsimul)[1]*"<<percentsimul<<")"<<endl;
outR<<"points=selectionsimulation(c(statdata[i,3],statdata[i,4]),statsimul,2,1,npoints)"<<endl;
outR<<"poids=calculpoids(points[npoints,1],points[,1],npoints)"<<endl;
outR<<"dd=density(points[,2],bw=0.2,weights=poids/sum(poids))"<<endl;
outR<<"plot(dd)"<<endl;
outR<<"deltas[i]=dd$x[dd$y==max(dd$y)]"<<endl;
outR<<"}"<<endl;
outR<<"results=cbind(statdata,deltas)"<<endl;
outR<<"write.table(results,file='"<<bufi<<"_out2.txt',row.names=F,quote=F)"<<endl;
    outR.flush();
    outR.close();
    
    
    ofstream outdata(nomfou);  
    outdata<<"Sample\t J\t S\t H\t ABC_file_name\n";
    for (int u=0;u<nJ;u++){
        sprintf(nomfo,"ABC_J%d_S%d.txt",Jafaire[u],Safaire[u]);
        ifstream inf3;                                                         
        inf3.open(nomfo);
        if(!inf3){
            ofstream outt(nomfo,ios::app);
            outt<<"J\t Theta\t m\t delta\t S\t H\n";
            outt.flush();
            outt.close();
            outdata<<(u+1)<<"\t"<<Jdata[u]<<"\t"<<Species[u]<<"\t"<<Shannondata[u]<<"\t"<<nomfo<<endl;
            nbJS[u]=0;
            
        }
        else {
            inf3.close();
            int testJS=0;
            nbJS[u]=2;
            while (testJS==0){
                sprintf(nomfo,"ABC_J%d_S%d_%d.txt",Jafaire[u],Safaire[u],nbJS[u]);
                inf3.open(nomfo);
                if (!inf3){
                    outdata<<(u+1)<<"\t"<<Jdata[u]<<"\t"<<Species[u]<<"\t"<<Shannondata[u]<<"\t"<<nomfo<<endl;
                    ofstream outt(nomfo,ios::app);
                    outt<<"J\t Theta\t m\t delta\t S\t H\n";
                    outt.flush();
                    outt.close();
                    testJS=1;
                }
                else{
                    nbJS[u]++;
                }
            }
        }
    }
    outdata.flush();
    outdata.close();
    
    i=0;
    i2=0;
    cerr<<endl;
    cerr<<"DELTA INFERENCE ..."<<endl;
    while(i<nbsamptot){
    for (int u=0;u<nJ;u++){
        if (nbJS[u]==0){
            sprintf(nomfo,"ABC_J%d_S%d.txt",Jafaire[u],Safaire[u]);
        }
        else {
            sprintf(nomfo,"ABC_J%d_S%d_%d.txt",Jafaire[u],Safaire[u],nbJS[u]);
        }
        ofstream outt(nomfo,ios::app);
        JJsample=Jafaire[u];
        testsimul=0;
        while (testsimul==0){
        i2++;
        double theta=expl(thetamin+(thetamax-thetamin)*rg.Random());
        double I1=expl(mmin+(mmax-mmin)*rg.Random());
        double m1=I1/(I1+JJ1-1);
        double delta=deltamin+(deltamax-deltamin)*rg.Random();
        
        for (int k=0;k<2000;k++){
            abondancest[0][k]=0;
        }
        
        indic[0]=0;
        indic[1]=0;
        indic[2]=0;
        
        //cerr<<"s";
        testsimul=simularbrebis(spoidsindividus1,stemps,theta,I1,JJ1,abondancest,indic);
        //cerr<<"s";
        if (testsimul==1){
        
        if ((indic[0]>(speciestarget))&&(indic[0]<(speciestargetmax))){
            densitedependance(dprobamultinom,dabondmetacom,dcomsample,dcomsample2,abondancest[0],delta,indic,JJsample,mforward);
            if (indic[1]>1){
                
                stat[0]=JJsample;
                stat[1]=theta;
                stat[2]=m1;
                stat[3]=delta;
                calculstat(cemptynoeud,carbre,cspeciestott,cspeciestot,cpartnerspec,cfather,ctemps,ctemps2,ctempsordre,cdescendants,cmatricedistance,cmatricenoeud,cpoidsnoeud,cniveaunoeud,cniveaunoeuddate,cniveaunoeuddatenorm,cdistanceroot,cabondnoeud,ccaractere,ccaractere2,tree,abondancest,dec,indic,stat);
                for (int j=0;j<6;j++){
                    outt<<stat[j]<<"\t";
                }
                outt<<endl;
            }
        }
        else{
            testsimul=0;
        }
        }
        }
    outt.flush();
    outt.close();
    }
            i++;
            if (i%100==0){
                cerr<<i<<" ";
            }
    }
    
    }
 
    
    
    
    
    
    return 1;
}
