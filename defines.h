//------Prevent multiple includes of de.h----------------------------------
#ifndef _DE_H
#define _DE_H

//------General constants--------------------------------------------------
#define MAXDIM  3        // maximum number of dimensions i.e. parameters.
#define MAXPOP  40       // number of random vectors to be stored. Watch
// out! gi_D must be <= 33 because w_index writes
// up to location 3*gi_D-1.
#define MAXCOST   1        // maximum number of objectives to be minimized
#define MAXCONST  1        // maximum number of constraints

typedef double floatT;   ///floating point type. can be float, double, currently doesn't support long double and other types, due to Matlab engine limitation.

#define PRECISION 10    ///precision of the parameters

#define BUFSIZE 512  ///for file names

///file names
#define FEBio_FILENAME "./27Nov.feb"   ///this is the name of the .feb file. The program won't change this file. Instead it creates copies named ./26Aug_2_0, ./26Aug_2_1 for different threads.

#define MATLAB_OUTPUT "./Matlab_output.dat"
#define IDEALIZED_DATA "./idealized_data.log"

///the directory of FEBio, on Yu's iMac
#ifdef yuMac
#define RUN_FEBIO_COMMAND "/Applications/FEBio2.3.0/bin/FEBio2 -nosplash > null"
#endif

#ifdef LenovoDebian
///the directory of FEBio, on Yu's lenovo laptop

#define RUN_FEBIO_COMMAND "/home/yuquan/febio-2.4.0-packaged/bin/febio2.lnx64 -nosplash > null"
#endif

#define NUM_THREADS 4

#define STORE_CRASH_FILE_TEXT

////#define USE_IDEALIZED_DATA

#endif // _DE_H
