//
//  engin.h
//  differential evolution
//
//  Created by Yu Quan on 1/8/15.
//  Copyright (c) 2015 National Univerisity of Singapore. All rights reserved.
//

#ifndef differential_evolution_engin_h
#define differential_evolution_engin_h

#include "defines.h"
#include "file_processing.h"
#include <iostream>
#include <random>


#define WRITE_BINARY(STREAM,X) ( STREAM.write(reinterpret_cast<const char*> (&(X)), sizeof(X)) )
#define READ_BINARY(STREAM,X) ( STREAM.read(reinterpret_cast<char*> (&(X)), sizeof(X)) )

extern char CRASH_SAVE[90000];
extern unsigned long long byteSize;

//------Typedefs---------------------------------------------------
struct mystruct
//*************************************
//** Definition of population member
//*************************************
{
    floatT fa_vector[MAXDIM];         //parameter vector
    floatT fa_cost[MAXCOST];          //vector of objectives (costs)
    floatT fa_constraint[MAXCONST];   //vector of constraints
    
};

typedef struct mystruct t_pop;

template <class ElementT>
std::ostream& print_array( std::ostream& stream, const ElementT array[],
                          size_t array_size, const char dilimiter = ' ');

std::ostream& operator<< (std::ostream& stream, const t_pop& pop);


struct output_control {
    int screen_output = 0;
    char output_filename[BUFSIZE] = "./out.out";
    char dump_filename[BUFSIZE] = "./dump.dump";
    char crash_filename[BUFSIZE] = "./crash.crash";
    int output_frequency = -1;
    int dump_frequency = -1;  ///default is no dump.
    int crash_save = 0;
};


class DE{
    friend std::ostream& operator<< (std::ostream&, const DE&);
    
private:
    ///optimization parameters
    int i_strategy;
    int b_bs_flag;
    int b_bound_constr;
    
    int i_D;
    int i_NP;
    floatT f_weight;
    floatT f_cross;
    
    floatT fa_minbound[MAXDIM];
    floatT fa_maxbound[MAXDIM];
    
    ///stopping criteria
    int i_genmax;
    floatT valueToReach;
    floatT worstBestDiff;
    floatT coeffOfVariation;
    
    ///storage of population members
    t_pop ta_pop[2*MAXPOP];
    
    ///counters
    int i_gen;
    long l_nfeval;
    
    
    ///use int to represent positions, instead of pointers.
    ///So that the positions can be easily stored into files.
    int i_pta_current, i_pta_old, i_pta_new, i_pta_swap;
    
    ///temp storages
    ///t_pop* pta_current;
    t_pop t_best;
    ///t_pop* pta_old, *pta_new, *pta_swap;
    t_pop t_bestit;
    
    ///random generator
    std::mt19937_64 generator;
    
    ///random number [0,1)
    floatT genrand();
    
    
public:
    ///objective function
    t_pop (*objective_function) (const int, t_pop&, long *const, t_pop *const, const int);
    
public:
    struct output_control text_output;
    struct output_control binary_output;
    
    
public:
    DE() { }
    
    ///constructor, copies every member
    DE(const int i_strategy, const int b_bs_flag, const int b_bound_constr,
       const int i_D, const int i_NP, const floatT f_weight,
       const floatT f_cross, const floatT _fa_minbound[], const floatT _fa_maxbound[],
       const int i_genmax, const floatT valueToReach, const floatT worstBestDiff,
       const floatT coeffOfVariation, const t_pop _ta_pop[], const int i_gen,
       const long l_nfeval, const int i_pta_current, const int i_pta_old,
       const int i_pta_new, const int i_pta_swap, const t_pop t_best, const t_pop t_bestit,
       t_pop (*const objective_function)(const int, t_pop&, long *const, t_pop *const, const int)):
    i_strategy(i_strategy), b_bs_flag(b_bs_flag), b_bound_constr(b_bound_constr),
    i_D(i_D), i_NP(i_NP), f_weight(f_weight), f_cross(f_cross),i_genmax(i_genmax),
    valueToReach(valueToReach), worstBestDiff(worstBestDiff),
    coeffOfVariation(coeffOfVariation), i_gen(i_gen), l_nfeval(l_nfeval),
    i_pta_current(i_pta_current), t_best(t_best), i_pta_old(i_pta_old), i_pta_new(i_pta_new),
    i_pta_swap(i_pta_swap), t_bestit(t_bestit),
    objective_function(objective_function)
    {
        for (int i = 0; i < MAXDIM; ++i) {
            fa_minbound[i] = _fa_minbound[i];
            fa_maxbound[i] = _fa_maxbound[i];
        }
        
        for (int i = 0; i < 2*MAXPOP; ++i) {
            ta_pop[i] = _ta_pop[i];
        }
    }
    
    ///constructor, to a initialized state. Some members are not copied; they are initialized.
    DE(const int i_strategy, const int b_bs_flag, const int b_bound_constr,
       const int i_D, const int i_NP, const floatT f_weight, const floatT f_cross,
       const floatT _fa_minbound[], const floatT _fa_maxbound[], const int i_genmax,
       const floatT valueToReach, const floatT worstBestDiff, const floatT coeffOfVariation, t_pop (*const objective_function)(const int, t_pop&, long *const, t_pop *const, const int), const long i_seed):
    i_strategy(i_strategy), b_bs_flag(b_bs_flag), b_bound_constr(b_bound_constr),
    i_D(i_D), i_NP(i_NP), f_weight(f_weight), f_cross(f_cross), i_genmax(i_genmax),
    valueToReach(valueToReach), worstBestDiff(worstBestDiff), coeffOfVariation(coeffOfVariation) {
        for (int i = 0; i < MAXDIM; ++i) {
            fa_minbound[i] = _fa_minbound[i];
            fa_maxbound[i] = _fa_maxbound[i];
        }
        i_gen = 0;l_nfeval = 0;
        generator.seed(i_seed);
    }
    
    DE(std::ifstream& stream);
    
    bool terminate(char* const to_print = nullptr) const;
    
    void optimize();
    
    floatT best_cost() { return t_best.fa_cost[0]; }
    
    
    ///clean will re-assign i_gen and l_neval to 0 so that it can be used for a fresh start of an optimization. (saves the effort to re-declare another copy of DE object which costs.
    void clean();
};

std::ostream& operator<< (std::ostream& stream, const DE& de);


///utility functions used only by core.h and core.cpp

/**C*F****************************************************************
 **
 ** Function       :void sort (t_pop ta_ary[], int i_len)
 **
 ** Author         :Rainer Storn
 **
 ** Description    :Shell-sort procedure which sorts array ta_ary[] according
 **                 to ta_ary[].fa_cost[0] in ascending order.
 **
 ** Parameters     :ta_ary[]    (I/O)   population array
 **                 i_len        (I)    length of array to be sorteds
 **
 ** Postconditions :ta_ary[] will be sorted in ascending order (according to fa_cost[0])
 **
 ***C*F*E*************************************************************/
void sort_t_pop (t_pop ta_ary[], int i_len);



void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, std::mt19937_64* generator);
/**C*F****************************************************************
 **
 ** Function       :void permute(int ia_urn2[], int i_urn2_depth)
 **
 ** Author         :Rainer Storn
 **
 ** Description    :Generates i_urn2_depth random indices ex [0, i_NP-1]
 **                 which are all distinct. This is done by using a
 **                 permutation algorithm called the "urn algorithm"
 **                 which goes back to C.L.Robinson.
 **
 ** Parameters     :ia_urn2       (O)    array containing the random indices
 **                 i_urn2_depth  (I)    number of random indices (avoided index included)
 **                 i_NP          (I)    range of indices is [0, i_NP-1]
 **                 i_avoid       (I)    is the index to avoid and is located in
 **                                      ia_urn2[0].
 **
 ** Preconditions  :# Make sure that ia_urn2[] has a length of i_urn2_depth.
 **                 # i_urn2_depth must be smaller than i_NP.
 **
 ** Postconditions :# the index to be avoided is in ia_urn2[0], so fetch the
 **                   indices from ia_urn2[i], i = 1, 2, 3, ..., i_urn2_depth.
 **
 ***C*F*E*************************************************************/


void  assigna2b(int i_D, floatT fa_a[], floatT fa_b[]);
/**C*F****************************************************************
 **
 ** Function       :void  assigna2b(int i_D, float fa_a[], float fa_b[])
 **
 ** Author         :Rainer Storn
 **
 ** Description    :Assigns i_D-dimensional vector fa_a to vector f_b.
 **
 ** Parameters     :i_D     (I)     size of vectors
 **                 fa_a[]  (I)     source vector
 **                 fa_b[]  (I)     destination vector
 ***C*F*E*************************************************************/

#endif
