#include "core.h"

#include "file_processing.h"
#include "eval.h"
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <gsl/gsl_statistics.h>


floatT DE::genrand() {
    floatT x;
    ///#pragma omp critical(rand)
    x = (floatT)generator()/(floatT)generator.max();
    return x;
}



DE::DE(std::ifstream& stream) {
    char* line;
    
    line = line_wo_comment(stream); i_strategy = atoi(line);
    line = line_wo_comment(stream); b_bs_flag = atoi(line);
    line = line_wo_comment(stream); b_bound_constr = atoi(line);
    line = line_wo_comment(stream); i_genmax = atoi(line);
    line = line_wo_comment(stream); /*i_refresh =*/ atoi(line);
    line = line_wo_comment(stream); i_D = atoi(line);
    line = line_wo_comment(stream); i_NP = atoi(line);
    line = line_wo_comment(stream); f_weight = atof(line);
    line = line_wo_comment(stream); f_cross = atof(line);
    line = line_wo_comment(stream); valueToReach = atof(line);
    line = line_wo_comment(stream); worstBestDiff = atof(line);
    line = line_wo_comment(stream); coeffOfVariation = atof(line);
    
    line = line_wo_comment(stream);

    long long i_seed;
    if (memcmp(line, "time", 4) == 0)
        i_seed = std::chrono::system_clock::now().time_since_epoch().count();
    else
        i_seed = atol(line);
    generator.seed(i_seed);
    
    for (int i = 0; i < i_D; ++i) {
        line = line_wo_comment(stream); fa_minbound[i] = atof(line);
    }
    
    for (int i = 0; i < i_D; ++i) {
        line = line_wo_comment(stream); fa_maxbound[i] = atof(line);
    }
    
    i_gen = 0;
    l_nfeval = 0;
    
    objective_function = &evaluate;
    
    free(line);
}

void DE::optimize() {
    std::ofstream output_text, dump_text, dump_binary;
    
    if (text_output.output_frequency >= 1)
        output_text.open(text_output.output_filename, std::ios::out);
    
    if (text_output.dump_frequency >= 1)
        dump_text.open(text_output.dump_filename, std::ios::out);
    
    output_text << "iteration   best cost  parameters   worst cost   CoeffOfVariation" << '\n' << std::flush;

    
    
#define URN_DEPTH   5   //4 + one index to avoid
    int   i, j, k;
    int   i_r1, i_r2, i_r3, i_r4;
    const int i_refresh = 1;
    int   ia_urn2[URN_DEPTH];
    floatT f_jitter, f_dither;
    t_pop t_tmp, t_origin;
    int restarting_loop = 0; ///for restart algorithm
    
    char buffer[99];
    
    ///if initialization has not completed, redo from the start.
    if (i_gen == 0) {
        l_nfeval = 0;
    }
    
    if (l_nfeval == 0) {
        assert(i_gen == 0);
        
        i_pta_current = 0;
        
        //------Initialization-----------------------------
        for (j=0; j<i_D; j++)
        {
            ta_pop[0].fa_vector[j] = fa_minbound[j]+genrand()*(fa_maxbound[j] - fa_minbound[j]);
        }
        ta_pop[0]      = objective_function(i_D,ta_pop[0],&l_nfeval,&ta_pop[0],i_NP);

        t_best  = ta_pop[0];
        
        for (i=1; i<i_NP; i++)
        {
            for (j=0; j<i_D; j++)
            {
                ta_pop[i].fa_vector[j] = fa_minbound[j]+genrand()*(fa_maxbound[j] - fa_minbound[j]);
                
            }
        }
        
        for (i = 1; i < i_NP; ++i) {
            
            ta_pop[i] = objective_function(i_D,ta_pop[i],&l_nfeval,&ta_pop[0],i_NP);
            
            if (left_vector_wins(ta_pop[i],t_best) == true)
            {
                t_best = ta_pop[i];
            }
        }
        
        ///output part of initialization iteration
        if (text_output.output_frequency == 1)
        {
            output_text << std::setw(6) << i_gen << std::setprecision(6) << std::setw(12) <<
            t_best.fa_cost[0] << std::setw(12) << t_best.fa_vector[0] << "  " << std::setw(12) << t_best.fa_vector[1] << '\n' << std::flush;
        }
        
        if (text_output.screen_output)
        {
            std::cout << "iteration " << std::setw(6) << i_gen <<"  best cost  " << std::setprecision(6) << std::setw(12) <<
            t_best.fa_cost[0] << "  parameters  " << std::setw(12) << t_best.fa_vector[0] << "  " << std::setw(12) << t_best.fa_vector[1] << '\n';
        }
        
        if (text_output.dump_frequency == 1)
        {
            dump_text << "after iteration " << i_gen << '\n' << (*this) << std::flush;
        }
        
        if (binary_output.dump_frequency == 1)
        {
            dump_binary.open(binary_output.dump_filename, std::ios::binary);
            WRITE_BINARY(dump_binary, *this);
            dump_binary.close();
        }
        
        i_gen++; ///now i_gen == 1
        ///end of initialization iteration
        
        t_bestit  = t_best;
        
        //---assign pointers to current ("old") and new population---
        
        i_pta_old = 0;
        i_pta_new = &ta_pop[i_NP] - ta_pop;
        
    }  ///if (l_nfeval == 0)
    
    else ///means restart from last time
    {
        ++i_gen;
        objective_function = &evaluate;
        restarting_loop = 1;
    }
    
    //------Iteration loop--------------------------------------------
    
    //Note that kbhit() needs conio.h which is not always available under Unix.
    do
    {
        //----computer dithering factor (if needed)-----------------
        f_dither = f_weight + genrand()*(1.0 - f_weight);

        //----start of loop through ensemble------------------------
        for (i=0; i<i_NP; i++)
        {

            ///do for every population member
            permute(ia_urn2,URN_DEPTH,i_NP,i, &generator); //Pick 4 random and distinct

            i_r1 = ia_urn2[1];                 //population members
            i_r2 = ia_urn2[2];
            i_r3 = ia_urn2[3];
            i_r4 = ia_urn2[4];
            
            //========Choice of strategy=======================================================
            //---classical strategy DE/rand/1/bin-----------------------------------------
            
            
            ///computes based on random generator, the parameters, for all populations. This differentiates the different strategies.
            
            if (i_strategy == 1)
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                
                do
                {                            // add fluctuation to random target
                    t_tmp.fa_vector[j] = ta_pop[i_pta_old +i_r1].fa_vector[j] + f_weight*(ta_pop[i_pta_old +i_r2].fa_vector[j]-ta_pop[i_pta_old +i_r3].fa_vector[j]);
                    
                    j = (j+1)%i_D;
                    k++;
                }while((genrand() < f_cross) && (k < i_D));
                
                if (b_bound_constr)
                    assigna2b(i_D,ta_pop[i_pta_old +i_r1].fa_vector,t_origin.fa_vector);
            }
            //---DE/local-to-best/1/bin---------------------------------------------------
            else if (i_strategy == 2)
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                do
                {                            // add fluctuation to random target
                    t_tmp.fa_vector[j] = t_tmp.fa_vector[j] + f_weight*(t_bestit.fa_vector[j] - t_tmp.fa_vector[j]) +
                    f_weight*(ta_pop[i_pta_old +i_r2].fa_vector[j]-ta_pop[i_pta_old +i_r3].fa_vector[j]);
                    
                    j = (j+1)%i_D;
                    k++;
                }while((genrand() < f_cross) && (k < i_D));
                
                if (b_bound_constr)
                    assigna2b(i_D,t_tmp.fa_vector,t_origin.fa_vector);
            }
            //---DE/best/1/bin with jitter------------------------------------------------
            else if (i_strategy == 3)
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                do
                {                            // add fluctuation to random target
                    f_jitter = (0.0001*genrand()+f_weight);
                    t_tmp.fa_vector[j] = t_bestit.fa_vector[j] + f_jitter*(ta_pop[i_pta_old +i_r1].fa_vector[j]-ta_pop[i_pta_old +i_r2].fa_vector[j]);
                    
                    j = (j+1)%i_D;
                    k++;
                }while((genrand() < f_cross) && (k < i_D));
                if (b_bound_constr)
                    assigna2b(i_D,t_tmp.fa_vector,t_origin.fa_vector);
            }
            //---DE/rand/1/bin with per-vector-dither-------------------------------------
            else if (i_strategy == 4)
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                do
                {                            // add fluctuation to random target
                    t_tmp.fa_vector[j] = ta_pop[i_pta_old +i_r1].fa_vector[j] +
                    (f_weight + genrand()*(1.0 - f_weight))*
                    (ta_pop[i_pta_old +i_r2].fa_vector[j]-ta_pop[i_pta_old +i_r3].fa_vector[j]);
                    
                    j = (j+1)%i_D;
                    k++;
                }while((genrand() < f_cross) && (k < i_D));
                if (b_bound_constr)
                    assigna2b(i_D,t_tmp.fa_vector,t_origin.fa_vector);
            }
            //---DE/rand/1/bin with per-generation-dither---------------------------------
            else if (i_strategy == 5)
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                do
                {                            // add fluctuation to random target
                    t_tmp.fa_vector[j] = ta_pop[i_pta_old +i_r1].fa_vector[j] + f_dither*(ta_pop[i_pta_old +i_r2].fa_vector[j]-ta_pop[i_pta_old +i_r3].fa_vector[j]);
                    
                    j = (j+1)%i_D;
                    k++;
                }while((genrand() < f_cross) && (k < i_D));
                if (b_bound_constr)
                    assigna2b(i_D,t_tmp.fa_vector,t_origin.fa_vector);
            }
            //---variation to DE/rand/1/bin: either-or-algorithm--------------------------
            else
            {
                assigna2b(i_D,ta_pop[i_pta_old +i].fa_vector,t_tmp.fa_vector);
                j = (int)(genrand()*i_D); // random parameter
                k = 0;
                if (genrand() < 0.5) //Pmu = 0.5
                {//differential mutation
                    do
                    {                            // add fluctuation to random target
                        t_tmp.fa_vector[j] = ta_pop[i_pta_old +i_r1].fa_vector[j] + f_weight*(ta_pop[i_pta_old +i_r2].fa_vector[j]-ta_pop[i_pta_old +i_r3].fa_vector[j]);
                        
                        j = (j+1)%i_D;
                        k++;
                    }while((genrand() < f_cross) && (k < i_D));
                }
                else
                {//recombination with K = 0.5*(F+1) --> F-K-Rule
                    do
                    {                            // add fluctuation to random target
                        t_tmp.fa_vector[j] = ta_pop[i_pta_old +i_r1].fa_vector[j] + 0.5*(f_weight+1.0)*
                        (ta_pop[i_pta_old +i_r2].fa_vector[j]+ta_pop[i_pta_old +i_r3].fa_vector[j] -
                         2*ta_pop[i_pta_old +i_r1].fa_vector[j]);
                        
                        j = (j+1)%i_D;
                        k++;
                    }while((genrand() < f_cross) && (k < i_D));
                }
                if (b_bound_constr)
                    assigna2b(i_D,ta_pop[i_pta_old +i_r1].fa_vector,t_origin.fa_vector);
            }//end if (gi_strategy ...
            
            
            ///if parameter is out of bound, then do this:
            if (b_bound_constr)
                for (j=0; j<i_D; j++) //----boundary constraints via random reinitialization-------
                {                      //----and bounce back----------------------------------------
                    if (t_tmp.fa_vector[j] < fa_minbound[j])
                    {
                        t_tmp.fa_vector[j] = fa_minbound[j]+genrand()*(t_origin.fa_vector[j] - fa_minbound[j]);
                    }
                    if (t_tmp.fa_vector[j] > fa_maxbound[j])
                    {
                        t_tmp.fa_vector[j] = fa_maxbound[j]+genrand()*(t_origin.fa_vector[j] - fa_maxbound[j]);
                    }
                }

            ///computes the cost function
            //------Trial mutation now in t_tmp-----------------
            t_tmp = objective_function(i_D,t_tmp,&l_nfeval,&ta_pop[i_pta_old +0],i_NP);  // Evaluate mutant in t_tmp[]
            restarting_loop = 0;
            
            if (b_bs_flag)
            {
                ta_pop[i_pta_new +i]=t_tmp; //save new vector, selection will come later
            }
            else
            {
                if (left_vector_wins(t_tmp,ta_pop[i_pta_old +i]) == true)
                {
                    ta_pop[i_pta_new +i]=t_tmp;              // replace target with mutant
                    
                    if (left_vector_wins(t_tmp,t_best) == true)// Was this a new minimum?
                    {                               // if so...
                        t_best = t_tmp;             // store best member so far
                    }                               // If mutant fails the test...
                }                                  // go to next the configuration
                else
                {
                    ta_pop[i_pta_new +i]=ta_pop[i_pta_old +i];              // replace target with old value
                }
            }//if (b_bs_flag)
        }
        // End mutation loop through pop.
        
        if (b_bs_flag)
        {
            sort_t_pop (ta_pop+i_pta_old, 2*i_NP); //sort array of parents + children
            t_best = ta_pop[i_pta_old +0];
        }
        else
        {
            i_pta_swap = i_pta_old;
            i_pta_old  = i_pta_new;
            i_pta_new  = i_pta_swap;
        }//if (b_bs_flag)
        
        i_pta_current = i_pta_old;
        t_bestit = t_best;
        
        
        //======Output Part=====================================================
        
        if ((text_output.output_frequency >= 1) && (i_gen%text_output.output_frequency == 0))
        {
            output_text << std::setw(6) << i_gen << std::setprecision(6) << std::setw(12) <<
            t_best.fa_cost[0] << std::setw(12) << t_best.fa_vector[0] << "  " << std::setw(12) << t_best.fa_vector[1];
        }
        
        if (text_output.screen_output)
        {
            std::cout << "iteration " << std::setw(6) << i_gen <<"  best cost  " << std::setprecision(6) << std::setw(12) <<
            t_best.fa_cost[0] << "  parameters  " << std::setw(12) << t_best.fa_vector[0] << "  " << std::setw(12) << t_best.fa_vector[1];
        }
        
        if ((text_output.dump_frequency >= 1) && (i_gen%text_output.dump_frequency == 0))
        {
            dump_text << "after iteration " << i_gen << '\n' << (*this) << std::flush;
        }
        
        if ((binary_output.dump_frequency >= 1) && (i_gen%binary_output.dump_frequency == 0))
        {
            dump_binary.open(binary_output.dump_filename, std::ios::binary);
            WRITE_BINARY(dump_binary, *this);
            dump_binary.close();
        }
        
        if (text_output.crash_save) {
            //std::ofstream save_stream("crash.crash", std::ios::out | std::ios::binary);
            //WRITE_BINARY(save_stream, *this);
            //save_stream.close();
            byteSize = sizeof(*this);
            memcpy(CRASH_SAVE, reinterpret_cast<const char*> (this), byteSize);
        }
        
        ++i_gen;
        
    } while (!terminate(buffer) && ( ((text_output.output_frequency >= 1) && ((i_gen-1)%text_output.output_frequency == 0))? bool(output_text<<buffer<<std::flush):1));
}





std::ostream& operator<< (std::ostream& stream, const DE& de) {
    ///follow the sequence of that in the class header.
    stream << de.i_strategy << '\n';
    stream << de.b_bs_flag << '\n';
    stream << de.b_bound_constr << '\n';
    
    stream << de.i_D << '\n';
    stream << de.i_NP << '\n';
    stream << de.f_weight << '\n';
    stream << de.f_cross << '\n';
    
    
    print_array(stream, de.fa_minbound, MAXDIM); stream << '\n';
    print_array(stream, de.fa_maxbound, MAXDIM); stream << '\n';
    
    stream << de.i_genmax << '\n';
    stream << de.valueToReach << '\n';
    stream << de.worstBestDiff << '\n';
    stream << de.coeffOfVariation << '\n';
    
    ///print_array(stream, de.ta_pop, 2*MAXPOP, '\n'); stream << '\n';
    ///the line above is not working. unknown reason
    
    for (int i = 0; i < de.i_NP; ++i)
        stream  << i << '\n'  << de.ta_pop[i] << '\n';
    
    stream << de.i_gen << '\n';
    stream << de.l_nfeval << '\n';
    stream << de.i_pta_current << '\n';
    stream << de.i_pta_old << '\n';
    stream << de.i_pta_new << '\n';
    stream << de.i_pta_swap << '\n';
    
    stream << de.t_best << '\n';
    stream << de.t_bestit << '\n';
    
    stream << de.generator << '\n';
    ///c++11 provides this function. The random generator class can be output and retrived precisely, in terms of the state of the generator.
    
    return stream;
}




bool DE::terminate(char* const to_print) const{
    int (*function_pt) (const t_pop&, const t_pop&) = &left_vector_wins;
    
    floatT min_value = std::min_element(ta_pop + i_pta_old, ta_pop + i_pta_old + i_NP, function_pt)->fa_cost[0];
    
    floatT max_value = std::max_element(ta_pop + i_pta_old, ta_pop + i_pta_old + i_NP, function_pt)->fa_cost[0];
    
    ///computes coefficient of variation, using the gsl library
    floatT arr_f_cost[i_NP];
    for (int i = 0; i < i_NP; ++i)
        arr_f_cost[i] = ta_pop[i_pta_old+i].fa_cost[0];
    floatT current_coeffOfVariationans = (gsl_stats_sd( (double*)arr_f_cost, 1, i_NP) / std::abs(gsl_stats_mean( (double*)arr_f_cost, 1, i_NP)));
    
    if (text_output.screen_output)
        std::cout << "   max value  " << max_value <<"  coeffOfVariation  " << current_coeffOfVariationans<< '\n';
    
    if (to_print != nullptr)
        sprintf(to_print, "   %17.10f %17.10f\n", max_value, current_coeffOfVariationans);

    return (min_value < valueToReach || (max_value - min_value) < worstBestDiff || current_coeffOfVariationans < coeffOfVariation);
}

void DE::clean() {
    i_gen = 0; l_nfeval = 0;
}


template <class ElementT>
std::ostream& print_array( std::ostream& stream, const ElementT array[],
                          size_t array_size, const char dilimiter) {
    for (int i = 0; i < array_size; ++i) {
        stream << array[i] << dilimiter;
    }
    stream << '\n';
    return stream;
}

std::ostream& operator<< (std::ostream& stream, const t_pop& pop) {
    print_array(stream, pop.fa_vector, MAXDIM);
    print_array(stream, pop.fa_cost, MAXCOST);
    print_array(stream, pop.fa_constraint, MAXCONST);
    return stream;
}


void sort_t_pop (t_pop ta_ary[], int i_len)
{
    int   done;
    int   step, bound, i, j;
    t_pop temp;
    
    step = i_len;  //array length
    while (step > 1)
    {
        step /= 2;	//halve the step size
        do
        {
            done   = true;
            bound  = i_len - step;
            for (j = 0; j < bound; j++)
            {
                i = j + step + 1;
                if (ta_ary[j].fa_cost[0] > ta_ary[i-1].fa_cost[0])
                {
                    temp     = ta_ary[i-1];
                    ta_ary[i-1] = ta_ary[j];
                    ta_ary[j]   = temp;
                    done = false; //if a swap has been made we are not finished yet
                }  // if
            }  // for
        } while (done == false);   // while
    } //while (step > 1)
} //end of sort()



void permute(int ia_urn2[], int i_urn2_depth, int i_NP, int i_avoid, std::mt19937_64* generator)
{
    /*
     static long long seed = std::chrono::system_clock::now().time_since_epoch().count();
     srand(seed);
     static std::mt19937_64 generator( rand() );*/
    
    int  i, k, i_urn1, i_urn2;
    int  ia_urn1[MAXPOP] = {0};      //urn holding all indices
    
    k      = i_NP;
    i_urn1 = 0;
    i_urn2 = 0;
    for (i=0; i<i_NP; i++) ia_urn1[i] = i; //initialize urn1
    
    i_urn1 = i_avoid;                  //get rid of the index to be avoided and place it in position 0.
    while (k >= i_NP-i_urn2_depth)     //i_urn2_depth is the amount of indices wanted (must be <= NP)
    {
        ia_urn2[i_urn2] = ia_urn1[i_urn1];      //move it into urn2
        ia_urn1[i_urn1] = ia_urn1[k-1]; //move highest index to fill gap
        k = k-1;                        //reduce number of accessible indices
        i_urn2 = i_urn2 + 1;            //next position in urn2
        floatT U01 = (floatT((*generator)()) / floatT((*generator).max()));
        i_urn1 = (int)(U01*k);    //choose a random index
    }
}


void  assigna2b(int i_D, floatT fa_a[], floatT fa_b[])
{
    int j;
    for (j=0; j<i_D; j++)
    {
        fa_b[j] = fa_a[j];
    }
}