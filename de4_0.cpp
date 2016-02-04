/***************************************************************
 **                                                            **
 **            DE-Engine v4.0, Rainer Storn, 2004              **
 **                                                            **
 ***************************************************************/

/***************************************************************
 **                                                            **
 **            Modified by Yu Quan, 2015                       **
 **                                                            **
 ***************************************************************/

#include <iostream>

#include "eval.h"
#include "core.h"
#include "defines.h"
#include "file_processing.h"
#include "gpc.h"
#include <unistd.h>
#include <csignal>

DE* GLOBAL_DE_PTR = nullptr;
char CRASH_SAVE[90000];
unsigned long long byteSize;

void grid_evaluate(floatT min1, floatT int1, floatT max1, floatT min2, floatT int2, floatT max2);

void my_handler(int i) {
    if (GLOBAL_DE_PTR->binary_output.crash_save) {
        std::ofstream crash_file(GLOBAL_DE_PTR->binary_output.crash_filename, std::ios::out | std::ios::binary);
        crash_file.write(CRASH_SAVE, byteSize);
        crash_file.close();
        std::cerr << "termination handler has saved the status of the algorithm into file " << GLOBAL_DE_PTR->binary_output.crash_filename << '\n';
    }
    exit(1);
}


int main(int argc, char *argv[])
{
    ///the default file names. to be changed if provided by the user.
    char de_input_filename[100] = "in.dat";
    
    char restart_filename[100];
    
    int print = 1;  ///to print iteration results in terminal
    int restart = 0; ///restart from a previously incompleted optimization.
    
    ///the structures for output files. See 'core.h'
    struct output_control text_output;
    text_output.screen_output = 1;
    strcpy(text_output.output_filename, "out.dat");    ///default output file name
    text_output.output_frequency = 1;
    text_output.crash_save = 1;
    

    
    struct output_control binary_output;   ///the default values says no dumping
    
    ///parse shell command
    char option;
    while ((option = getopt(argc, argv, "i:o:p:r:b:t:f:c:g")) != -1)
        switch (option) {
            case 'i':   ///intput filename
                strcpy(de_input_filename, optarg);
                break;
           
            case 'o':   ///output filename
                strcpy(text_output.output_filename, optarg);
                break;
                
            case 'p':   ///1 to print iteration results in terminal. 0 otherwise
                text_output.screen_output = atoi(optarg);
                break;
                
            case 'r':
                restart = 1;
                strcpy(restart_filename, optarg);
                break;
            case 'b':   ///binary dump file. If not provided, no binary dump will be performed.
                strcpy(binary_output.dump_filename, optarg);
                binary_output.dump_frequency = 1;
                break;
                
            case 't':   ///text dump file. If not provided, no text dump will be performed.
                strcpy(text_output.dump_filename, optarg);
                text_output.dump_frequency = 1;
                break;
                
            case 'f':   ///dump frequency. Applies to both binary and text dump.
                text_output.dump_frequency = atoi(optarg);
                binary_output.dump_frequency = atoi(optarg);
                break;
                
            case 'c':   ///crash filename. if not set the default will be used.
                strcpy(binary_output.crash_filename, optarg);
                break;
                
            case 'g':
                grid_evaluate(0.2,0.2/60.0,0.5, 0.1,20.0/60.0,14);
                return 0;
                break;
                
                
            case '?':
                if (optopt == 'i' || optopt == 'o' || optopt == 'p' || optopt == 'b'|| optopt == 't' || optopt == 'f' || optopt == 'c')
                    std::cerr << "Option -"<<char(optopt)<<" requires an argument.\n";
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                return 1;
                
            default:
                abort();
        }
    std::ifstream in(de_input_filename, std::ios::in);
    
    DE de(in);
    
    if (restart) {
        std::ifstream restart_stream(restart_filename, std::ios::in | std::ios::binary);
        READ_BINARY(restart_stream, de);
    }
    
    GLOBAL_DE_PTR = &de;
    
    de.text_output = text_output;
    de.binary_output = binary_output;
    
    ///std::set_terminate(myterminate);
    signal(SIGABRT, my_handler);
    signal(SIGFPE, my_handler);
    signal(SIGILL, my_handler);
    signal(SIGINT, my_handler);
    signal(SIGSEGV, my_handler);
    signal(SIGTERM, my_handler);
    
    de.optimize();
    
    return 0;
}


void grid_evaluate(floatT min1, floatT int1, floatT max1, floatT min2, floatT int2, floatT max2) {
    ///below is for evaluating objective function at specified grid
    
    std::cout << "1\n";
    t_pop t;
    floatT x1, x2;
    std::ofstream file("grid_evaluate2.dat", std::ios::out);
    long i = 0;
    
    for (x1 = min1; x1 <= max1; x1 += int1) {
        for (x2 = min2; x2 <= max2; x2 += int2) {
            t.fa_vector[0] = x1; t.fa_vector[1] = x2;
            evaluate(2,t,&i,&t,20);
            file << t.fa_vector[0] << "   " << t.fa_vector[1] << "   " << t.fa_cost[0] << '\n' << std::flush;
        }
    }
    
    file.close();
}


