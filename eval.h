//
//  eval.h
//  differential evolution/Users/yuquan/Downloads/OneDrive/Dropbox/de_terminal/source/gpc.c
//
//  Created by Yu Quan on 31/7/15.
//  Copyright (c) 2015 National Univerisity of Singapore. All rights reserved.
//

#ifndef differential_evolution_eval_h
#define differential_evolution_eval_h

#include <string>
#include "core.h"

extern DE* GLOBAL_DE_PTR;

///this class is a shell to change the parameters in the *.feb file.
class FEBio_input_config {
private:
    ///the position of the parameters in the file
    FileSizeT E_position;
    FileSizeT T_position;
    
    ///the file name
    const std::string file_name;
public:
    
    ///the default constructor
    FEBio_input_config(const std::string file_name = std::string(""));
    
private:
    ///this opens the file for reading and find the position of the paramters and store them.
    void initialize();
    
public:
    
    ///this is the function call that modifies the parameter values in the actual file.
    void operator() (floatT E, floatT T) const;
};




floatT compare( const floatTArray3D& FE_outcome);


//------objective function---------------------------------------
t_pop evaluate(int i_D, t_pop& t_tmp, long *l_nfeval, t_pop *tpa_array, int i_NP);

/**C*F****************************************************************
 **
 ** Function       :t_pop evaluate(int i_D, t_pop t_tmp, long *l_nfeval,
 **                                t_pop *tpa_array, int i_NP)
 **
 ** Parameters     :i_D         (I)    number of parameters
 **                 t_tmp       (I)    parameter vector
 **                 l_nfeval   (I/O)   counter for function evaluations
 **                 tpa_array   (I)    pointer to current population (not needed here)
 **                 i_NP   	   (I)    number of population members (not needed here)
 **
 ** Return Value   :TRUE if trial vector wins, FALSE otherwise.
 **
 ***C*F*E*************************************************************/


int left_vector_wins(const t_pop& t_trial, const t_pop& t_target);
/**C*F****************************************************************
 **
 ** Function       :int left_vector_wins(t_pop t_trial, t_pop t_target)
 **
 ** Author         :Rainer Storn
 **
 ** Description    :Selection criterion of DE. Decides when the trial
 **                 vector wins over the target vector.
 **
 ** Parameters     :t_trial    (I)   trial vector
 **                 t_target   (I)   target vector
 **
 ** Return Value   :TRUE if trial vector wins, FALSE otherwise.
 **
 ***C*F*E*************************************************************/

#endif
