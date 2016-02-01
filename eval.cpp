#include "eval.h"

#include <iomanip>
#include <iostream>

#include <boost/spirit/include/support_multi_pass.hpp>

#include "gpc.h"
#include <gsl/gsl_statistics.h>
#include "file_processing.h"


gpc_polygon* P_TRUE_DATA = nullptr;
gpc_polygon* P_FE = nullptr;


///this class is a shell to change the parameters in the *.feb file.
FEBio_input_config::FEBio_input_config(const std::string file_name): E_position(0), T_position(0), file_name( file_name )
{if(file_name.length()) initialize( );}

void FEBio_input_config::initialize() {
    std::string prior_E = "<material id=\"1\"";
    std::string E_string = "<E>";
    std::string T_string = "<loadpoint>1,";
    std::ifstream file_stream_local( file_name, std::ios::in | std::ios::binary);
    check_openfile(file_stream_local);
    
    typedef boost::spirit::multi_pass< std::istreambuf_iterator<char> > boost_it;
    boost_it begin_of_file = boost::spirit::make_default_multi_pass( std::istreambuf_iterator<char>(file_stream_local));
    boost::spirit::multi_pass< std::istreambuf_iterator<char> > end_of_file = boost::spirit::make_default_multi_pass( std::istreambuf_iterator<char>() );
    
    boost_it search_result = std::search( begin_of_file,
                                         end_of_file, prior_E.begin(), prior_E.end());
    search_result = std::search( search_result, end_of_file, E_string.begin(), E_string.end());
    E_position = std::distance( begin_of_file, search_result) + E_string.length();
    
    search_result = std::search( search_result, end_of_file, T_string.begin(), T_string.end());
    T_position = std::distance( begin_of_file, search_result) + T_string.length();
    file_stream_local.close();
    
    return;
}
///this is the function call that modifies the parameter values in the actual file.
void FEBio_input_config::operator() (floatT E, floatT T) const{
    char postfix[7];
    ///const int thread_num = omp_get_thread_num();
    const int thread_num = 0;
    std::string new_string(file_name);
    sprintf(postfix, "_%d", thread_num);
    new_string.insert(file_name.length() - 4,postfix);
    
    std::ofstream file_stream( new_string, std::ios::out | std::ios::in);
    check_openfile(file_stream);
    
    file_stream.seekp(E_position);
    file_stream << std::fixed << std::scientific << std::setprecision(PRECISION) << E;
    
    
    
    file_stream.seekp(T_position);
    file_stream << std::setprecision(PRECISION) << T;
    
    
    file_stream.close();
    return;
}


struct to_compare_x {
    bool operator() (gpc_vertex v1, gpc_vertex v2) {
        return v1.x < v2.x;
    }
};

struct to_compare_y {
    bool operator() (gpc_vertex v1, gpc_vertex v2) {
        return v1.y < v2.y;
    }
};


///this objective function is f(x,y), where [x,y] is the translation vector to the polygon P_FE, and f is the negative of the area(intersection).
t_pop evaluate2(const int i_D, t_pop& t_tmp, long *const l_nfeval, t_pop *const tpa_array, const int i_NP) {
    ++(*l_nfeval);
    
    ///translate the polygon
    gpc_translate(P_FE, t_tmp.fa_vector[0], t_tmp.fa_vector[1]);
    
    ///intersect the two polygons.
    gpc_polygon p_int;
    gpc_polygon_clip(GPC_INT, P_FE, P_TRUE_DATA, &p_int);
    
    ///translate back.
    gpc_translate(P_FE, -t_tmp.fa_vector[0], -t_tmp.fa_vector[1]);
    
    ///the cost is the negative of the area of intersection
    t_tmp.fa_cost[0] = -gpc_polygon_area(&p_int);
    
    gpc_free_polygon( &p_int );
    return t_tmp;
}


floatT compare( const floatTArray3D& FE_outcome) {
    static bool first_time = 1;///second entry onwards, this will become 0
    
    ///for pminx, pminy
    static struct to_compare_x cmp_x;
    static struct to_compare_y cmp_y;
    
    ///Matlab data and idealized data is read only once, and kept even control left the function
    static const struct Matlab_output_data Matlab_output_data = read_Matlab_output( MATLAB_OUTPUT );
    static const floatTArray3D true_data_idealized = read_FEBio_log(IDEALIZED_DATA, Matlab_output_data.n_boundary_nodes, Matlab_output_data.boundary_node_index);
    
    ///number of time steps, in FE model
    static const int num_time_step = true_data_idealized.shape()[0];  ///this is 10, for now
    
    ///this is the size of polygons for idealized data and FE_outcome. true_data it's different.
    static const int num_poly_pts = true_data_idealized.shape()[1];
    
    ///to be used to store iris contour "true data" from the OCT video
    static std::vector<floatTArray2D> true_data_left(num_time_step);
    static std::vector<floatTArray2D> true_data_right(num_time_step);
    
    
    static int hole_flag = 0;  ///used by all gpc_polygon
    
    ///stores true data as in the format required by gpc_polygon. Outer array is for 10 time steps. This structure does the storage. The next two structures are essentially addresses of the storage.
    static std::vector< std::vector<gpc_vertex> > vertices_true_data(num_time_step);
    
    ///an intermediate gpc structure. Outer array is for 10 time steps
    static std::vector<gpc_vertex_list> array_v_list_true_data(num_time_step);
    
    ///the final product of the previous two. The polygon structure for true data for 10 time steps. These never change in the entire program; hence it is static and only needed to be assigned value once.
    static std::vector<gpc_polygon> array_p_true_data(num_time_step);
    
    static floatT registration_min_bound[2] = {-0.5, -0.25};
    static floatT registration_max_bound[2] = {0.5, 0.25};
    
    static DE registration(1, 0, 1, 2, 20, 0.85, 1, registration_min_bound, registration_max_bound, 100, -999,1e-6,1e-5,  &evaluate2,0);
    
    ///The block below should be ran only when the function is executed for the first time
    if (first_time) {
        int local_num_polypts;
        read_true_data( "true_data.dat", &true_data_left[0], &true_data_right[0], num_time_step);
        
        for (int i = 0; i < num_time_step; ++i)
        {
            ///using true data (OCT video) or idealized data (from FE model). The MACRO is found in the file 'define.h'. The two blocks looks nearly the same; only the right hand side of the assignment expressions are different.
#ifdef USE_IDEALIZED_DATA
            local_num_polypts = num_poly_pts;
            vertices_true_data[i].resize(local_num_polypts);
            for (int j = 0; j < local_num_polypts; ++j) {
                vertices_true_data[i][j].x = true_data_idealized[i][j][0];
                vertices_true_data[i][j].y = true_data_idealized[i][j][2];
            }
#else
            local_num_polypts = true_data_right[i].shape()[0] - 1;
            vertices_true_data[i].resize(local_num_polypts);
            for (int j = 0; j < local_num_polypts; ++j) {
                vertices_true_data[i][j].x = (true_data_right[i])[j][0];
                vertices_true_data[i][j].y = (true_data_right[i])[j][2];
            }
#endif
            array_v_list_true_data[i] = {.num_vertices = local_num_polypts, .vertex = &vertices_true_data[i][0]};
            
            array_p_true_data[i] = {.num_contours = 1, .hole = &hole_flag, .contour = &array_v_list_true_data[i]};
            
            floatT pminx = std::min_element(array_p_true_data[i].contour -> vertex, array_p_true_data[i].contour -> vertex + array_p_true_data[i].contour->num_vertices, cmp_x) -> x;
            floatT pminy = std::min_element(array_p_true_data[i].contour -> vertex, array_p_true_data[i].contour -> vertex + array_p_true_data[i].contour->num_vertices, cmp_y) -> y;
            
                gpc_translate(&array_p_true_data[i], -pminx, -pminy);

            
        }  ///end for (i)
        
    }     ////end if (first_time)
    
    
    ///the rest are executed each time of function evaluation
    std::vector<gpc_vertex> vertices_FE(num_poly_pts);
    gpc_vertex_list v_list_FE = {.num_vertices = num_poly_pts, .vertex = &vertices_FE[0]};
    gpc_polygon p_int;
    gpc_polygon p_FE = {.num_contours = 1, .hole = &hole_flag, .contour = &v_list_FE};
    
    std::vector<floatT> cost(num_time_step);
    for (int i = 0; i < num_time_step; ++i) {
        for (int j = 0; j < num_poly_pts; ++j) {
            vertices_FE[j].x = FE_outcome[i][j][0];
            vertices_FE[j].y = FE_outcome[i][j][2];
        }
        floatT pminx = std::min_element(p_FE.contour -> vertex, p_FE.contour -> vertex + p_FE.contour->num_vertices, cmp_x) -> x;
        floatT pminy = std::min_element(p_FE.contour -> vertex, p_FE.contour -> vertex + p_FE.contour->num_vertices, cmp_y) -> y;
        
        
        gpc_translate(&p_FE, -pminx, -pminy);
        
        P_FE = &p_FE;
        P_TRUE_DATA = &array_p_true_data[i];
        ////add the image registration optimization here
        ///something like  best_cost = optimize(p1,p2)...
        ///gpc_polygon_clip(GPC_INT, &p_FE, &array_p_true_data[i], &p_int);
        
        //FILE* p1file = fopen("p1file.dat", "w");
        //FILE* p2file = fopen("p2file.dat", "w");
        
        //gpc_write_polygon(p1file, 1, P_FE);
        //gpc_write_polygon(p2file, 1, P_TRUE_DATA);
        
        //fclose(p1file); fclose(p2file);

        registration.clean(); registration.objective_function = &evaluate2;
        registration.optimize();
        cost[i] = floatT(1.0 + 2.0*registration.best_cost() / (gpc_polygon_area(&p_FE) + gpc_polygon_area(&array_p_true_data[i])));
    }
    
    first_time = 0;
    return 100.0*((floatT)gsl_stats_mean(&cost[0], 1, num_time_step));
}




//------objective function---------------------------------------
t_pop evaluate(const int i_D, t_pop& t_tmp, long *const l_nfeval, t_pop *const tpa_array, const int i_NP)

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
{
    static bool first_time = 1;
    ///iteration number and population number
    const int pop_num = (*l_nfeval)%i_NP + 1;
    const int iteration_num = (*l_nfeval)/i_NP;
    
    ///#pragma omp atomic
    (*l_nfeval)++;  //increment function evaluation count
    
    ///static variables are created only once, and are shared between threads.
    ///it is very important not to modify these variables.
    static const char original_filename[] = FEBio_FILENAME;
    static std::string _string(original_filename);
    static const FEBio_input_config change_parameters(_string);
    
    static char* command[NUM_THREADS];
    static char* log_filename[NUM_THREADS];
    
    static const struct Matlab_output_data Matlab_output_data = read_Matlab_output( MATLAB_OUTPUT );
    
    ///#pragma omp critical(first_eval)
    if (first_time) {   ///i_D == -1 means it's a restart
        for (int i = 0; i < NUM_THREADS; ++i) {
            ///initialize files for the threads.
            char copied_filename[100];
            strcpy(copied_filename, original_filename);
            insert_postfix(copied_filename, i);
            file_copy(copied_filename, original_filename);
            
            ///initialize commands to run FEBio
            command[i] = (char*)malloc(200);
            strcpy(command[i], RUN_FEBIO_COMMAND);
            strcat(command[i], " -i ");
            strcat(command[i], copied_filename);
            strcat(command[i], " \n");
            
            ///initialize file names of log files
            log_filename[i] = (char*)malloc(100);
            strcpy(log_filename[i],copied_filename);
            strcpy(&log_filename[i][strlen(log_filename[i])-3],"log");
        }
    }
    
    ///#pragma omp critical(std_cout)
    {
        std::cout << "evaluating "<< pop_num <<"th member in " << iteration_num << "th iteration...\n"<< "parameters = " << t_tmp.fa_vector[0] << "  " << t_tmp.fa_vector[1] << "\n";
    }
    const int thread_num = 0;
    ///#pragma omp critical(extra)
    ///edit the FEBio input file, just change the relevant parameters.
    change_parameters( t_tmp.fa_vector[0], t_tmp.fa_vector[1]);
    ////change_parameters( 0.0200, 0.1000);
    
    system(command[thread_num]);
    
    ///passing the time taken by FEBio in seconds
    /*if (FE_running_time) {
     time_t after; time(&after);
     *FE_running_time = difftime(after, before);
     }*/
    ///read the FEBio log file for the boundary_nodes
    
    
    ///check if log file contains "E R R O R   T E R M I N A T I O N" at the end of the file
    FileSizeT fileSize = get_filesize( log_filename[thread_num]);
    std::ifstream stream(log_filename[thread_num], std::ios::in);
    stream.seekg( fileSize - 60);
    char end_of_file[55];
    stream.read(end_of_file, 54);
    stream.close();
    char* result = nullptr;
    result = strstr(end_of_file, " E R R O R ");
    
    ///if it does contain this, then this particular member is invalid, assign it a maximal cost.
    if (result != nullptr) {
        std::cout << "5\n";
        t_tmp.fa_cost[0] = 100.0;
        return t_tmp;
    }
    
    const floatTArray3D FE_outcome = read_FEBio_log(log_filename[thread_num],Matlab_output_data.n_boundary_nodes, Matlab_output_data.boundary_node_index);
    
    ///compare the boundary_nodes, using Matlab Engine, with the true data, and produce a cost
    ///#pragma omp critical(extra)
    t_tmp.fa_cost[0] = compare(FE_outcome);
    
    ///#pragma omp critical(std_cout)
    std::cout << pop_num << "th cost = " << t_tmp.fa_cost[0] << "\n\n";
    
    first_time = 0;
    return(t_tmp);
}



int left_vector_wins(const t_pop& t_trial, const t_pop& t_target)
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
{
    ///---trial wins against target even when cost is equal.-----
    return (t_trial.fa_cost[0] <= t_target.fa_cost[0]);
}
