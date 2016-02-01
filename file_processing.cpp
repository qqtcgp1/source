//
//  file_processing.cpp
//  
//
//  Created by Yu Quan on 14/9/15.
//
//

#include "file_processing.h"
#include <iostream>

FileSizeT get_filesize(const char *const filename) {
    std::ifstream file_stream(filename,  std::ios::ate | std::ios::in | std::ios::binary);
    check_openfile(file_stream);
    FileSizeT answer = file_stream.tellg();
    file_stream.close();
    return answer;
}


char* get_filestring(const char* const filename) {
    const FileSizeT filesize = get_filesize(filename);
    char *const filestring = new char[ (FileSizeT)ceil(1.2*filesize) ];
    ///char * filestring = (char*) malloc((FileSizeT)ceil(1.2*filesize) * sizeof(char));
    FILE *file = fopen( filename, "rb");
    if (!file) {std::cerr << "File could not be opened\n"; exit(EXIT_FAILURE);}
    fread(filestring, 1, filesize, file);
    fclose(file);
    return filestring;
}

char* get_filestring(const char *const filename, FileSizeT * file_size) {
    const FileSizeT filesize = get_filesize(filename);
    char *const filestring = new char[ (FileSizeT)ceil(1.2*filesize) ];
    ///char * filestring = (char*)malloc((FileSizeT)ceil(1.2*filesize));
    FILE *file = fopen( filename, "rb");
    if (!file) {std::cerr << "File could not be opened\n"; exit(EXIT_FAILURE);}
    fread(filestring, 1, filesize, file);
    fclose(file);
    *file_size = filesize;
    return filestring;
}

char* begin_indicator(const char *const indicator_str) {
    int string_lengh = strlen(indicator_str);
    char *const result = new char[ string_lengh + 5];
    result[0] = '<';
    strcpy( result+1, indicator_str);
    result[ string_lengh + 1 ] = '>';
    result[ string_lengh + 2 ] = '\0';
    return result;
}

char* end_indicator(const char *const indicator_str) {
    int string_lengh = strlen(indicator_str);
    char *const result = new char[ string_lengh + 6];
    result[0] = '<'; result[1] = '/';
    strcpy( result+2, indicator_str);
    result[ string_lengh + 2 ] = '>'; result[ string_lengh + 3 ] = '\0';
    return result;
}

const char* into_section(const char* main_string, const char *const indicator_str) {
    char* indicator = begin_indicator(indicator_str);
    const char* result = strstr(main_string, indicator) + strlen(indicator);
    delete [] indicator;
    return result;
}

const char* end_section(const char* main_string, const char *const indicator_str) {
    char* indicator = end_indicator(indicator_str);
    const char* result = strstr(main_string, indicator) + strlen(indicator);
    delete [] indicator;
    return result;
}

floatT read_attribute( const char *& string, const char * attribute ) {
    string = strstr(string, attribute) + strlen( attribute);
    char* c_pt;
    floatT d = strtod(string, &c_pt);
    string = c_pt;
    return d;
}

char* line_wo_comment(std::istream& stream) {
    char* line = (char*)malloc(1000);
    char* hex = nullptr;
    while (stream.getline(line,999))
    {
        if (line[0] != '#' && (!isnewline(line[0])) && (line[0] != '\0')) {
            ///check and erase end-of-line comments
            if ((hex = (char*)memchr(line, '#', strlen(line))) != nullptr)
                (*hex) = '\0'; ///this discards anything  after and including '#'.
            return line;
        }
        else continue;
    }
    free(line);
    return nullptr;
}


void set_extension(char* const filename, const char* const extension) {
    if (filename[strlen(filename) - strlen(extension) - 1] == '.')
        return;
    else strcat(strcat(filename, "."),extension);
}


void insert_postfix(char file_name[], int postfix) {
    int length = strlen(file_name);
    char postfix_char[5];
    sprintf(postfix_char, "_%d", postfix);
    char extension[5];
    strcpy(extension, file_name + length - 4); ///here extension includes the dot in front.
    strcpy(file_name + length - 4, postfix_char);
    strcat(file_name, extension);
}

void file_copy(const char* const destintation_file, const char* const source_file) {
    FileSizeT source_file_size;
    char* source_file_string = get_filestring(source_file, &source_file_size);
    FILE* file;
    file = fopen(destintation_file, "w");
    fwrite(source_file_string, 1, source_file_size, file);
    fclose(file);
}

void check_openfile( std::ofstream& stream ) {
    if (stream.fail()) {
        std::cerr << "File could not be opened.\n";
        exit(EXIT_FAILURE);
    }
}

void check_openfile( std::ifstream& stream ) {
    if (stream.fail()) {
        std::cerr << "File could not be opened.\n";
        exit(EXIT_FAILURE);
    }
}


template <class Iterator_type>
bool skip_line( Iterator_type* string, const int num_lines,
               const char dilimiter, const int max_line_length) {
    int j;
    for (int i = 0; i < num_lines; ++i) {
        for (j = 0; j < max_line_length; ++j, ++(*string)) {
            if ((**string) == dilimiter) {++(*string);break;}
            if ((**string) == '\0') return EXIT_FAILURE;
        }
        if (j==max_line_length) return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


void read_true_data(const char* const file_name , floatTArray2D* true_data_left, floatTArray2D* true_data_right, const int num_steps) {
    
    double a,b;
    char line[100];
    char* endpt;
    char* endpt2;
    int step_num, num_points;
    FILE* file = fopen(file_name, "r");
    
    if (file == nullptr) {
        std::cerr << "file could not be opened\n";
        exit(1);
    }
    
 
    for (int i = 0; i < num_steps; ++i) {
        fscanf(file, " step %d left, %d points\n", &step_num, &num_points);
        assert(step_num == i+1);
        true_data_left[i].resize(boost::extents[num_points][3]);
        for (int j = 0; j < num_points; ++j) {
            ///fscanf(file, " %f %f \n", &((true_data_left[i])[j][0]), &((true_data_left[i])[j][2]));
            ///fscanf(file, "  %f %f", &a, &b);   ////f*ck don't know why fscanf doesn't work!!!!!
            fgets(line, 99, file);
            a = strtod(line, &endpt);
            b = strtod(endpt, &endpt2);
            (true_data_left[i])[j][0] = a;
            (true_data_left[i])[j][2] = b;
            ////std::cout << a << ' ' << b << ' ' << feof(file) << ' ' << ferror(file) << '\n';
        }
        
        fscanf(file, " step %d right, %d points\n", &step_num, &num_points);
        assert(step_num == i+1);
        true_data_right[i].resize(boost::extents[num_points][3]);
        for (int j = 0; j < num_points; ++j) {
            ///fscanf(file, " %f %f \n", &((true_data_right[i])[j][0]), &((true_data_right[i])[j][2]));
            ///fscanf(file, " %f %f", &a, &b);
            fgets(line, 99, file);
            a = strtod(line, &endpt);
            b = strtod(endpt, &endpt2);
            (true_data_right[i])[j][0] = a;
            (true_data_right[i])[j][2] = b;
            ////std::cout << a << ' ' << b << ' ' << feof(file) << ' ' << ferror(file) << '\n';
        }
    }
    
    fclose(file);
}




const struct Matlab_output_data read_Matlab_output( const char* const file_name ) {
    ///read a file into a unmodifiable string matlab_output_file_manipulate
    const char * matlab_output_file_manipulate = get_filestring( file_name );
    
    const int n_boundary_nodes = read_attribute( matlab_output_file_manipulate, "Number of boundary nodes:");
    const int n_nodes = read_attribute( matlab_output_file_manipulate, "Number of nodes:");
    const int n_elements = read_attribute( matlab_output_file_manipulate, "Number of elements:");
    
    std::vector<int> boundary_node_index(n_boundary_nodes);
    floatTArray2D initial_nodes(boost::extents[n_nodes][3]);
    
    matlab_output_file_manipulate = into_section(matlab_output_file_manipulate, "boundary node set");
    
    
    ///the strtof function reads in. the p_end ptr is modified to point to the white space after each number. however when the next text expression is not a valid numeric expression, the p_end ptr is unmodified by strtof, and this is when the loop stops.
    std::vector<int>::iterator int_vector_iterator; char* temp_char_pt;
    for (int_vector_iterator = boundary_node_index.begin(); int_vector_iterator != boundary_node_index.end(); ++int_vector_iterator) {
        *int_vector_iterator = strtol( matlab_output_file_manipulate, &temp_char_pt, 10) ;
        matlab_output_file_manipulate = temp_char_pt;}
    
    ///reads into the <Nodes> section
    matlab_output_file_manipulate = end_section(matlab_output_file_manipulate, "boundary node set");
    matlab_output_file_manipulate = into_section(matlab_output_file_manipulate, "Geometry");
    matlab_output_file_manipulate = into_section(matlab_output_file_manipulate, "Nodes");
    
    
    ///reads the coordinates into this vector inital_nodes
    ///for each line, search using memchr for character '"' followed by '>'. if '"' is not within 20 characters of the search, then the end </Nodes> has been reached, stop the loop (in which case temp_const_char_pt == nullptr due to the behavior of memchr function.
    for (floatTArray2D::iterator iterator_i = initial_nodes.begin(); iterator_i != initial_nodes.end(); ++iterator_i) {
        matlab_output_file_manipulate = (const char*)memchr(matlab_output_file_manipulate, '"', 20);
        matlab_output_file_manipulate=(const char*)memchr(matlab_output_file_manipulate, '>', 10);
        ++matlab_output_file_manipulate; ///move on to the next character of '>', which is usually a whitespace
        
        ///scans in the three numeric items
        *iterator_i->begin() = strtod( matlab_output_file_manipulate, &temp_char_pt);
        matlab_output_file_manipulate = temp_char_pt + 1;///p_end points to the comma after the number. we move to the next character after the comma
        *(iterator_i->begin()+1) = strtod( matlab_output_file_manipulate, &temp_char_pt);
        matlab_output_file_manipulate = temp_char_pt + 1;
        *(iterator_i->begin()+2) = strtod( matlab_output_file_manipulate, &temp_char_pt);
        
        ///search for the newline character to prepare for the reading of next line
        matlab_output_file_manipulate=(const char*)memchr(matlab_output_file_manipulate, '\n', 40);
    }
    
    
    floatTArray2D initial_boundary_nodes(boost::extents[n_boundary_nodes][3]);
    floatTArray2D::iterator temp_floatT_it = initial_boundary_nodes.begin();
    for (int i = 0; i < n_boundary_nodes; ++i, ++temp_floatT_it)
        *temp_floatT_it = initial_nodes[ boundary_node_index[i]-1 ];
        
        ///Matlab_output_data to_return = { .n_nodes = n_nodes, .n_elements = n_elements,
        /// .n_boundary_nodes = n_boundary_nodes, .boundary_node_index = boundary_node_index,
        /// .initial_nodes = initial_nodes, .initial_boundary_nodes = initial_boundary_nodes };
        
        return Matlab_output_data { n_nodes, n_boundary_nodes,
n_elements, boundary_node_index, initial_nodes, initial_boundary_nodes };
        }
    
const floatTArray3D read_FEBio_log(const char* const file_name , const int n_boundary_nodes, const std::vector<int> boundary_node_index) {

    FileSizeT filesize = 0;
    
    ///read a file into a unmodifiable string matlab_output_file_start
    const char *const FEBio_logfile_start = get_filestring( file_name, &filesize );
    const char * FEBio_logfile_manipulate = FEBio_logfile_start;
    
    ///scan n_steps. search from end to beginning
    std::reverse_iterator< const char* > reverse_begin( FEBio_logfile_manipulate + (filesize-1)), reverse_end(FEBio_logfile_manipulate);

    ///reverse string of "Number of time steps completed"
    std::string temp_stdstring = "= petS";
    
    std::reverse_iterator< const char* > search_result = std::search(reverse_begin, reverse_end, temp_stdstring.begin(), temp_stdstring.end());
    
    const char * temp_const_char_pt = (search_result).base();
    ///temp_const_char_pt = (const char*)memchr(temp_const_char_pt, ':', 100) + 1;
    
    char* temp_char_pt;
    const int n_steps = (int)strtol( temp_const_char_pt, &temp_char_pt, 10);
    /////this is the number of simulation steps!!!!!!!!!!! And is also the size1 of the main data array, boundary_nodes
    ///n_steps should be 10, for now
    
    floatTArray3D boundary_nodes(boost::extents[n_steps][n_boundary_nodes][3]);
    ///boundary_nodes.resize(boost::extents[n_steps + 1][n_boundary_nodes][3]);
    
    char temp_char_arr[20];
    
    floatTArray3D::iterator iterator_loop1= boundary_nodes.begin();
    floatTArray3D::subarray<2>::type::iterator iterator_loop2 =iterator_loop1->begin();

    for (int i = 1, diff, last_index;
         iterator_loop1 != boundary_nodes.end(); ++i, ++iterator_loop1) {
        
        ///move the control to step i.
        sprintf(temp_char_arr, "Step = %d", i);
        if (!(FEBio_logfile_manipulate = strstr( FEBio_logfile_manipulate, temp_char_arr ))) break;
        assert( !skip_line(&FEBio_logfile_manipulate, 2) );
        
        last_index = 0;
        iterator_loop2 = iterator_loop1->begin();
        for (std::vector<int>::const_iterator it = boundary_node_index.begin(); it != boundary_node_index.end(); ++it, ++iterator_loop2) {
            if ((diff = (*it) - last_index)>0) {
                assert(!skip_line( &FEBio_logfile_manipulate, diff));}
            
            else {
                std::reverse_iterator< const char* > temp_reverse(FEBio_logfile_manipulate);
                assert(!skip_line(&temp_reverse,1-diff));
                FEBio_logfile_manipulate = (temp_reverse).base()+1;
            }
            
            assert( (*it) == strtol(FEBio_logfile_manipulate, &temp_char_pt,10) );
            FEBio_logfile_manipulate = temp_char_pt;
            
            (*iterator_loop2)[0] = strtod(FEBio_logfile_manipulate, &temp_char_pt);
            FEBio_logfile_manipulate = temp_char_pt + 1;
            (*iterator_loop2)[1] = strtod(FEBio_logfile_manipulate, &temp_char_pt);
            FEBio_logfile_manipulate = temp_char_pt + 1;
            (*iterator_loop2)[2] = strtod(FEBio_logfile_manipulate, &temp_char_pt);
            last_index = (*it);
        }
    }

    delete[] FEBio_logfile_start;   ///this is neccessary to prevent memory leak when running many loops containing this function. Very important!!!!!
    
    return boundary_nodes;
    
}

std::string print_time() {
    time_t rawtime;
    time(&rawtime);
    return asctime(localtime(&rawtime));
}