//
//  file_processing.h
//  
//
//  Created by Yu Quan on 23/8/15.
//
//

#ifndef _file_processing_h
#define _file_processing_h

#include <boost/multi_array.hpp>
#include <cstring>
#include <fstream>
#include <vector>
#include "defines.h"


///isnewline() function
#define isnewline(X) (isspace(X) && (!isblank(X)))

///Macros; extract space characters until the next non-space character in stream
#define IGNORE_SPACE(X) \
{ while (isspace( (X).peek()))   \
(X).get(); }

///extract non-space characters until the next space character in stream
#define IGNORE_NON_SPACE(X) \
{ while (!isspace( (X).peek()))   \
(X).get(); }

///Macros for binary read/write.
///#define WRITE_BINARY(STREAM,X) ( STREAM.write(reinterpret_cast<const char*> (&(X)), sizeof(X)) )

///#define READ_BINARY(STREAM,X) ( STREAM.read(reinterpret_cast<char*> (&(X)), sizeof(X)) )

///typedef for boost::multi_array
typedef boost::multi_array<floatT,3> floatTArray3D;
typedef boost::multi_array<floatT,2> floatTArray2D;

typedef unsigned long long int FileSizeT;  ///this type has the limit of at least 1.8E19. Enough to represent a file size that is millions of TB :D

///find the size of a file in bytes.
FileSizeT get_filesize(const char *const filename);

///get c-style string from a file
char* get_filestring(const char *const filename);

///this version also retrives the file size
char* get_filestring(const char *const filename, FileSizeT* file_size);


///these are for xml input format
char* begin_indicator(const char *const indicator_str);

char* end_indicator(const char *const indicator_str);

const char* into_section(const char* main_string, const char *const indicator_str);

const char* end_section(const char* main_string, const char *const indicator_str);

///reads a numerical property
floatT read_attribute( const char *& string, const char * attribute );

///move the stream cursor to ignore the comments (xml only)
char* line_wo_comment(std::istream& stream);


///for file names.
void set_extension(char* const filename, const char* const extension);

///insert the postfix according to the int, _0, _1, ..., to the file name, before the extension. abc.txt becomes abc_0.txt, if postfix = 0.
void insert_postfix(char file_name[], int postfix);

///file copy
void file_copy(const char* const destintation_file, const char* const source_file);



///provide exception if specified input file failed to open.
void check_openfile( std::ofstream& stream );
void check_openfile( std::ifstream& stream );

///skip a particular number of lines. Moves the iterator to the correct position.
template <class Iterator_type>
bool skip_line( Iterator_type* string, const int num_lines = 1,
               const char dilimiter = '\n', const int max_line_length = 9999);

///for particular files
///struct of Matlab_output file
struct Matlab_output_data {
    int n_nodes;
    int n_boundary_nodes;
    int n_elements;
    
    std::vector<int> boundary_node_index;
    floatTArray2D initial_nodes;
    floatTArray2D initial_boundary_nodes;
};

///Matlab_output file
const struct Matlab_output_data read_Matlab_output( const char* const file_name );

///FEBio log file
const floatTArray3D read_FEBio_log( const char* file_name , const int n_boundary_nodes, const std::vector<int> boundary_node_index);

void read_true_data(const char* const file_name , floatTArray2D* true_data_left, floatTArray2D* true_data_right, const int num_steps);

///the time and date in text form
std::string print_time();

#endif
