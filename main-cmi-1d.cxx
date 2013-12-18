   /***************************************************************************
                          main.cxx  -  description
                             -------------------
    begin                : 20/01/2012 -8:00
    copyright            : (C) 2012 by Claudio Juan Tessone
    email                : tessonec@ethz.ch
 *********************  ******************************************************/


   
#include "cmutual.h"
#include "base1d.h"

#include <util/tools.h>





int main(int argc, char *argv[]){

  using namespace CTGlobal;

  initialize_program( argc,  argv );
  
  std::vector< double > complete_data_1 = LocalUtils::load_datafile_1d (data_1_filename);
  CTGlobal::query_output("MSG", "'"+data_1_filename+ "' n_lines = " + CTTools::to_string( complete_data_1.size() ) );
  std::vector< double > complete_data_2 = LocalUtils::load_datafile_1d (data_2_filename);
  CTGlobal::query_output("MSG", "'"+data_2_filename+ "' n_lines = " + CTTools::to_string( complete_data_2.size() ) );
  
  long  n_points = complete_data_1.size() ;
  if(n_points != complete_data_2.size() ){
        CTGlobal::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
  } 
  
     std::vector<double> data_1(moving_window_size);
     std::vector<double> data_2(moving_window_size);
     std::vector<double> data_2_lagged(moving_window_size);


     std::fstream* file_abs = 0;
     std::fstream* file_abs_nans = 0;
     file_abs = LocalUtils::create_output_file( CTGlobal::basename, "ABS", "cmi.dat");
     file_abs_nans = LocalUtils::create_output_file( CTGlobal::basename, "ABS", "NONNANS.cmi.dat");

     Mutual::fill_vec_gamma(n_points);
  for(long starting_point=0;starting_point < n_points-moving_window_size-tau_max ; starting_point += 1){
//      std::cerr << starting_point << std::endl;

     std::copy(
              complete_data_1.begin()+starting_point ,
              complete_data_1.begin()+starting_point +moving_window_size,
              data_1.begin() );
     std::copy(
              complete_data_2.begin()+starting_point ,
              complete_data_2.begin()+starting_point +moving_window_size,
              data_2.begin() );
     std::copy(
              complete_data_2.begin()+starting_point +tau_max,
              complete_data_2.begin()+starting_point +moving_window_size +tau_max,
              data_2_lagged.begin()  );

     ConditionalMutualInformation1D mi( data_1,data_2,data_2_lagged,kth_order);
     double cmi = mi.get_mi() ;
     *file_abs << cmi << std::endl;  // "\t" ;
     *file_abs_nans << mi.get_valid_points() << std::endl;  //<< "\t" ;
  }
  return EXIT_SUCCESS;

}
