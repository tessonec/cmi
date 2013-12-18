   /***************************************************************************
                          main.cxx  -  description
                             -------------------
    begin                : 20/01/2012 -8:00
    copyright            : (C) 2012 by Claudio Juan Tessone
    email                : tessonec@ethz.ch
 *********************  ******************************************************/


   
#include "base2d.h"
// :::~ mutual.h VA ABAJO DE base2d
#include "cmutual.h" 

#include <util/tools.h>


#include <fstream>
#include <algorithm>

using namespace std;


void calc_mutual_information(ConditionalMutualInformation2D *mi, std::fstream *file_out, std::fstream *file_out_nans, std::string metric = "NONE", long starting_point = -1, long lag = -1){
        using namespace CTGlobal;
      
//         std::cerr << metric << ": (" << starting_point<< ", " << lag  << ") " << std::endl;
        *file_out << mi->get_mi() << "\t" ;
        *file_out_nans << mi->get_valid_points() << "\t" ;
        
        if( lag >= tau_max-tau_step ){
          *file_out << std::endl;
          *file_out_nans << std::endl ;
          file_out->flush();file_out_nans->flush();
        }
        
}



/*
int main(int argc, char *argv[]){

  using namespace CTGlobal;

  initialize_program( argc,  argv );
  
  std::vector<double> complete_data_x = load_datafile(data_x_filename);
  std::vector<double> complete_data_y = load_datafile(data_y_filename);
  
  long  n_points = complete_data_x.size() ;
  if(n_points != complete_data_y.size() ){
        CTGlobal::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
  }
  n_points --;
  
  
      MutualInformation mi(complete_data_x,complete_data_y,kth_order);
    for(long lag=0;lag< tau_max; lag+= tau_step){
      if(store_tau_values)
        *f_tau_values_out << lag << std::endl;

      std::cout << mi.get_mi(lag) << "\n" ;
    }
  return EXIT_SUCCESS;

}


*/

int main(int argc, char *argv[]){

  using namespace CTGlobal;

  initialize_program( argc,  argv );
  
  std::vector< Cartesian2DCoordinate > complete_data_1 = LocalUtils::load_datafile_2d(data_1_filename);
  CTGlobal::query_output("MSG", "'"+data_1_filename+ "' n_lines = " + CTTools::to_string( complete_data_1.size() ) );
  std::vector< Cartesian2DCoordinate > complete_data_2 = LocalUtils::load_datafile_2d(data_2_filename);
  CTGlobal::query_output("MSG", "'"+data_2_filename+ "' n_lines = " + CTTools::to_string( complete_data_2.size() ) );
  
  long  n_points = complete_data_1.size() ;
  if(n_points != complete_data_2.size() ){
        CTGlobal::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
  }
  n_points --;
  Mutual::fill_vec_gamma(n_points);  
  
  std::vector< Cartesian2DCoordinate > data_1 (moving_window_size);
  std::vector< Cartesian2DCoordinate > data_2 (moving_window_size);
  std::vector< Cartesian2DCoordinate > data_2_lagged (moving_window_size);
//   std::cerr << "1" << std::endl;  
  if(store_tau_values){
    for(long lag=0;lag< tau_max; lag+= tau_step)
        *f_tau_values_out << lag << std::endl;
    f_tau_values_out->flush();
  }
//   std::cerr << "1" << std::endl;  
  
  if(store_time_window_values){
       for(long starting_point=0;starting_point < n_points-moving_window_size +tau_max; starting_point += moving_window_step)
         *f_time_window_values_out << starting_point << std::endl;
       f_time_window_values_out->flush();  
  }

  

  fstream* file_cartesian = 0;
  fstream* file_cartesian_nans = 0;
  fstream* file_abs = 0;
  fstream* file_abs_nans = 0;
  fstream* file_x = 0;
  fstream* file_x_nans = 0;
  fstream* file_y = 0;
  fstream* file_y_nans = 0;
  fstream* file_angle = 0;
  fstream* file_angle_nans = 0;
  
  if( CTGlobal::metric == "L2" || CTGlobal::metric == "ALL" ){
     file_cartesian = LocalUtils::create_output_file( CTGlobal::basename, "L2", "cmi.dat");
     file_cartesian_nans = LocalUtils::create_output_file( CTGlobal::basename, "L2", "NONNANS.cmi.dat");
  }
  
  if( CTGlobal::metric == "LINF" || CTGlobal::metric == "ALL" ){
     file_abs = LocalUtils::create_output_file( CTGlobal::basename, "LINF", "cmi.dat");
     file_abs_nans = LocalUtils::create_output_file( CTGlobal::basename, "LINF", "NONNANS.cmi.dat");
  }
  
  if( CTGlobal::metric == "X" || CTGlobal::metric == "ALL" ){
     file_x = LocalUtils::create_output_file( CTGlobal::basename, "X", "cmi.dat");
     file_x_nans = LocalUtils::create_output_file( CTGlobal::basename, "X", "NONNANS.cmi.dat");
  }
  if( CTGlobal::metric == "Y" || CTGlobal::metric == "ALL" ){
     file_y = LocalUtils::create_output_file( CTGlobal::basename, "Y", "cmi.dat");
     file_y_nans = LocalUtils::create_output_file( CTGlobal::basename, "Y", "NONNANS.cmi.dat");
  }
  if( CTGlobal::metric == "ANGLE" || CTGlobal::metric == "ALL" ){
     file_angle = LocalUtils::create_output_file( CTGlobal::basename, "ANGLE", "cmi.dat");
     file_angle_nans = LocalUtils::create_output_file( CTGlobal::basename, "ANGLE", "NONNANS.cmi.dat");
  }
  
  for(long starting_point=0;starting_point < n_points-moving_window_size-tau_max ; starting_point += moving_window_step){
    std::copy( complete_data_1.begin()+starting_point,complete_data_1.begin()+ starting_point + moving_window_size, data_1.begin() );  
    std::copy( complete_data_2.begin()+starting_point,complete_data_2.begin()+ starting_point + moving_window_size, data_2.begin() );  
    
    for(long lag=0;lag< tau_max; lag+= tau_step){
       std::cerr << metric << ": (" << starting_point<< ", " << lag  << ") " << std::endl;

       std::copy( complete_data_2.begin()+starting_point+lag,complete_data_2.begin()+ starting_point + moving_window_size+lag, data_2_lagged.begin() );  
         
       if( file_cartesian != 0 ){
         ConditionalMutualInformation2D mi ( data_1,data_2,data_2_lagged,kth_order) ;
         calc_mutual_information( &mi, file_cartesian, file_cartesian_nans, "L2", starting_point, lag);
       }

       if( file_abs != 0 ){
         MutualInformation2DLInf mi ( data_1,data_2,data_2_lagged,kth_order) ;
         calc_mutual_information( &mi, file_abs, file_abs_nans, "LINF", starting_point, lag);
       }
    
       if( file_x != 0 ){
         MutualInformation2DX mi ( data_1,data_2,data_2_lagged,kth_order) ;
         calc_mutual_information( &mi, file_x, file_x_nans, "X", starting_point, lag);
       }
    
       if( file_y != 0 ){
         MutualInformation2DY mi ( data_1,data_2,data_2_lagged,kth_order) ;
         calc_mutual_information( &mi, file_y, file_y_nans, "Y", starting_point, lag);
       }
    
       if( file_angle != 0 ){
         MutualInformation2DAngle mi ( data_1,data_2,data_2_lagged,kth_order) ;
         calc_mutual_information( &mi, file_angle, file_angle_nans, "ANGLE", starting_point, lag);
       }
    }
  }
  return EXIT_SUCCESS;

}

