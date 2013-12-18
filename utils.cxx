#include "utils.h"

using namespace std;

std::vector<double> Mutual::vec_digamma;


bool is_nan(double var){
    volatile double d = var;
    return d != d;
}

void Mutual::fill_vec_gamma(long int n_points){
  vec_digamma.resize(1+n_points);
  int foo;
  for(long i=0;i<1+n_points;i++){
    vec_digamma[i] = digamma(i,&foo);
  }
}

void read_next(std::istream &fin, double &ret){
  std::string foo;
  std::stringstream os;
  fin >> foo;
  std::transform(foo.begin(), foo.end(), foo.begin(), ::tolower);
  
//   std::cout << "'"<< foo.c_str()   << "' ";
  
  if ( foo.substr(0,3).c_str()  == string("nan") || foo.substr(0,4).c_str()  == string("-nan") ){
//     std::cout << "here: " ;
    ret = sqrt(-1);
    return;
  }
  os << foo;
  os >> ret;
  
//   std::cout << "'"<< ret << "' ";
  
}


std::vector<double> LocalUtils::load_datafile_1d (std::string fname){
  std::vector<double> vec_ret;
  std::fstream fin(fname.c_str(), ios::in );
  long n_lines = 0;
  double foo;
  while( ! fin.eof() ){
    n_lines ++; 
    read_next(fin, foo) ;
  }
  fin.close();
  
  vec_ret.resize( --n_lines  );

  std::fstream fin2(fname.c_str(), ios::in );
  for( long iline =0;iline < n_lines;iline++){
       read_next( fin2, foo );
       
       vec_ret.at( iline )=foo ;
  }
//    std::cout  << vec_ret.size() << std::endl;
  
  return vec_ret;
}



std::fstream *LocalUtils::create_output_file(std::string basename, std::string app, std::string extension ){
    std::stringstream fout_name;
    fout_name << basename << "_";
    fout_name << app ;
    fout_name << "." << extension;
    return new std::fstream( fout_name.str().c_str(), std::ios::out) ;
}



vector< Cartesian2DCoordinate > LocalUtils::load_datafile_2d(std::string fname){
  std::vector< Cartesian2DCoordinate > vec_ret;
  std::fstream fin(fname.c_str(), ios::in );
  long n_lines = 0;
  double foo_x, foo_y;
  while( ! fin.eof() ){
    n_lines ++; 
    read_next( fin , foo_x);
    read_next( fin , foo_y);
  }
  fin.close();
  vec_ret.resize( -- n_lines  );
  
  
  std::fstream fin2(fname.c_str(), ios::in );
  
  for( long iline =0;iline < n_lines;iline++){
    read_next( fin2, foo_x );
    read_next( fin2, foo_y );
    vec_ret.at( iline ) =  Cartesian2DCoordinate(foo_x, foo_y) ;
  }
  return vec_ret;
  
}



vector< Cartesian3DCoordinate > LocalUtils::load_datafile_3d(std::string fname){
  std::vector< Cartesian3DCoordinate > vec_ret;
  std::fstream fin(fname.c_str(), ios::in );
  long n_lines = 0;
  double foo_x, foo_y,foo_z;
  while( ! fin.eof() ){
    n_lines ++; 
    read_next( fin , foo_x);
    read_next( fin , foo_y);
    read_next( fin , foo_z);
  }
  fin.close();
  vec_ret.resize( -- n_lines  );
  
  
  std::fstream fin2(fname.c_str(), ios::in );
  
  for( long iline =0;iline < n_lines;iline++){
    read_next( fin2, foo_x );
    read_next( fin2, foo_y );
    read_next( fin2, foo_z );
    vec_ret.at( iline ) =  Cartesian3DCoordinate(foo_x, foo_y, foo_z) ;
  }
  return vec_ret;
  
}





void LocalUtils::query_output(std::string type, std::string msg , std::string oc, int indent){
  for(int i=0;i<indent;i++) std::cerr << " ";
  std::cerr << oc[0] << " ";
  std::cerr << "          " << " - "<< type << " ";
  std::cerr << oc[1] << " ";
  std::cerr << msg << std::endl;
}


