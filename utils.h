
#include <util/tools.h>

#include <cmath>
#include <fstream>
#include <algorithm>
#include <locale>
#include <vector>
#include <sstream>

#include "digamma.h"

typedef std::pair<double,double> Cartesian2DCoordinate ;


struct Cartesian3DCoordinate {
  double x;
  double y;
  double z;
  
  Cartesian3DCoordinate(double xxx=0, double yyy=0, double zzz=0){
    x=xxx;y=yyy;z=zzz;
  }
  ~Cartesian3DCoordinate(){}
};





bool is_nan(double var);
void read_next(std::istream &fin, double &ret);

namespace Mutual {
  extern std::vector<double> vec_digamma ;

  void fill_vec_gamma(long int n_points);

};


namespace LocalUtils {
  void query_output(std::string type, std::string msg , std::string oc = "<>", int indent = 0);
  
  std::fstream *create_output_file(std::string basename, std::string app, std::string extension = "dat");
  
  std::vector<double> load_datafile_1d (std::string fname);
  std::vector< Cartesian2DCoordinate > load_datafile_2d(std::string fname);
  std::vector< Cartesian3DCoordinate > load_datafile_3d(std::string fname);
  
};  



   






