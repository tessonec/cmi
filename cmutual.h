

#include "digamma.h"

#include "utils.h"

#include <vector>

#include <cmath>
#include <cstdlib>

#include <algorithm>


class ConditionalMutualInformation1D {
  
private:
  long kth_order;
  long n_points;
  
  std::vector<double> &data_x;
  std::vector<double> &data_y;
  std::vector<double> &data_z;
  
  std::vector<double> distances_i;
  
  double accumulator_x;
  double accumulator_y;
  double accumulator_z;
  double count_accs;
  double count_nans;
  double count_nonnans;
  
  void reset_counters();

  double get_kth_distance(long i);
  
  virtual double distance(long i, long j);
  virtual double metric( double pt1, double pt2 ){
    return fabs(pt1-pt2);
  }
public:
  
    ConditionalMutualInformation1D( std::vector<double > &dx,std::vector<double > &dy, std::vector<double > &dz, long kth);
  
  ~ConditionalMutualInformation1D(){}
  
  void iterate_once( long i );
  
  long int get_data_length(){ return n_points; }
  
  
  double get_mi();
  
  
  double get_valid_points(){
    return count_nonnans/(count_nonnans+count_nans);
  }
};




class ConditionalMutualInformation2D {
  
protected:
  long kth_order;
  long n_points;
  
  std::vector< Cartesian2DCoordinate > &data_x;
  std::vector< Cartesian2DCoordinate > &data_y;
  std::vector< Cartesian2DCoordinate > &data_z;
  
  std::vector<double> distances_i;
  
  double accumulator_x;
  double accumulator_y;
  double accumulator_z;
  double count_accs;
  double count_nans;
  double count_nonnans;
  
  void reset_counters();

  double get_kth_distance(long i);
  
  virtual double distance(long i, long j);
  virtual double metric( Cartesian2DCoordinate &pt1, Cartesian2DCoordinate &pt2 );

public:
  
    ConditionalMutualInformation2D( std::vector< Cartesian2DCoordinate > &dx,std::vector< Cartesian2DCoordinate > &dy, std::vector< Cartesian2DCoordinate > &dz, long kth);
  
  ~ConditionalMutualInformation2D(){}
  
  void iterate_once( long i );
  
  long int get_data_length(){ return n_points; }

  double get_mi();
  
  double get_valid_points(){
    return count_nonnans/(count_nonnans+count_nans);
  }
};



class MutualInformation2DLInf : public ConditionalMutualInformation2D {

protected:  
  virtual double metric( Cartesian2DCoordinate &pt1, Cartesian2DCoordinate &pt2 );

public:
  
  MutualInformation2DLInf( std::vector< Cartesian2DCoordinate > &dx, std::vector< Cartesian2DCoordinate > &dy,  std::vector< Cartesian2DCoordinate > &dz, long kth):
        ConditionalMutualInformation2D(dx,dy,dz,kth){};
  ~MutualInformation2DLInf(){}  
};


class MutualInformation2DX : public ConditionalMutualInformation2D {

protected:  
  virtual double metric( Cartesian2DCoordinate &pt1, Cartesian2DCoordinate &pt2 );

public:
  
  MutualInformation2DX( std::vector< Cartesian2DCoordinate > &dx,std::vector< Cartesian2DCoordinate > &dy, std::vector< Cartesian2DCoordinate > &dz, long kth):
        ConditionalMutualInformation2D(dx,dy,dz,kth){};
  ~MutualInformation2DX(){}
  
};


class MutualInformation2DY : public ConditionalMutualInformation2D {

protected:  
  virtual double metric( Cartesian2DCoordinate &pt1, Cartesian2DCoordinate &pt2 );

public:
  
  MutualInformation2DY( std::vector< Cartesian2DCoordinate > &dx,std::vector< Cartesian2DCoordinate > &dy, std::vector< Cartesian2DCoordinate > &dz, long kth):
        ConditionalMutualInformation2D(dx,dy,dz,kth){};
  ~MutualInformation2DY(){}
};



class MutualInformation2DAngle : public ConditionalMutualInformation2D {

protected:  
  virtual double metric( Cartesian2DCoordinate &pt1, Cartesian2DCoordinate &pt2 );

public:
  
  MutualInformation2DAngle( std::vector< Cartesian2DCoordinate > &dx,std::vector< Cartesian2DCoordinate > &dy, std::vector< Cartesian2DCoordinate > &dz, long kth):
        ConditionalMutualInformation2D(dx,dy,dz,kth){};
  ~MutualInformation2DAngle(){}  
};




class ConditionalMutualInformation3D {
  
protected:
  long kth_order;
  long n_points;
  
  std::vector< Cartesian3DCoordinate > &data_x;
  std::vector< Cartesian3DCoordinate > &data_y;
  std::vector< Cartesian3DCoordinate > &data_z;
  
  std::vector<double> distances_i;
  
  double accumulator_x;
  double accumulator_y;
  double accumulator_z;
  double count_accs;
  double count_nans;
  double count_nonnans;
  
  void reset_counters();

  double get_kth_distance(long i);
  
  virtual double distance(long i, long j);
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
    ConditionalMutualInformation3D( std::vector< Cartesian3DCoordinate > &dx,std::vector< Cartesian3DCoordinate > &dy, std::vector< Cartesian3DCoordinate > &dz, long kth);
  
  ~ConditionalMutualInformation3D(){}
  
  void iterate_once( long i );
  
  long int get_data_length(){ return n_points; }

  double get_mi();
  
  double get_valid_points(){
    return count_nonnans/(count_nonnans+count_nans);
  }
};



class MutualInformation3DLInf : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DLInf( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DLInf(){}  
};



class MutualInformation3DX : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DX( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DX(){}  
};




class MutualInformation3DY : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DY( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DY(){}  
};




class MutualInformation3DZ : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DZ( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DZ(){}  
};




class MutualInformation3DL2XY : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DL2XY( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DL2XY(){}  
};




class MutualInformation3DL2XZ : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DL2XZ( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DL2XZ(){}  
};




class MutualInformation3DL2YZ : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DL2YZ( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DL2YZ(){}  
};




class MutualInformation3DLInfXY : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DLInfXY( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DLInfXY(){}  
};




class MutualInformation3DLInfXZ : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DLInfXZ( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DLInfXZ(){}  
};




class MutualInformation3DLInfYZ : public ConditionalMutualInformation3D {

protected:  
  virtual double metric( Cartesian3DCoordinate &pt1, Cartesian3DCoordinate &pt2 );

public:
  
  MutualInformation3DLInfYZ( std::vector< Cartesian3DCoordinate > &dx, std::vector< Cartesian3DCoordinate > &dy,  std::vector< Cartesian3DCoordinate > &dz, long kth):
        ConditionalMutualInformation3D(dx,dy,dz,kth){};
  ~MutualInformation3DLInfYZ(){}  
};



