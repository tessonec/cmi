#include "mutual.h"


MutualInformation1D::MutualInformation1D(std::vector< double >& dx, std::vector< double >& dy, long int kth):
 data_x(dx),data_y(dy) {
    kth_order = kth;
    n_points = data_x.size() ;
    if(n_points != data_y.size() ){
        LocalUtils::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
    }

    distances_i.resize(n_points);
    
}



double MutualInformation1D::distance(long i, long j){
  
  double dx = metric( data_x.at(i), data_x.at(j) );
  if(is_nan(dx) ) return dx ;

  double dy = metric( data_y.at(i), data_y.at(j) );
  if(is_nan(dy) ) return dy ;

  return dx>dy?dx:dy;
}




double MutualInformation1D::get_kth_distance(long int i){

   distances_i.clear();
   distances_i.reserve(n_points);
   for( long j=0;j<n_points;j++){
      double dij = distance(i,j);
      if( is_nan(dij) )
        continue;
      distances_i.push_back( dij );
   }
   
   if( distances_i.size() <  kth_order+1 ) return -1;
   partial_sort(distances_i.begin(), distances_i.begin()+kth_order+1, distances_i.end());

   return distances_i.at(kth_order);   
  
}

void MutualInformation1D::iterate_once(long int i)
{
   double kth_dist  = get_kth_distance(i);
   if(kth_dist  == -1)
     return;
   
   double n_x = 0;
   double n_y = 0;
   
   for( long j= 0;j < n_points; j++ ){
     double dx = metric( data_x.at(i) , data_x.at(j) );
     if(is_nan(dx) ){
       count_nans+=1;
       continue;  
     }
     double dy = metric( data_y.at(i) , data_y.at(j) );
     if(is_nan(dy) ){
       count_nans+=1;
       continue;  
     }

     if( dx < kth_dist  )
       n_x += 1.0;
     if( dy < kth_dist  )
       n_y += 1.0;

     count_nonnans+=1;
   }
//    accumulator_xz += digamma(n_xz, &foo); // the particle is already counted otherwise, +1
//    accumulator_yz += digamma(n_yz , &foo); // the particle is already counted otherwise, +1
    accumulator_x += Mutual::vec_digamma[n_x]; // the particle is already counted otherwise, +1
    accumulator_y += Mutual::vec_digamma[n_y]; // the particle is already counted otherwise, +1
   

   count_accs+=1;
}


void MutualInformation1D::reset_counters(){
    accumulator_x = 0.;
    accumulator_y = 0.;
    count_accs = 0;
    count_nans = 0;

    count_nonnans = 0;
}

double MutualInformation1D::get_mi()
{
    reset_counters();
    
    for(long i=0;i<n_points;i++)
      iterate_once(i);
    
//     return digamma(kth_order,&foo) - ( accumulator_xz + accumulator_yz - accumulator_z ) / count_accs  ;
    return Mutual::vec_digamma[kth_order] - ( accumulator_x + accumulator_y ) / count_accs  +  Mutual::vec_digamma[n_points] ;
}


//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------






MutualInformation2D::MutualInformation2D(std::vector< Cartesian2DCoordinate >& dx, std::vector< Cartesian2DCoordinate  >& dy, long int kth):
 data_x(dx),data_y(dy) {
    kth_order = kth;
    n_points = data_x.size() ;
    if(n_points != data_y.size() ){
        LocalUtils::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
    }
    n_points --;
    distances_i.resize(n_points);
    
    
}



double MutualInformation2D::distance(long i, long j){
  double dx = metric( data_x.at(i), data_x.at(j) );
  if(is_nan(dx) ) return dx ;

  double dy = metric( data_y.at(i), data_y.at(j) );
  if(is_nan(dy) ) return dy ;

  return dx>dy?dx:dy;
}




double MutualInformation2D::get_kth_distance(long int i){

   distances_i.clear();
   distances_i.reserve(n_points);
   for( long j=0;j<n_points;j++){
      double dij = distance(i,j);
      if( is_nan(dij) )
        continue;
      distances_i.push_back( dij );
   }
   if( distances_i.size() <  kth_order+1 ) return -1;
   partial_sort(distances_i.begin(), distances_i.begin()+kth_order+1, distances_i.end());

   return distances_i.at(kth_order);   
  
}

void MutualInformation2D::iterate_once(long int i)
{

   double kth_dist  = get_kth_distance(i);
   if(kth_dist  == -1)
     return;
   
   double n_x = 0;
   double n_y = 0;
   
   for( long j= 0;j < n_points; j++ ){
     double dx = metric( data_x.at(i), data_x.at(j) );
     if(is_nan(dx) ){
       count_nans +=1. ;
       continue;
     }
     double dy = metric( data_y.at(i), data_y.at(j) );
     if(is_nan(dy) ){
       count_nans +=1. ;
       continue;
     }
     
     count_nonnans += 1;
     if( dx < kth_dist )
       n_x += 1.0;
     if( dy < kth_dist )
       n_y += 1.0;
     
   }
//    accumulator_xz += digamma(n_xz, &foo); // the particle is already counted otherwise, +1
//    accumulator_yz += digamma(n_yz , &foo); // the particle is already counted otherwise, +1
//    accumulator_z += digamma(n_z, &foo); // the particle is already counted otherwise, +1
    accumulator_x += Mutual::vec_digamma[n_x]; // the particle is already counted otherwise, +1
    accumulator_y += Mutual::vec_digamma[n_y]; // the particle is already counted otherwise, +1
   

   count_accs+=1;
}


void MutualInformation2D::reset_counters(){
    accumulator_x = 0.;
    accumulator_y = 0.;
    
    count_accs = 0;
    count_nans = 0;
    count_nonnans = 0;
}

double MutualInformation2D::get_mi()
{
    reset_counters();
//    int foo;
    
    for(long i=0;i<n_points;i++)
      iterate_once(i);
    
//     return digamma(kth_order,&foo) - ( accumulator_xz + accumulator_yz - accumulator_z ) / count_accs  ;
    return Mutual::vec_digamma[kth_order] - ( accumulator_x + accumulator_y ) / count_accs +  Mutual::vec_digamma[n_points] ;
}


double MutualInformation2D::metric ( Cartesian2DCoordinate& pt1, Cartesian2DCoordinate& pt2 )
{
//     std::cerr << "METRIC" << std::endl;  

  // This metric is the norm L2
  double dx = pt1.first - pt2.first ;
  double dy = pt1.second - pt2.second ;
  
  return dx*dx + dy*dy ;
}

double MutualInformation2DLInf::metric ( Cartesian2DCoordinate& pt1, Cartesian2DCoordinate& pt2 )
{
  double dx = fabs( pt1.first - pt2.first );
  double dy = fabs( pt1.second - pt2.second );
  
  return dx > dy? dx:dy;
}

double MutualInformation2DX::metric ( Cartesian2DCoordinate& pt1, Cartesian2DCoordinate& pt2 )
{
  return fabs( pt1.first - pt2.first );
}

double MutualInformation2DY::metric ( Cartesian2DCoordinate& pt1, Cartesian2DCoordinate& pt2 )
{
  return fabs( pt1.second - pt2.second );
}



double MutualInformation2DAngle::metric ( Cartesian2DCoordinate& pt1, Cartesian2DCoordinate& pt2 )
{
  double dx = atan2( pt1.first, pt1.second );
  double dy = atan2( pt2.first, pt2.second );
  return fabs( dx-dy );
}




//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------
//:::~ (CT ) -----------------------------------------------------------------------------------------------






MutualInformation3D::MutualInformation3D(std::vector< Cartesian3DCoordinate >& dx, std::vector< Cartesian3DCoordinate  >& dy, long int kth):
 data_x(dx),data_y(dy) {
    kth_order = kth;
    n_points = data_x.size() ;
    if(n_points != data_y.size() ){
        LocalUtils::query_output("ERR", "data sizes do not coincide... exiting");
        std::exit(EXIT_FAILURE);
    }
    n_points --;
    distances_i.resize(n_points);
    
    
}



double MutualInformation3D::distance(long i, long j){
  double dx = metric( data_x.at(i), data_x.at(j) );
  if(is_nan(dx) ) return dx ;

  double dy = metric( data_y.at(i), data_y.at(j) );
  if(is_nan(dy) ) return dy ;

  return dx>dy?dx:dy;
}




double MutualInformation3D::get_kth_distance(long int i){

   distances_i.clear();
   distances_i.reserve(n_points);
   for( long j=0;j<n_points;j++){
      double dij = distance(i,j);
      if( is_nan(dij) )
        continue;
      distances_i.push_back( dij );
   }
   if( distances_i.size() <  kth_order+1 ) return -1;
   partial_sort(distances_i.begin(), distances_i.begin()+kth_order+1, distances_i.end());

   return distances_i.at(kth_order);   
  
}

void MutualInformation3D::iterate_once(long int i)
{

   double kth_dist  = get_kth_distance(i);
   if(kth_dist  == -1)
     return;
   
   double n_x = 0;
   double n_y = 0;
   
   for( long j= 0;j < n_points; j++ ){
     double dx = metric( data_x.at(i), data_x.at(j) );
     if(is_nan(dx) ){
       count_nans +=1. ;
       continue;
     }
     double dy = metric( data_y.at(i), data_y.at(j) );
     if(is_nan(dy) ){
       count_nans +=1. ;
       continue;
     }
     
     count_nonnans += 1;
     if( dx < kth_dist )
       n_x += 1.0;
     if( dy < kth_dist )
       n_y += 1.0;
     
   }
//    accumulator_xz += digamma(n_xz, &foo); // the particle is already counted otherwise, +1
//    accumulator_yz += digamma(n_yz , &foo); // the particle is already counted otherwise, +1
//    accumulator_z += digamma(n_z, &foo); // the particle is already counted otherwise, +1
    accumulator_x += Mutual::vec_digamma[n_x]; // the particle is already counted otherwise, +1
    accumulator_y += Mutual::vec_digamma[n_y]; // the particle is already counted otherwise, +1
   

   count_accs+=1;
}


void MutualInformation3D::reset_counters(){
    accumulator_x = 0.;
    accumulator_y = 0.;
    
    count_accs = 0;
    count_nans = 0;
    count_nonnans = 0;
}

double MutualInformation3D::get_mi()
{
    reset_counters();
//    int foo;
    
    for(long i=0;i<n_points;i++)
      iterate_once(i);
    
//     return digamma(kth_order,&foo) - ( accumulator_xz + accumulator_yz - accumulator_z ) / count_accs  ;
    return Mutual::vec_digamma[kth_order] - ( accumulator_x + accumulator_y ) / count_accs +  Mutual::vec_digamma[n_points] ;
}


double MutualInformation3D::metric ( Cartesian3DCoordinate& pt1, Cartesian3DCoordinate& pt2 )
{
//     std::cerr << "METRIC" << std::endl;  

  // This metric is the norm L2
  double dx = pt1.x- pt2.x ;
  double dy = pt1.y - pt2.y ;
  double dz = pt1.z - pt2.z ;
  
  return dx*dx + dy*dy + dz*dz ;
}

double MutualInformation3DLInf::metric ( Cartesian3DCoordinate& pt1, Cartesian3DCoordinate& pt2 )
{
  double dx = fabs( pt1.x- pt2.x );
  double dy = fabs( pt1.y- pt2.y );
  double dz = fabs( pt1.z- pt2.z );
  if(dx > dy )
    return dx > dz? dx:dz;
  else
    return dy > dz? dy:dz;
}


