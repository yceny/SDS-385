// [[Rcpp::depends(RcppEigen)]]
# include <RcppEigen.h>
# include <Rcpp.h>
# include <math.h>
# include <list>

using namespace Rcpp;
using Eigen::SparseMatrix;
using Eigen::MappedSparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXi;
typedef Eigen::VectorXd Vd;
typedef Eigen::VectorXi Vi;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix< double > SpMat;
typedef MSpMat::InnerIterator InIterMat;
typedef List list;
typedef Eigen::Triplet<double> T;


double euclidean_distance(double lon1, double lat1,double lon2,double lat2){
  double s = pow((lon1 - lon2),2) + pow((lat1 - lat2),2);
    double ed = sqrt(s);
    return(ed);
}


// [[Rcpp::export]]
SpMat mymain(float alpha, NumericVector lat, NumericVector lon){
	std::list<T> L;
	
	int length = lat.size();
	int index = 0;
	for (int i = 0; i < length - 1; i++ ){
   	 for (int j = i+1; j < length; j++){
            double lon1 = lon[i];
            double lon2 = lon[j];
            double lat1 = lat[i];
            double lat2 = lat[j];
            double dist = euclidean_distance(lon1,lat1,lon2,lat2);
            dist = exp(-dist/alpha);
            if (dist > 0.01 ){
				L.push_back(T(index, i, dist));
				L.push_back(T(index, j, -dist));
				index++;
            }
        
    	}
     }
     int nrows = L.size()/2;
     int ncols = length;
     SpMat D(nrows, ncols);
	 D.setFromTriplets(L.begin(), L.end());
	 return(D);
}



/*** R

#To run this script, use:
#library(RcppEigen)
#sourceCpp("D.cpp")

#prepare lat and lon here or make sure they are loaded in R

D = mymain(alpha, lat, lon)

*/
