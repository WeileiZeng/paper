#ifndef templ_H
#define templ_H
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
// #include <itpp/base/gf2mat.h>
#include <sstream>
#include <stdio.h>


//! Return the postion of the min element in the vector by absolute value
template<class T>
int min_abs_index(const itpp::Vec<T> &in)
{
  int minindex = 0;
  for (int i = 1; i < in.length(); i++)
    if (abs(in[i]) < abs(in[minindex]))
      minindex = i;
  return minindex;
}

template<class T>
int weight(const itpp::Vec<T> &in)
{
  int w = 0;
  for (int i = 1; i < in.length(); i++)
    if (in[i])
      w++;
  return w;
}

template<class T>
itpp::Sparse_Mat<T> kronecker( const itpp::Mat<T>  X, const itpp::Mat<T>  Y){
  itpp::Sparse_Mat<T> result(X.rows() * Y.rows(), X.cols() * Y.cols());
  itpp::Mat<T> cY(Y.rows(),Y.cols());
  for (int i=0;i<Y.rows(); i++)
    for (int j=0;j<Y.cols(); j++)
      cY(i,j)=Y(i,j);
  for (int i = 0; i < X.rows(); i++)
    for (int j = 0; j < X.cols(); j++){
      if(X(i,j))
	result.set_submatrix(i * Y.rows(), j * Y.cols(), cY);
    }
  result.compact();
  return result;
}

template<class T>
itpp::Sparse_Mat<T> concat(const itpp::Sparse_Mat<T> X, const itpp::Sparse_Mat<T> Y){
  // combine the cols of the two matrices 
  it_assert(X.rows()==Y.rows(),"row counts must be equal");
  itpp::Sparse_Mat<T> ans(X.rows(),X.cols()+Y.cols());
  int j;
  for(j=0;j<X.cols();j++)
    ans.set_col(j,X.get_col(j));
  for(int i=0;i<Y.cols();j++, i++)
    ans.set_col(j,Y.get_col(i));
  return ans;
}


/* permutation routines */
void triv_perm(itpp::Vec<int> & p, const int n){ /* return identity permutation */
  p.set_length(n);
  for(int i=0;i<n;i++)
    p[i]=i;
}

template<class T>
void permute(itpp::Vec<T> & v, const itpp::Vec<int> & p, bool inv = false){  
  /* permute elements of vector v according to p; no change to p.  
   * Optionally, use inverse of p if inv !=0 */
  int n=length(p);  
  it_assert(length(v) == n,"permute(): dimensions do not match");
  itpp::Vec<T> ans(n);
  for (int j = 0; j < n; j++) {
    if (inv == 0) 
      ans.set(j,v[p[j]]);    
    else 
      ans.set(p[j], v[j]);    
  }
  v=ans;
}

#define BIT(x,i) (((x)>>(i))&1)



#endif /*  templ_H */
