/*************************************************
 * qhpc.C
 * utilities to generate QHPC codes.             
 * Author: Leonid Pryadko 
 * ***********************************************/
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
//#include "templ.H"
#include "util.h"

using namespace itpp;
using namespace std;

bmat do_cycl(int r, int n, const bvec & vals){
  /* make an r (rows) by n (cols)  circulant mat from vals */
  bmat mat(r,n);
  mat.clear();
  for(int i=0;i<vals.size();i++){
    for(int j=0;j<r;j++){
      mat(j,(j+i)%n) = vals[i];      
    }
  }
  return mat;
}
#if 1

Sparse_Mat<bin> kronecker( const GF2mat  X, const GF2mat  Y){
  Sparse_Mat<bin> result(X.rows() * Y.rows(), X.cols() * Y.cols());
  Mat<bin> cY(Y.rows(),Y.cols());
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


Sparse_Mat<bin> concat(const Sparse_Mat<bin> X, const Sparse_Mat<bin> Y){
  // combine the cols of the two matrices 
  it_assert(X.rows()==Y.rows(),"row counts must be equal");
  Sparse_Mat<bin> ans(X.rows(),X.cols()+Y.cols());
  int j;
  for(j=0;j<X.cols();j++)
    ans.set_col(j,X.get_col(j));
  for(int i=0;i<Y.cols();j++, i++)
    ans.set_col(j,Y.get_col(i));
  return ans;
}
#endif 


#if 0

void rand_perm(PermutationMatrix<Dynamic,Dynamic> & perm, const long n, long debug){
  long i;
  perm.resize(n);
  perm.setIdentity();
  //  if(debug) cout << "rand_perm: n="<< n<< endl;
  for (i = 0; i < n; i++) { /* Knuth algorithm for random permutation */ 
    long j = (i+1)*rand() / (RAND_MAX+1);
    //    if ((i!=j)&&(j>=0) && (j<n))
    swap (perm.indices()[i],perm.indices()[j]);
  }
  return  ;
}
#endif 

/* produce code and save generators */
int main(int argc, char *argv[]){
  Parser p; 
  p.set_silentmode(false); /* shut down echo */

  p.init(argc,argv); /* initialize the parser class on cmd line input */
  if((argc==1)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
    cout<< argv[0]<< ": generate hypergraph product code\n"
	<< "\tusage: qhpc param=val\n"
	<< "\tg1=\"[1 1]\" n1=3 r1=3: specify H1; g2,n2,r2: H2\n"
	<< "\tqhpc=1: generate QHPC, outc=1: save LDPC code\n"
        << "\tdebug=1: bitmap for aux information\n"
	<< "\tstrings: fin1=\"\" fin2=\"\": base fname for reading H1, H2 (.dat)\n"
	<< "\tfout=\"try\": base out fname (try_G.dat, try_P.dat, try.it)\n" 
	<< endl;
    exit (-1);
  }
  int debug=1; p.get(debug,"debug");
  int copy=0;
  if(debug==0) p.set_silentmode(true); /* shut down echo */
  // if(dbg>0) debug|=dbg; else debug^=-dbg;

  LDPC_Parity H1,H2,HH; /* input check matrices */
  GF2mat_sparse mG,mP; /* explicit check and generator matrices */
  int n1=3, r1=3, n2=0, r2=0;  /* default params of the small codes */	 
  string fin1="", fin2="", fout="try"; 

  p.get(fin1,"fin1"); 
  if(fin1.length()){/* read H from the file */
    H1.load_alist(fin1+".dat");   
    n1=H1.get_nvar();
    r1=H1.get_ncheck();
    if(debug&4) // H1 gen info
      cout << "loading H1 from " << fin1+".dat" << endl;
  }
  else{
    p.get(n1,"n1"); p.get(r1,"r1"); 
    bvec g1="1 1";     
    p.get(g1,"g1"); 
    if(g1.length()>0){
      GF2mat_sparse hh1=GF2mat(do_cycl(r1,n1,g1)).sparsify();
      GF2mat_sparse_alist hl1; hl1.from_sparse(hh1);
      H1=LDPC_Parity(hl1);
      if(debug&4) // H1 gen info
	cout << "generating circ H1 from g1=" << g1 << endl;
    }
  }
  if(((H1.get_nvar()<40)&&(debug&1))||(debug&128)) // output H1
    cout << "H1=\n"<< H1.get_H()<< endl;
  if(H1.get_nvar()==0)
    ERROR("failed to create H1 matrix, exiting");
  

  p.get(fin2,"fin2");
  if(fin2.length()){/* read H from the file */
    H2.load_alist(fin2+".dat");
    if(debug&4) // H2 gen info
      cout << "loading H2 from " << fin2+".dat" << endl;
  }
  else{ /* see if generator is defined */
    bvec g2="";     
    p.get(g2,"g2"); 
    if(g2.length()>0){
      p.get(n2,"n2"); 
      p.get(r2,"r2"); 
      if(n2==0) /* same size */
	n2=n1;
      if(r2==0)
	r2=r1; 
      GF2mat_sparse hh2=GF2mat(do_cycl(r2,n2,g2)).sparsify();
      GF2mat_sparse_alist hl2; hl2.from_sparse(hh2);
      H2=LDPC_Parity(hl2);
      if(debug&4) // H2 gen info
	cout << "generating circ H2 from g2=" << g2 << endl;
    }
    else{ /* last resort: copy H1 */
      H2=LDPC_Parity(H1.export_alist());    
      copy=1;
      if(debug&4) // H2 gen info
	cout << "copying H2 from H1" << endl;
    }
  }
  if(!copy)
    if(((H2.get_nvar()<40)&&(debug&1))||(debug&128)) // output H2
      cout << "H2=\n"<< H2.get_H()<< endl;

  int qhpc=1; /* =1: doing QHPC, =0: cyclic binary */
  p.get(qhpc,"qhpc"); 
  
  int nc, nr, ns; /* cols, rows in G, rows in H = # syndromes */		  
  if(qhpc){ // G=Gx=(E*H1,H2^t*E), P=Gz=(H2*E,E*H1^t) 
    if(debug&4) cout << "generating QHPC code"<< endl;
    n1=H1.get_nvar(); r1=H1.get_ncheck();
    n2=H2.get_nvar(); r2=H2.get_ncheck();
    GF2mat EE=do_cycl(n2,n2,bvec("1"));
    GF2mat_sparse M1=kronecker(EE,H1.get_H(false));
    //    cout << "M1 "<< M1.rows()<< " by "<< M1.cols()<< endl;
    EE=do_cycl(r1,r1,bvec("1"));
    GF2mat_sparse M2=kronecker(H2.get_H(true),EE);
    //    cout << "M2 "<< M2.rows()<< " by "<< M2.cols()<< endl;
    mG=concat(M1,M2);  /* Gx generator matrix */
  
    EE=do_cycl(n1,n1,bvec("1"));
    M1=kronecker(H2.get_H(false),EE);
    //    cout << "M1 "<< M1.rows()<< " by "<< M1.cols()<< endl;
    EE=do_cycl(r2,r2,bvec("1"));
    M2=kronecker(EE,H1.get_H(true));
    // cout << "M2 "<< M2.rows()<< " by "<< M2.cols()<< endl;
    mP=concat(M1,M2);
    nc=n1*n2+r1*r2; /* expected */
    nr=n2*r1; 
    ns=n1*r2; 
  }
  else{ /* use g1 for Gx and g2 for Gz */
    if(debug&4) cout << "using G=H1, P=H2"<< endl;
    n1=H1.get_nvar(); r1=H1.get_ncheck();
    n2=H2.get_nvar(); r2=H2.get_ncheck();
    it_assert(n1==n2,"column size should be equal");
    nc=n1; /* expected */
    nr=r1;
    ns=r2;
    mG=H1.get_H();
    mP=H2.get_H();
  }
  it_assert(nc==mG.cols() && nc==mP.cols(), "expect nc vs actual");
  it_assert(nr==mG.rows(), "expect nr vs actual");
  it_assert(ns==mP.rows(), "expect ns vs actual");

  if(((nc<40)&&(debug&1))||(debug&256)){ // output G,P  
    cout << "mG=\n"<< mG<< endl;
    cout << "mP=\n"<< mP<< endl;   
  }
  if(((nc<=1000)&&(debug&1))||(debug&512)){ // compute rank 
    int rrG=GF2mat(mG).row_rank();
    int rrP=GF2mat(mP).row_rank();
    cout << " rank mG="<< rrG
	 << " rank mP="<< rrP << " k="<< nc-rrG-rrP
	 << " [["<< nc<< ","<< nc-rrG-rrP<< "]]"
	 << endl;
  }
  GF2mat_sparse MM=mG*mP.transpose(); 
  if(MM.nnz()>0){ 
    it_warning("G*Pt non-zero!");
    }
  GF2mat_sparse_alist hl; 
  hl.from_sparse(mP);
  HH=LDPC_Parity(hl);
  if(debug&2) HH.display_stats();
  p.get(fout,"fout"); /* out file base name */
  int outc=0; /* do not output code by default */
  p.get(outc,"outc");
  if(outc){
    LDPC_Code C1(&HH);      
    if(debug&2)  cout << C1<< endl; /* display code stats */
    C1.save_code(fout+".it");
  }
  HH.save_alist(fout+"_P.dat");

  hl.from_sparse(mG);
  hl.write(fout+"_G.dat");

  return 0;
}
