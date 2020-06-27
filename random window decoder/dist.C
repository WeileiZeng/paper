/************************************************************************ 
 * routines to calculate distance of an LDPC code (quantum or classical) 
 * author: Leonid Pryadko
 ************************************************************************/
#include <itpp/itcomm.h>
#include <sstream>
using namespace std;
using namespace itpp;

void addto(bvec &x, Sparse_Vec<bin> y){  /* x+=y */
  int i, nnz=y.nnz();
  for(int p=0;p<nnz;p++){
    i=y.get_nz_index(p);
    x[i]+=y.get_nz_data(p);
  }
}

int prep_neis(int z0, vector<int> & nei, const bvec & v, const bvec & s,
	      LDPC_Parity & P, const LDPC_Parity & G){ 
  /* given current value bits and syndrome bits, see which bits in the
     vicinity of the check z0 are unhappy.
     Return number of entries to process (length of nei) */
  int cnt=0;
  Sparse_Vec<bin> row=P.get_row(z0); // list of non-zero nodes for syndrome z0
  //  cout<< "z0="<<z0<< " row="<< row.full()<< endl;
  for(int i=0;i<row.nnz();i++){     
    if(row.get_nz_data(i) && !v[row.get_nz_index(i)]){ /* never flip a bit twice */
      nei.push_back(row.get_nz_index(i));
      cnt++;
    }
  }
  return cnt;
}


int start_rec(int w, int wmax, bvec & v, bvec & s,
	      LDPC_Parity & P, LDPC_Parity & G){
  /* recursive function. return: -1: fail; 0: success; 1: termination */
  int res=0, all_zero=1;
  //  cout <<"w="<< w<< "    v="<< v<<" s="<< s << " Pv="<< P.get_H()*v<< endl;
  // if(w==1) cout <<"w="<< w<< "    v="<< v<<" s="<< s<< endl;
  for(int i=0;i<s.size(); i++)  /* check bad syndromes */
    if(s[i]){ /* found a syndrome to fix */
      if(w>=wmax)
	return -1; /* failed to find a cluster */
      all_zero=0;
      vector<int> nei;
      prep_neis(i,nei,v,s,P,G);
      //      cout << " row="<< i<< " nei=" << nei<< endl;
      for(unsigned int p=0;p< nei.size();p++){
	int j=nei[p];
	v[j]=1;
	bvec s1=s;
	addto(s1,P.get_col(j)); 
	res=start_rec(w+1,wmax,v,s1,P,G);
	if(res==1)
	  return 1; /* just get out fast */
	v[j]=0; /* clean up */
      }      
    }
  if(all_zero){
    it_assert(w==wmax," looking for min weight cluster!");
    /* todo: check for G; see what to return in that case */
    //    cout << "found vec="<< v<< endl;
    //    cout << "      s="<< s<< endl;
    return 1; /* success */
  }
  return 0; /* keep going */
}

int do_dist_clus(int wmax, LDPC_Parity & P, LDPC_Parity & G){
  int nc=P.get_nvar(), ns=P.get_ncheck(); 
  bvec v(nc),s(ns); v.zeros(); s.zeros(); // value and syndrome vecs 
  for(int w=2;w<=wmax;w++){ // cluster weight loop 
    //    cout << "# starting w="<< w<< endl;
    for(int i=0;i<nc-w;i++){ // starting bit loop 
      v.zeros();s.zeros();
      v[i]=1;
      //      Sparse_Vec<bin> col=mP.get_col(i);
      addto(s,P.get_col(i)); //  s+=col;
      int done=start_rec(1,w,v,s,P,G);
      if(done==1){
	//	cout << "distance="<< weight(v)<< endl;
	//	cout << "prod="<< P.get_H()*v<< endl;
	return w;
      }
    }
  }
  return -wmax; /* failed up to wmax */
}
#ifdef DEBUG
int main(int argc, char **argv){
  Parser p; 
  p.set_silentmode(false); /* shut down echo */

  p.init(argc,argv); /* initialize the parser class on cmd line input */
  string finG="try_G.dat", finP="try_P.dat"; 
  int stab=0; p.get(stab,"stab"); // classical; stab=1 for quantum 
  p.get(finP,"finP"); if(stab) p.get(finG,"finG");
  LDPC_Parity mP, mG; /* check matrix, optional generator for stab codes */
  mP.load_alist(finP); 
  // cout << "# loaded P("<<ns<<","<<nc<<") matrix"<< endl;
  // cout << "P=\n"<< mP.get_H()<< endl;
  if(stab){
    mG.load_alist(finG);
    //    nr=mG.get_ncheck();
    it_assert(mP.get_nvar()==mG.get_nvar(),"col count in G and P must match");
    //    cout << "# loaded G("<<nr<<","<<nc<<") matrix"<< endl;
  }
  //  Sparse_Mat<bin> P=mP.get_H(), G=mG.get_H();

  int wmax=5; p.get(wmax,"wmax");
  int d=do_dist_clus(wmax,mP,mG);
  if(d>0)
    cout << "success d="<< d<< endl;
  else 
    cout << "failed up to wmax="<< -d<< endl;
    
  return 0;
}
#endif /* DEBUG */
