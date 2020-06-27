/************************************************************************ 
 * random window algorithm to calculate distance of a code
 * (quantum or classical) author: Leonid Pryadko
 * first version: just classical binary and
 * 	quantum weakly self-orthogonal codes
 ************************************************************************/
#include <itpp/itcomm.h>
#include <sstream>
#include "templ.H"
using namespace std;
using namespace itpp;

typedef struct {
  int debug;
  int seed; 
  int mode; /* 0: binary; (the rest is not implemented)
	       1: quant weakly-self-orthogonal; 2: quant CSS; 3: quant stab */
  int dtry; /* how many tries in the random window */    
} dist_params_t;

dist_params_t q={1,0,-1,100};

int do_gen_matr(GF2mat & H, /* input: orig parity check; 
			     * output: H with permutted cols */
		GF2mat & P, /* output: systematic form of H, generator G */
		GF2mat & G, ivec & perm){  /* output: permutation */
  /* Return rank of H == number of rows of P. */ 
  int r, n=H.cols(); 
  GF2mat T;
  r = H.T_fact(T, P, perm);
  H.permute_cols(perm,false); /* they are now in the same order as for U */
  P.set_size(r,n,true); 
  GF2mat M=P.get_submatrix(0,0,r-1,r-1);
  P=M.inverse()*P; /* now P=(I,A) */
  G=P.get_submatrix(0,r,r-1,n-1).transpose(); 
  G=G.concatenate_horizontal(eye_b(n-r)); /* now G=(AT,I) systematic */
  return r;
}

/* random window algorithm from given generator matrix. */
int do_dist_rand(const GF2mat & G0, const GF2mat & P0,
		 /* input: generator matrix for codewords; 
		  * for quantum codes, check matrix in systematic form. */ 		 
		 int dmin, int rmax){  
  /* If found w smaller than dmin, 
   * or found w equal to dmin but rank(G0)>= rmax, we are not interested
   */
  int r=G0.rows(), n=G0.cols();   /* dimensions of generator matrix */
  int wmin=n+1;
  GF2mat G=G0,T,U; 
  //  cout << G*H.transpose()<< endl;
  for(int l=0;l<q.dtry;l++){ 
    ivec perm, ind = sort_index(randu(n)); // random permutation 
    G.permute_cols(ind,false); 
    r = G.T_fact(T, U, perm);  /* gauss elimination trying for these cols */
    U.set_size(r,n,true); 
    GF2mat M=U.get_submatrix(0,0,r-1,r-1);
    G=M.inverse()*U; /* now G=(I,A) systematic form of G */
    G.permute_cols(ind,true); /* undo row permutation */
    G.permute_cols(perm,true); /* now G has the same col order as orig P0 */
    //    cout << "l="<<l<<" "<< G*H.transpose()<< endl;
    //    cout << "l="<<l<<" "<< G<<endl;
    for(int i=0;i<r;i++){
      bvec x=G.get_row(i), y;      
      int w=weight(x);
      if (w<wmin){
	y=x;
	if(q.mode){ /* quantum code: check for linear dependence with P */
	  for(int j=0;j<P0.rows();j++){
	    if (x[j])
	      x+=P0.get_row(j);
	  }
	  if(!any(x))
	    continue;  /* orig row was trivial */
	}
	if(w<wmin){
	  wmin=w; /* weight-one encoded vectors */	      
	  if (q.debug&2)
	    cout<<"# l="<< l<< " i="<< i <<" w="<< w<< " vec=\n"<<y<< endl;	  
	}
      }
    }
  }
  if((wmin<dmin)||((wmin==dmin)&&(r>rmax))){
    return (-wmin); /* not interesting */
  }  
  return wmin;
}

#ifdef DEBUG
int main(int argc, char **argv){
  Parser p; 
  Real_Timer timer;
  
  p.init(argc,argv);
  p.set_silentmode(false); /* shut down echo */
  if((argc==1)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
    cout<< argv[0]<< ": calc distance of a code using random window algorithm\n"
	<< "\tusage: "<< argv[0]<< " param=val\n"
        << "\tdebug=1: bitmap for aux information\n"
	<< "\tstring: fin=\"try\": base name for reading the code (try.dat) \n"
	<< "\t mode=0: classical binary, 1: weakly self-orthogonal\n"
	<< "\t\t 2: CSS code, 3: stabilizer code (not implemented).\n"
	<< "\t default: fin=try dtry=100 seed=-1 mode=-1 (auto)\n"
	<< endl;
    exit (-1);
  }
  p.set_silentmode(true); /* shut down echo */
  p.get(q.debug,"debug"); 
  if(q.debug){
    cout << "# debug="<< q.debug << endl;
    p.set_silentmode(false);
  }    
  if(q.debug) cout<< "# "; q.dtry=100;    p.get(q.dtry,"dtry");
  if(q.debug) cout<< "# "; q.seed=-1;     p.get(q.seed,"seed");
  if(q.debug) cout<< "# "; q.mode=-1;     p.get(q.mode,"mode");
  if(q.seed==-1)
    RNG_randomize();
  else
    RNG_reset(q.seed);

  if(q.debug) cout<< "# "; string fin="try"; p.get(fin,"fin");
  GF2mat_sparse_alist alist;

  alist.read(fin+".dat");
  GF2mat H=alist.to_sparse(false); /* do not transpose */
  int r=H.rows(), n=H.cols();
  if(q.debug)
    cout << "# calculating the distance of a code with ("
	 << r << " x " << n << ") check matrix  "<< endl;
  if(q.debug&8) // type out matrix 
    cout << H<< endl;
 
  if(q.mode!=0){  /* check the mode */
    GF2mat prod=H*H.transpose();
    if(prod.is_zero()){
      if(q.debug)
	cout<< "# weakly-self-orthogonal code"<< endl;
      if (q.mode==-1)
	q.mode=1;
    }
    else{
      if(q.mode==-1)
	q.mode=0;
      else{
	cout<<"ERROR: rows of H are not self-orthogonal!"<< endl;
	exit(-1);
      }
    }
  }
  if(q.debug)
    cout<< "# working in mode="<< q.mode<< endl;
  if(q.debug)
    timer.tic(); 
  GF2mat P,G;
  ivec perm;
  int rank=do_gen_matr(H,P,G,perm); 
  if(q.debug&16) // type out matrix 
    cout << "# after permutation H\n"<< H<< endl;
  int wmin=1, rmax=n;
  int w=do_dist_rand(G,P,wmin,rmax);  
  if (q.mode)
    cout << " result: [[" << n<<","<< n-2*rank << ","<< w<< "]] code"<< endl;
  else
    cout << " result: [" << n<<","<< n-rank << ","<< w<< "] code"<< endl;
  if(q.debug){
    cout << "# "; timer.toc_print();
  }
  return 0;
}
#endif /* DEBUG */
