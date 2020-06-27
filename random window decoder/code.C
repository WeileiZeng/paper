// Generate some example LDPC codes
#include <itpp/itcomm.h>
#include <stdio.h>
#include "dist.C"
#include "templ.H"
using namespace itpp;
using namespace std;
int main(int argc, char **argv){
  Parser p; 
  p.init(argc,argv);
  p.set_silentmode(false); /* shut down echo */
  if((argc==0)||(strcmp(argv[1],"--help")==0)||(strcmp(argv[1],"-h")==0)){
    cout<< argv[0]<< ": generate random (rw,cw)-regular code\n"
	<< "\tusage: code param=val\n"
        << "\tdebug=1: bitmap for aux information\n"
	<< "\tstring: fout=\"try\": base name for saving\n"
	<< "default: n=28 cw=3 rw=4 ntry=1 dmin=3 wmax=7 seed=0\n"
	<< endl;
    exit (-1);
  }
  int debug=1; p.get(debug,"debug");
  if(debug==0) p.set_silentmode(true); /* shut down echo */
  int n=28, cw=3, rw=4; /* default dim, col and row weights */
  p.get(n,"n");  p.get(cw,"cw");  p.get(rw,"rw");
  int ntry=1, dmin=3;  p.get(ntry,"ntry");  p.get(dmin,"dmin");
  int wmax=7; p.get(wmax,"wmax");
  int seed=0; p.get(seed,"seed");
  if(seed==-1)
    RNG_randomize();
  else
    RNG_reset(seed);

  string fout="try";   p.get(fout,"fout");
  
  cout << "# generating random ("<< cw<<","<<rw<<") code n="<< n<< endl;
  LDPC_Parity_Regular H;      
  LDPC_Parity G; /* dummy */
  GF2mat_sparse_alist alist;

  int dmax=0, d=0, rr=0; 
  for(int i=0; i< ntry; i++){
    H.generate(n, cw, rw,
	       "rand", // random unstructured matrix
	       "200 6"); // optimize girth
    rr=GF2mat(H.get_H()).row_rank();    

    //H.display_stats();
    if (n<=40){
      d=do_dist_clus(wmax,H,G);
      if(d>dmax){
	cout << "# H=\n"<< H.get_H() << endl;
	alist.from_sparse(H.get_H());
	dmax=d;
	cout << "found d="<< d<< endl;
      }
      if(dmax>=dmin)
	break;
    }
    else
      break;
    cout << " rank H="<< rr
	 << " code: ["<< n<< ","<< n-rr<<","<<d<<"]"	 << endl;
  }
  d=dmax;
  //  H.initialize(r,n);
  cout << " saving code: ["<< n<< ","<< n-rr<<","<<d<<"]"	 << endl;
  
  G=LDPC_Parity(alist);
  if(n<=40)
    cout << "# H=\n"<< G.get_H() << endl;
  //LDPC_Code C1(&G);      
  //  C1.save_code(fout+".it");
  //      it_file f(str);
  G.save_alist(fout+".dat");

}
