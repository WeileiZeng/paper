#include <stdio.h>
#include <stdlib.h>
#include <itpp/itcomm.h>
#include <sstream>
using namespace std;
using namespace itpp;
#include "mmio.h"
#include "util.h"



/* read Matrix Market sparse integer matrix to GF2mat_sparse */
GF2mat_sparse GF2mat_sparse_mm_read(string fname, bool transpose){
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   

 
    if ((f = fopen(fname.c_str(), "r")) == NULL) 
      ERROR("can't open file %s",fname.c_str());

    if (mm_read_banner(f, &matcode) != 0)
        ERROR("Could not process Matrix Market banner.");

    if (!(mm_is_matrix(matcode) && mm_is_sparse(matcode) && 
	  mm_is_integer(matcode) && mm_is_general(matcode) )){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
	ERROR("input file %s",fname.c_str());
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
      ERROR("Cannot read size in input file %s",fname.c_str());

    GF2mat_sparse sbmat(M, N, nz);

    for (int i = 0; i < nz; i++) {
      long ir,ic,iv; 
      fscanf(f, "%ld %ld %ld\n", &ir , &ic, &iv);
      /* TODO: check if iv=1 */
      sbmat.set_new(ir-1,ic-1, bin(iv));   /* adjust from 1-based to 0-based */
    }
  
    sbmat.compact();

    if (transpose) {
      return sbmat.transpose();
    }
    else {
      return sbmat;
    }
}
