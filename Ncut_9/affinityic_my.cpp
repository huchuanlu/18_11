/*================================================================
* function w = affinityic(emag,ephase,pi,pj,sp_ind,sigma)
* Input:
*   emag = edge strength at each pixel
*   ephase = edge phase at each pixel
*   [pi,pj] = index pair representation for MALTAB sparse matrices
*   sigma = sigma for IC energy
* Output:
*   w = affinity with IC at [pi,pj]
*

% test sequence
f = synimg(10);
[i,j] = cimgnbmap(size(f),2);
[ex,ey,egx,egy] = quadedgep(f);
a = affinityic(ex,ey,egx,egy,i,j)
show_dist_w(f,a);

* Stella X. Yu, Nov 19, 2001.
*=================================================================*/



#include <typeinfo>
# include "mex.h"
# include "math.h"

void mexFunction(
    int nargout,
    mxArray *out[],
    int nargin,
    const mxArray *in[]
)
{
    /* declare variables */
    int nr, nc, np, total,pi_length,pj_length,nsp;
    int  i, j,k, ix, iy, jx, jy, ii, jj, iip1, jjp1, iip2, jjp2, step;
    double sigma, di, dj, a, z, maxori, phase1, phase2, slope;
// 	int *ir, *jc;
    mwIndex*ir,*jc;
	unsigned int *pi, *pj,*sp_ind_middle;
	double *emag, *ephase, *w;
   // printf("Line 43");
    
    /* check argument */
    if (nargin<5) {
        mexErrMsgTxt("Four input arguments required");
    }
    if (nargout>1) {
        mexErrMsgTxt("Too many output arguments");
    }

    /* get edgel information */
	nr = mxGetM(in[0]);
	nc = mxGetN(in[0]);
//     printf("row = %d, col = %d",nr,nc);
	if ( nr*nc ==0 || nr != mxGetM(in[1]) || nc != mxGetN(in[1]) ) {
	    mexErrMsgTxt("Edge magnitude and phase shall be of the same image size");
	}
    emag = mxGetPr(in[0]);
    ephase = mxGetPr(in[1]);
    np = nr * nc;
	pi_length = mxGetM(in[2]);
    pj_length = mxGetM(in[3]);
//     printf("%d\n",pj_length);  //FX
    /* get new index pair */
    if (!mxIsUint32(in[2]) | !mxIsUint32(in[3])) {
        mexErrMsgTxt("Index pair shall be of type UINT32");
    }
  //   printf("Line 68");
	// Change by me
    nsp = mxGetM(in[4]);
  

 //  printf("Line 75");

    pi = (unsigned int*)mxGetData(in[2]);
    pj = (unsigned int*)mxGetData(in[3]);    
    sp_ind_middle = (unsigned int*)mxGetData(in[4]);

  //   printf("Line 81");
    /* create output */
    out[0] = mxCreateSparse(pi_length,1,pi_length,mxREAL);
    if (out[0]==NULL) {
	    mexErrMsgTxt("Not enough memory for the output matrix");
	}
	w = mxGetPr(out[0]);
	ir = mxGetIr(out[0]);
	jc = mxGetJc(out[0]);
	
  //   printf("Line 91");
    /* find my sigma */
	if (nargin<6) {
	    sigma = 0;
    	for (k=0; k<nsp; k++) { 
    	    if (emag[sp_ind_middle[k]]>sigma) { sigma = emag[sp_ind_middle[k]]; }
    	}
    	sigma = sigma / 6;
    	
	} else {
	    sigma = mxGetScalar(in[5]);
	}
    //sigma=0.06;
    //printf("sigma = %6.5f\n",sigma);
	a = 0.5 / (sigma * sigma);
	
//    printf("Line 107");
    /* computation */ 
    total = 0;
    jc[0] = 0;
    for (j=0; j<pi_length; j++) {            

	//printf("%d / %d\n",j,np);        //FX
   //     jc[j] = total;
        jx = sp_ind_middle[pj[j]] / nr; /* col */
        jy = sp_ind_middle[pj[j]] % nr; /* row */
     //printf("Line116");   
      
    //         printf("%d\n",pi[j]); //FX
    //         printf("%d\n",sp_ind_middle[pi[j]]);    
            i = sp_ind_middle[pi[j]];
          
            if (pi[j]==pj[j]) {
                maxori = 1;
            
            } else {
         //printf("Line 125\n");
        // printf("i = %d£¬j = %d",i,j);
                ix = i / nr; 
                iy = i % nr;

                /* scan */            
                di = (double) (iy - jy);
                dj = (double) (ix - jx);
            
                maxori = 0.;
	            phase1 = ephase[sp_ind_middle[pj[j]]];

	               
                /* sample in i direction */
                if (abs(di) >= abs(dj)) {  
            	    slope = dj / di;
            	    step = (iy>=jy) ? 1 : -1;
            	
              	    iip1 = jy;
            	    jjp1 = jx;
	
	               
	                for (ii=0;ii<abs(di);ii++){
	                    iip2 = iip1 + step;
	                    jjp2 = (int)(0.5 + slope*(iip2-jy) + jx);
	  	  
	                    phase2 = ephase[iip2+jjp2*nr];
               
                        //if (1){
	                    if (phase1 != phase2) {
	                        z = (emag[iip1+jjp1*nr] + emag[iip2+jjp2*nr]);
	                        if (z > maxori){
	                            maxori = z;
	                        }
	                    } 
	             
	                    iip1 = iip2;
	                    jjp1 = jjp2;
	                    phase1 = phase2;
	                }
	            
	            /* sample in j direction */    
                } else { 
	                slope = di / dj;
	                step =  (ix>=jx) ? 1: -1;

    	            jjp1 = jx;
	                iip1 = jy;	           
	    
	 
	                for (jj=0;jj<abs(dj);jj++){
	                    jjp2 = jjp1 + step;
	                    iip2 = (int)(0.5+ slope*(jjp2-jx) + jy);
	  	  
	                    phase2 = ephase[iip2+jjp2*nr];
	     
                        //if (1){
	                    if (phase1 != phase2){
	                        z = (emag[iip1+jjp1*nr] + emag[iip2+jjp2*nr]);
	                        if (z > maxori){ 
	                            maxori = z; 
	                        }
	                        
	                    }
	  
	                    iip1 = iip2;
	                    jjp1 = jjp2;
	                    phase1 = phase2;
	                }
                }            
            
                 maxori = 0.5 * maxori;
                 maxori = exp(-maxori * maxori * a);
                
                 //maxori=maxori*a;
                 //maxori=maxori*0.06*255;
                
                 //maxori=exp(-maxori/0.2)*15;
                            
//                 if (maxori>0.4){
//                     maxori=0;
//                 }else{
//                     maxori=1;
//                 }
                

            }       
		    ir[total] = j;
          // printf("Line 212\n");  
            if (maxori<0.01){
                maxori=0.001;
            }
            // printf("Line 216\n");
            w[total] = maxori;
		   //  printf("Line   218\n");
		    total = total + 1;
		
		
    } /* j */
//    printf("Line    223\n");
   jc[1] = total;
    
}  

