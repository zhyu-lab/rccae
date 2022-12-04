/* 
 Forward Backward HMM algorithm


 Usage  
 -------
 //forward only
 [alpha, gamma, loglik] = Forward_Backward_Algorithm(PI, A, L, 0) 

 //forward and backward 	
 [alpha, gamma, loglik, beta] = Forward_Backward_Algorithm(PI, A, L, 1)

 //forward and backward, with xi_summed
 [alpha, gamma, loglik, beta, xi_summed] = Forward_Backward_Algorithm(PI, A, L, 2)

 Inputs
 -------

 PI            Initial proabilities (N x 1) : Pr(x_1 = i) , i=1,...,N

 A             State transition probabilities matrix Pr(x_{k} = i| x_{k - 1} = j) such
               sum_{x_k}(A) = 1 => sum(A , 2) = 1, sum of row equals 1

 L        Time indexed Likelihood matrix Pr(z_k | x_k = i) (N x K), i=1,...,N, k=1,...,K. 
               Extracted from B matrix such that B = Pr(z | x) (M x N), sum(B , 1) = 1 and B(z_k , :)' = L(: , k).

 filter        Optional flag. If filter = 0 => Just the alpha probabilities are computed, else (default filter = 1)
               the two passes alpha & beta.


 Ouputs
 -------
 alpha         alpha(i,t) = p(Q(t)=i | y(1:t)) (or p(Q(t)=i, y(1:t)) if scaled=0)
 
 beta          beta(i,t) = p(y(t+1:T) | Q(t)=i)*p(y(t+1:T)|y(1:t)) (or p(y(t+1:T) | Q(t)=i) if scaled=0)

 gamma         gamma(i,t) = p(Q(t)=i | y(1:T))

 loglik        loglik = log p(y(1:T))

 xi_summed     xi_summed(i,j) = sum_{t=}^{T-1} xi(i,j,t)
 
 
 to complie 
 mex -output Forward_Backward_Algorithm.dll Forward_Backward_Algorithm.c

*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

	void ForwardWithScale(int , int , int , double *, double *, double *, 
	double *, double *, double *, double *);
	
	void BackwardWithScale(int , int , double *, double *, double *, 
	double *);
			
	void ComputeGamma(int , int , double *, double *, double *);			
			
	void ComputeXi(int , int , double *, double *, 
	double *, double *, double *, double *, double *);

	void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{
	
	
	double *PI , *A, *L ;
	
	double *alpha, *beta, *gamma, *loglik, *xi_summed; 
	
	const int *dimsPI , *dimsA , *dimsL ;
	
	
	int K , N , numdimsPI , numdimsA  , numdimsL  , filter = 2;
	
	double  *scale, *tmp_xi;
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE INPUT ----------------------------*/	
	/*---------------------------------------------------------------*/
	
	if( (nrhs < 3) | (nrhs >5))
		
	{
		
		mexErrMsgTxt("3 or 4 input are requiered");
		
	}
	
	PI         = mxGetPr(prhs[0]);
    
  numdimsPI  = mxGetNumberOfDimensions(prhs[0]);
    
	dimsPI     = mxGetDimensions(prhs[0]);
	
	if ( (numdimsPI>2) & (dimsPI[1] > dimsPI[0]) )
	{
		
		mexErrMsgTxt("PI must be (N x 1)");			 
		
	}
	
	
	A         = mxGetPr(prhs[1]);
    
  numdimsA  = mxGetNumberOfDimensions(prhs[1]);
    
	dimsA     = mxGetDimensions(prhs[1]);
	
	if ( (numdimsA>2) & (dimsA[1] != dimsA[0]) )
	{
		
		mexErrMsgTxt("A must be (N x N)");			 
		
	}
	
	
	L         = mxGetPr(prhs[2]);
    
  numdimsL  = mxGetNumberOfDimensions(prhs[2]);
    
	dimsL     = mxGetDimensions(prhs[2]);
	
	if ( (numdimsL>2) & (dimsL[0] != dimsA[0]) )
	{
		
		mexErrMsgTxt("L must be (N x K)");			 
		
	}
	
  N         = dimsL[0];
	
	K         = dimsL[1];
		
	
	if (nrhs == 4)
		
	{
		filter = (int)mxGetScalar(prhs[3]);
		
	}
	
	
	/*---------------------------------------------------------------*/
	/*---------------------- PARSE OUTPUT ---------------------------*/	
	/*---------------------------------------------------------------*/
	
	plhs[0]       = mxCreateDoubleMatrix(N , K , mxREAL);
	
	alpha         = mxGetPr(plhs[0]);
	
	
	plhs[1]       = mxCreateDoubleMatrix(N , K , mxREAL);
    
	gamma         = mxGetPr(plhs[1]);
	
	
	plhs[2]       = mxCreateDoubleMatrix(1 , 1 , mxREAL);
    
	loglik        = mxGetPr(plhs[2]);
	
	
	if (filter > 0) { /*calculate beta*/
				plhs[3]       = mxCreateDoubleMatrix(N , K , mxREAL);
	    
				beta          = mxGetPr(plhs[3]);
	}
	

	
	if (filter > 1) { /*calculate xi_summed*/
				plhs[4]       = mxCreateDoubleMatrix(N , N , mxREAL);
			    
				xi_summed     = mxGetPr(plhs[4]);
	}
	

	
	/*---------------------------------------------------------------*/
	/*--------- Internal Tempory vector & matrices ------------------*/
	/*---------------------------------------------------------------*/

	scale         = (double *)mxMalloc(K*sizeof(double));
	
	tmp_xi				= (double *)mxMalloc(N*N*sizeof(double));
	/*---------------------------------------------------------------*/
	/*------------------------ MAIN CALL ----------------------------*/
	/*---------------------------------------------------------------*/
	
	ForwardWithScale(N, K, filter, PI, A, L, alpha, gamma, loglik, scale);
	
	if (filter > 0) /*calculate both alpha and beta, but xi_summed is not calculated*/
			BackwardWithScale(N, K, PI, A, L, beta);
			
	if (filter == 1) /*calculate gamma with alpha and beta*/
			ComputeGamma(N, K, alpha, beta, gamma);			
			
	if (filter > 1) /*calculate gamma and xi_summed based on xi*/
			ComputeXi(N, K, A, L, alpha, beta, gamma, tmp_xi, xi_summed);
	
	/*---------------------------------------------------------------*/
	/*------------------------ FREE MEMORY --------------------------*/
	/*---------------------------------------------------------------*/
	
	
 	mxFree(scale);
	
	mxFree(tmp_xi);
	
}

	
void ForwardWithScale(int N, int K, int filter, double *PI, double *A, double *L, 
	double *alpha, double *gamma, double *loglik, double *scale)
{
	int	i, j; 	/* state indices */
	int	t;	/* time index */

	double sum;	/* partial sum */

	/* 1. Initialization */

	scale[0] = 0.0;	
	for (i = 0; i < N; i++) {
		alpha[i] = PI[i]*L[i];
		scale[0] += alpha[i];
	}
	for (i = 0; i <N; i++) 
		alpha[i] /= scale[0]; 
	
	/* 2. Induction */

	for (t = 1; t < K; t++) {
		scale[t] = 0.0;
		for (j = 0; j <N; j++) {
			sum = 0.0;
			for (i = 0; i <N; i++)
				sum += alpha[i+(t-1)*N]*A[i+j*N]; 
				/*sum += alpha[t][i]* (phmm->A[i][j]); 
				sum   += At[j + iN]*alpha[j + kN1];*/
			
			alpha[j+t*N] = sum*L[j+t*N];	
			/*alpha[t+1][j] = sum*(phmm->B[j][O[t+1]]);*/
			scale[t] += alpha[j+t*N];
			/*scale[t+1] += alpha[t+1][j];*/
		}
		for (j = 0; j <N; j++) 
			alpha[j+t*N] /= scale[t];
			/*alpha[t+1][j] /= scale[t+1]; */
	}

	/* 3. Termination */
	*loglik = 0.0;

	for (t = 0; t < K; t++)
		*loglik += log(scale[t]);
		
	if (filter == 0) {/*forward only*/
				for (i = 0; i<N; i++) 
							for (j = 0; j<K; j++) 
										gamma[i+j*N] = alpha[i+j*N];
	}
}

void BackwardWithScale(int N, int K, double *PI, double *A, double *L, 
double *beta)	
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
        double sum;
 
 
        /* 1. Initialization */
        
        for (i = 0; i <N; i++) {
                beta[i+(K-1)*N] = 1.0;
                /*beta[T][i] = 1.0;*/
        }
        
        
        /* 2. Induction */
 
        for (t = K - 2; t >= 0; t--) {
        				sum = 0.0;
                for (i = 0; i <N; i++) {
                        beta[i+t*N] = 0.0;
                        for (j = 0; j <N; j++)
                                beta[i+t*N] += A[i+j*N]*L[j+(t+1)*N]*beta[j+(t+1)*N];
                        sum += beta[i+t*N];
                }
                for (i = 0; i <N; i++) 
                				beta[i+t*N] /= sum;
                
					}
        	
}

void ComputeGamma(int N, int K, double *alpha, double *beta, double *gamma)	
/*calculate gamma(:,T) by definition, normalise(alpha(:,T).*beta(:,T))*/
{
        int     i, j;   /* state indices */
        int     t;      /* time index */
        double sum;

	for (t = 0; t<K; t++) {
				  sum = 0.0;
        	for (i = 0; i<N; i++) {
             			gamma[i+t*N] = alpha[i+t*N]*beta[i+t*N];
             			sum += gamma[i+t*N];
       		}
        	for (i = 0; i <N; i++) 
								gamma[i+t*N] /= sum;
	} 	
		
}


void ComputeXi(int N, int K, double *A, double *L, double *alpha, double *beta, double *gamma, double *tmp_xi, double *xi_summed)	
{
	int i, j;
	int t;
	double sum;

	for (i = 0; i <N; i++) 
			for (j = 0; j <N; j++) 
				xi_summed[i+j*N] = 0.0;
				
	for (t = 0; t <= K - 2; t++) {
		sum = 0.0;	
		for (i = 0; i <N; i++) 
				for (j = 0; j <N; j++) {
					tmp_xi[i+j*N] = alpha[i+t*N]*A[i+j*N]*L[j+(t+1)*N]*beta[j+(t+1)*N];
					sum += tmp_xi[i+j*N];
				}
		/*calculate gamma with xi for 1 to T-1 (can not calculate gamma(:,T) with this formula) */
		for (i = 0; i <N; i++) {
				gamma[i+t*N] = 0.0;
				for (j = 0; j <N; j++) { 
						tmp_xi[i+j*N]  /= sum;
						xi_summed[i+j*N] += tmp_xi[i+j*N];
						gamma[i+t*N] += tmp_xi[i+j*N];
				}
		}
	}
	
	/*finally calculate gamma(:,T) by definition, normalise(alpha(:,T).*beta(:,T)) */
	sum = 0.0;
  for (i = 0; i<N; i++) {
   			gamma[i+(K-1)*N] = alpha[i+(K-1)*N]*beta[i+(K-1)*N];
    		sum += gamma[i+(K-1)*N];
  }
  for (i = 0; i <N; i++) 
				gamma[i+(K-1)*N] /= sum;
	
}
