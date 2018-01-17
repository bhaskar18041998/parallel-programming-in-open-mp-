#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<omp.h>
#include<time.h>
void MatMul1( int m, int n, int p, int b, double alpha, double beta, double *A, double *B, double *D )
{   int i,j,k;

    
    for( i = 0 ; i < m ; i++ ){
	
      for( k = 0 ; k < p ; k++ )
      {
        //D[i*p+k] *= beta ;
        for( j = 0 ; j < n ; j++ )
          D[i*p+k] += /*alpha */A[i*n+j]*B[j*p+k] ;
      }
    }


}
void MatMul2( int m, int n, int p, int b, double alpha, double beta, double *A, double *B, double *C )
{  																				 //MATRIX - MATRIX OPERATION..C=alpha*A*B+beta*C
		omp_set_num_threads(2);
	#pragma omp parallel shared(m,n,p,b,alpha,beta,A,B,C) num_threads(4)
	{
		double sum=0;
		int x,y,z,k,i,j;
		 #pragma omp for schedule(dynamic) collapse(2)
			for (  i=0; i<n; i+=b ){
			
		        for ( j=0; j<p; j+=b ){
				
		            for ( k=0; k<m; k+=b ){
					
		            /* normal multiply inside the blocks */
		                for (  x=k; x<k+b; x++ ){
						
		                    for ( y=j; y<j+b; y++ ){
								//sum=C[x*p+y];
							//	if(i==0)
								//	sum=beta*C[x*p+y];
		                        for ( z=i; z<i+b; z++ ){
					
									C[x*p+y] += /*alpha*/A[x*n+z]*B[z*p+y];
								}
						//	if(i==0){
							
						//	#pragma omp critical
							//	C[x*p+y]=sum;
						//	}
						//	else C[x*p+y]=sum;
							}
							
						}
					}
				}
			}
    }
}







//VECTOR MATRIX OPERATION..C=alpha*A*B+beta*C

/*	omp_set_num_threads(4);
	   #pragma omp parallel shared(A,B,C,m,n,p,b,alpha,beta) num_threads(4)  
		{
				int ii,jj,i,j,k,kk;double sum=0;
			  #pragma omp for schedule(dynamic) collapse(3)
			   for (kk = 0; kk < n; kk += b) {
					for (jj = 0; jj < p; jj += b) {
						 for (i = 0; i < m; i++) {
							 for (j = jj; j < ((jj + b) > p ? p : (jj+b)); j++) {
							 	sum = beta*C[i*p+j];
							 	for (k = kk;  k < ((kk + b) > n ? n : (kk + b)); k++) {
							 		sum += alpha*A[i*n+k]*B[k*p+j];
							    }
							 	#pragma omp critical
								 C[i*p+j] = sum;
							}
						}
					}
				}		
		}*/

int main(){
	double *A, *B, *C,*D ,alpha, beta ;int count=0;
    int m, n, p, b,row,col;
    printf("\nenter m,n,p,b\n");
    scanf("%d %d %d %d",&m,&n,&p,&b);
    

  

    A = (double *) malloc( sizeof(double) * m * n ) ;
    B = (double *) malloc( sizeof(double) * n * p ) ;
    C = (double *) malloc( sizeof(double) * m * p ) ;
    D = (double *) malloc( sizeof(double) * m * p ) ;
    
    if ( ( A == NULL ) || ( B == NULL ) || ( C == NULL ) )
    {
        printf( "Out of Memory\n" ) ;
        exit(1) ;
    }
    for( row = 0 ; row < m ; row++ )
        for( col = 0 ; col < p ; col++ )
        	C[row*p+col]=10.236;
    
    for( row = 0 ; row < m ; row++ )
        for( col = 0 ; col < p ; col++ )
        	D[row*p+col]=10.236;

    // m = n = p = 64 ;
    // b = 16 ;
    printf("\nenter 1 mat\n");
	for( row = 0 ; row < m ; row++ )
        for( col = 0 ; col < n ; col++ )
            A[row*n+col]=rand() ;
     printf("\nenter 2 mat\n");
	for( row = 0 ; row < n ; row++ )
        for( col = 0 ; col < p ; col++ )
            B[row*p+col]=rand() ; 

//	printf("\n 1 matrix\n");
//	for( row = 0 ; row < m ; row++ ){
//	    for( col = 0 ; col < n ; col++ )
   //         printf("%lf ",A[row*n+col]) ;
    //        printf("\n");
//	}
//		printf("\n 2 matrix\n");
//	for( row = 0 ; row < n ; row++ ){
//		for( col = 0 ; col < p ; col++ )
    //        printf("%lf  ",B[row*p+col]) ;
//			printf("\n");
//	}

	
	
	
	MatMul1(m,n,p,b,2,3,A,B,D);
	double t1=omp_get_wtime();
	MatMul2(m,n,p,b,2,3,A,B,C);
	double t2=omp_get_wtime();

 
printf("\n%0.20g\t %0.20g\t %0.20g",t1,t2,t2-t1);
	printf("\n");
	

	
/*	printf("\n multiplied mat real case:\n");
	for( row = 0 ; row < m ; row++ ){
		for( col = 0 ; col < p ; col++ )
            printf("%lf  ",D[row*p+col]) ;
			printf("\n");
	}
	
	printf("\n multiplied mat test case:\n");
	for( row = 0 ; row < m ; row++ ){
		for( col = 0 ; col < p ; col++ )
            printf("%lf  ",C[row*p+col]) ;
			printf("\n");
	}*/
		for( row = 0 ; row < m ; row++ ){
		for( col = 0 ; col < p ; col++ ){
		
        	if(C[row*p+col]!=D[row*p+col]){
			printf("   error ");printf("%d,%d",row,col);
			}
			else count++;
			
		}
		
	}
	printf("\ncount: %d ",count);

	

}
