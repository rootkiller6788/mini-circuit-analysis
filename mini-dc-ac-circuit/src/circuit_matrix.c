#include "circuit_matrix.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* Matrix Allocation L3 */
double* matrix_alloc(int rows,int cols){if(rows<=0||cols<=0)return NULL;return(double*)calloc((size_t)(rows*cols),sizeof(double));}
void matrix_free(double *m){free(m);}
void matrix_zero(double *A,int rows,int cols){if(!A)return;for(int i=0;i<rows*cols;i++)A[i]=0.0;}
void matrix_identity(double *A,int n){if(!A)return;for(int i=0;i<n*n;i++)A[i]=0.0;for(int i=0;i<n;i++)A[i*n+i]=1.0;}
void matrix_copy(const double *A,double *B,int rows,int cols){if(!A||!B)return;for(int i=0;i<rows*cols;i++)B[i]=A[i];}
void matrix_add(const double *A,const double *B,double *C,int rows,int cols){if(!A||!B||!C)return;for(int i=0;i<rows*cols;i++)C[i]=A[i]+B[i];}
void matrix_sub(const double *A,const double *B,double *C,int rows,int cols){if(!A||!B||!C)return;for(int i=0;i<rows*cols;i++)C[i]=A[i]-B[i];}

/* Matrix Multiply L3 O(m*k*n) */
void matrix_mul(const double *A,const double *B,double *C,int m,int k,int n){
  if(!A||!B||!C)return;
  for(int i=0;i<m;i++){for(int j=0;j<n;j++){double sum=0;for(int p=0;p<k;p++)sum+=A[i*k+p]*B[p*n+j];C[i*n+j]=sum;}}}
void matrix_vec_mul(const double *A,const double *x,double *y,int m,int n){
  if(!A||!x||!y)return;
  for(int i=0;i<m;i++){double sum=0;for(int j=0;j<n;j++)sum+=A[i*n+j]*x[j];y[i]=sum;}}
void matrix_transpose(const double *A,double *B,int rows,int cols){
  if(!A||!B)return;for(int i=0;i<rows;i++)for(int j=0;j<cols;j++)B[j*rows+i]=A[i*cols+j];}

/* Determinant via LU L3 */
int matrix_determinant(const double *A,int n,double *det){
  if(!A||!det||n<=0)return -1;
  double *LU=(double*)malloc((size_t)(n*n)*sizeof(double));if(!LU)return -1;
  for(int i=0;i<n*n;i++)LU[i]=A[i];int *pivot=(int*)malloc((size_t)n*sizeof(int));if(!pivot){free(LU);return -1;}
  int swaps=matrix_lu_decompose(LU,n,pivot);if(swaps<0){free(LU);free(pivot);return -1;}
  *det=1.0;for(int i=0;i<n;i++)*det*=LU[i*n+i];if(swaps%2)*det=-*det;
  free(LU);free(pivot);return 0;}

/* LU Decomposition L3 Doolittle */
int matrix_lu_decompose(double *A,int n,int *pivot){
  if(!A||!pivot||n<=0)return -1;
  for(int i=0;i<n;i++)pivot[i]=i;int swaps=0;
  for(int col=0;col<n;col++){
    int max_row=col;double max_val=fabs(A[col*n+col]);
    for(int row=col+1;row<n;row++){double v=fabs(A[row*n+col]);if(v>max_val){max_val=v;max_row=row;}}
    if(max_val<1e-15)return -1;
    if(max_row!=col){
      for(int j=0;j<n;j++){double t=A[col*n+j];A[col*n+j]=A[max_row*n+j];A[max_row*n+j]=t;}
      int tp=pivot[col];pivot[col]=pivot[max_row];pivot[max_row]=tp;swaps++;}
    double piv=A[col*n+col];
    for(int row=col+1;row<n;row++){A[row*n+col]/=piv;for(int j=col+1;j<n;j++)A[row*n+j]-=A[row*n+col]*A[col*n+j];}}
  return swaps;}

/* LU Solve L3 */
void matrix_lu_solve(const double *LU,const int *pivot,const double *b,double *x,int n){
  if(!LU||!pivot||!b||!x)return;
  double *y=(double*)malloc((size_t)n*sizeof(double));if(!y)return;
  for(int i=0;i<n;i++)y[i]=b[pivot[i]];
  for(int i=0;i<n;i++){for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}
  for(int i=n-1;i>=0;i--){double sum=y[i];for(int j=i+1;j<n;j++)sum-=LU[i*n+j]*x[j];x[i]=sum/LU[i*n+i];}
  free(y);}

/* Full Solve L3 */
int matrix_solve(double *A,const double *b,double *x,int n){
  if(!A||!b||!x||n<=0)return -1;int *pivot=(int*)malloc((size_t)n*sizeof(int));if(!pivot)return -1;
  int r=matrix_lu_decompose(A,n,pivot);if(r<0){free(pivot);return -1;}matrix_lu_solve(A,pivot,b,x,n);free(pivot);return 0;}

/* Matrix Inverse L3 */
int matrix_inverse(const double *A,double *Ainv,int n){
  if(!A||!Ainv||n<=0)return -1;
  double *LU=(double*)malloc((size_t)(n*n)*sizeof(double));if(!LU)return -1;
  for(int i=0;i<n*n;i++)LU[i]=A[i];int *pivot=(int*)malloc((size_t)n*sizeof(int));if(!pivot){free(LU);return -1;}
  int r=matrix_lu_decompose(LU,n,pivot);if(r<0){free(LU);free(pivot);return -1;}
  double *e=(double*)calloc((size_t)n,sizeof(double));double *col=(double*)malloc((size_t)n*sizeof(double));
  if(!e||!col){free(e);free(col);free(LU);free(pivot);return -1;}
  for(int j=0;j<n;j++){for(int i=0;i<n;i++)e[i]=0.0;e[j]=1.0;matrix_lu_solve(LU,pivot,e,col,n);for(int i=0;i<n;i++)Ainv[i*n+j]=col[i];}
  free(e);free(col);free(LU);free(pivot);return 0;}

/* Matrix Norms L3 */
double matrix_cond_est(const double *A,int n){
  double norm1=matrix_norm_1(A,n,n);if(norm1<=0)return INFINITY;
  double *Ainv=(double*)malloc((size_t)(n*n)*sizeof(double));if(!Ainv)return INFINITY;
  int r=matrix_inverse(A,Ainv,n);if(r<0){free(Ainv);return INFINITY;}
  double norm1_inv=matrix_norm_1(Ainv,n,n);free(Ainv);return norm1*norm1_inv;}
double matrix_norm_frobenius(const double *A,int rows,int cols){
  if(!A)return 0;double sum=0;for(int i=0;i<rows*cols;i++)sum+=A[i]*A[i];return sqrt(sum);}
double matrix_norm_1(const double *A,int rows,int cols){
  if(!A)return 0;double max_sum=0;
  for(int j=0;j<cols;j++){double col_sum=0;for(int i=0;i<rows;i++)col_sum+=fabs(A[i*cols+j]);if(col_sum>max_sum)max_sum=col_sum;}
  return max_sum;}
double matrix_norm_inf(const double *A,int rows,int cols){
  if(!A)return 0;double max_sum=0;
  for(int i=0;i<rows;i++){double row_sum=0;for(int j=0;j<cols;j++)row_sum+=fabs(A[i*cols+j]);if(row_sum>max_sum)max_sum=row_sum;}
  return max_sum;}

/* Eigenvalues 2x2 L3 */
int matrix_eigen_2x2(double a,double b,double c,double d,double *lr1,double *li1,double *lr2,double *li2){
  double trace=a+d;double det=a*d-b*c;double disc=trace*trace-4*det;
  if(disc>=0){if(lr1)*lr1=(trace+sqrt(disc))/2;if(li1)*li1=0;if(lr2)*lr2=(trace-sqrt(disc))/2;if(li2)*li2=0;return 1;}
  else{if(lr1)*lr1=trace/2;if(li1)*li1=sqrt(-disc)/2;if(lr2)*lr2=trace/2;if(li2)*li2=-sqrt(-disc)/2;return 0;}}

/* Eigenvalues Symmetric 3x3 L3 */
int matrix_eigen_sym_3x3(const double *A,double *eval){
  if(!A||!eval)return 0;
  double a=A[0],b=A[1],c=A[2],d=A[4],e=A[5],f=A[8];
  double p1=b*b+c*c+e*e;if(p1<1e-15){eval[0]=a;eval[1]=d;eval[2]=f;return 3;}
  double q=(a+d+f)/3.0;double p2=(a-q)*(a-q)+(d-q)*(d-q)+(f-q)*(f-q)+2*p1;
  double p=sqrt(p2/6.0);double r_num=a*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d)-q*(a*d+a*f+d*f-b*b-c*c-e*e)+2*q*q*q;
  double r=r_num/(2*p*p*p);if(r>1.0)r=1.0;if(r<-1.0)r=-1.0;
  double phi=acos(r)/3.0;eval[0]=q+2*p*cos(phi);eval[1]=q+2*p*cos(phi+2*M_PI/3);eval[2]=q+2*p*cos(phi+4*M_PI/3);return 3;}

/* Complex Matrix Operations L3 */
double complex* cmatrix_alloc(int rows,int cols){if(rows<=0||cols<=0)return NULL;return(double complex*)calloc((size_t)(rows*cols),sizeof(double complex));}
void cmatrix_free(double complex *m){free(m);}
void cmatrix_zero(double complex *A,int rows,int cols){if(!A)return;for(int i=0;i<rows*cols;i++)A[i]=0;}
void cmatrix_add(const double complex *A,const double complex *B,double complex *C,int rows,int cols){
  if(!A||!B||!C)return;for(int i=0;i<rows*cols;i++)C[i]=A[i]+B[i];}
void cmatrix_mul(const double complex *A,const double complex *B,double complex *C,int m,int k,int n){
  if(!A||!B||!C)return;
  for(int i=0;i<m;i++){for(int j=0;j<n;j++){double complex sum=0;for(int p=0;p<k;p++)sum+=A[i*k+p]*B[p*n+j];C[i*n+j]=sum;}}}
void cmatrix_vec_mul(const double complex *A,const double complex *x,double complex *y,int m,int n){
  if(!A||!x||!y)return;
  for(int i=0;i<m;i++){double complex sum=0;for(int j=0;j<n;j++)sum+=A[i*n+j]*x[j];y[i]=sum;}}

/* Complex LU L3 */
int cmatrix_lu_decompose(double complex *A,int n,int *pivot){
  if(!A||!pivot||n<=0)return -1;
  for(int i=0;i<n;i++)pivot[i]=i;int swaps=0;
  for(int col=0;col<n;col++){
    int max_row=col;double max_val=cabs(A[col*n+col]);
    for(int row=col+1;row<n;row++){double v=cabs(A[row*n+col]);if(v>max_val){max_val=v;max_row=row;}}
    if(max_val<1e-15)return -1;
    if(max_row!=col){for(int j=0;j<n;j++){double complex t=A[col*n+j];A[col*n+j]=A[max_row*n+j];A[max_row*n+j]=t;}int tp=pivot[col];pivot[col]=pivot[max_row];pivot[max_row]=tp;swaps++;}
    double complex piv=A[col*n+col];
    for(int row=col+1;row<n;row++){A[row*n+col]/=piv;for(int j=col+1;j<n;j++)A[row*n+j]-=A[row*n+col]*A[col*n+j];}}
  return swaps;}
void cmatrix_lu_solve(const double complex *LU,const int *pivot,const double complex *b,double complex *x,int n){
  if(!LU||!pivot||!b||!x)return;
  double complex *y=(double complex*)malloc((size_t)n*sizeof(double complex));if(!y)return;
  for(int i=0;i<n;i++)y[i]=b[pivot[i]];
  for(int i=0;i<n;i++){for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}
  for(int i=n-1;i>=0;i--){double complex sum=y[i];for(int j=i+1;j<n;j++)sum-=LU[i*n+j]*x[j];x[i]=sum/LU[i*n+i];}
  free(y);}
int cmatrix_solve(double complex *A,const double complex *b,double complex *x,int n){
  if(!A||!b||!x||n<=0)return -1;int *p=(int*)malloc((size_t)n*sizeof(int));if(!p)return -1;
  int r=cmatrix_lu_decompose(A,n,p);if(r<0){free(p);return -1;}cmatrix_lu_solve(A,p,b,x,n);free(p);return 0;}

/* Sparse Matrix CRS Format L8 */
void sparse_init(SparseMatrix_t *sp,int rows,int cols,int capacity){
  if(!sp)return;sp->n_rows=rows;sp->n_cols=cols;sp->nnz=0;sp->max_nnz=capacity;
  sp->values=(double*)malloc((size_t)capacity*sizeof(double));
  sp->col_idx=(int*)malloc((size_t)capacity*sizeof(int));
  sp->row_ptr=(int*)calloc((size_t)(rows+1),sizeof(int));}
void sparse_free(SparseMatrix_t *sp){if(!sp)return;free(sp->values);free(sp->col_idx);free(sp->row_ptr);sp->nnz=0;}
void sparse_insert(SparseMatrix_t *sp,int i,int j,double value){
  if(!sp||i<0||i>=sp->n_rows||j<0||j>=sp->n_cols||value==0.0)return;
  for(int k=sp->row_ptr[i];k<sp->nnz&&(sp->row_ptr[i]<=k);k++){
    if(sp->col_idx[k]==j&&k>=sp->row_ptr[i]&&(k<sp->row_ptr[i+1]||(i==sp->n_rows-1&&k<sp->nnz))){sp->values[k]=value;return;}}
  if(sp->nnz>=sp->max_nnz)return;int pos=sp->nnz;
  for(int k=sp->row_ptr[i];k<sp->nnz;k++){if(sp->col_idx[k]>j){pos=k;break;}}
  for(int k=sp->nnz;k>pos;k--){sp->values[k]=sp->values[k-1];sp->col_idx[k]=sp->col_idx[k-1];}
  sp->values[pos]=value;sp->col_idx[pos]=j;sp->nnz++;
  for(int r=i+1;r<=sp->n_rows;r++)sp->row_ptr[r]++;}
double sparse_get(const SparseMatrix_t *sp,int i,int j){
  if(!sp||i<0||i>=sp->n_rows||j<0||j>=sp->n_cols)return 0.0;
  int start=sp->row_ptr[i];int end=(i<sp->n_rows-1)?sp->row_ptr[i+1]:sp->nnz;
  for(int k=start;k<end;k++){if(sp->col_idx[k]==j)return sp->values[k];}return 0.0;}
void sparse_vec_mul(const SparseMatrix_t *A,const double *x,double *y){
  if(!A||!x||!y)return;
  for(int i=0;i<A->n_rows;i++){double sum=0;int start=A->row_ptr[i];int end=(i<A->n_rows-1)?A->row_ptr[i+1]:A->nnz;
    for(int k=start;k<end;k++)sum+=A->values[k]*x[A->col_idx[k]];y[i]=sum;}}
int sparse_solve(SparseMatrix_t *A,const double *b,double *x){
  if(!A||!b||!x)return -1;int n=A->n_rows;double *full=(double*)calloc((size_t)(n*n),sizeof(double));if(!full)return -1;
  for(int i=0;i<n;i++){int s=A->row_ptr[i];int e=(i<n-1)?A->row_ptr[i+1]:A->nnz;for(int k=s;k<e;k++)full[i*n+A->col_idx[k]]=A->values[k];}
  int r=matrix_solve(full,b,x,n);free(full);return r;}
int sparse_build_from_circuit(const Circuit_t *ckt,SparseMatrix_t *G,double *I_vec){
  if(!ckt||!G||!I_vec)return -1;int n=ckt->n_nodes-1;if(n<=0)return -1;
  for(int i=0;i<n;i++)I_vec[i]=0;
  int *nm=(int*)malloc((size_t)ckt->n_nodes*sizeof(int));if(!nm)return -1;
  int mi=0;for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)nm[i]=-1;else nm[i]=mi++;}
  for(int k=0;k<ckt->n_resistors;k++){const Resistor_t *r=&ckt->resistors[k];int ni=r->t[0].node_id,nj=r->t[1].node_id;
    if(r->R<=0.0)continue;double g=1.0/r->R;
    if(nm[ni]>=0)sparse_insert(G,nm[ni],nm[ni],sparse_get(G,nm[ni],nm[ni])+g);
    if(nm[nj]>=0)sparse_insert(G,nm[nj],nm[nj],sparse_get(G,nm[nj],nm[nj])+g);
    if(nm[ni]>=0&&nm[nj]>=0){sparse_insert(G,nm[ni],nm[nj],sparse_get(G,nm[ni],nm[nj])-g);sparse_insert(G,nm[nj],nm[ni],sparse_get(G,nm[nj],nm[ni])-g);}}
  for(int k=0;k<ckt->n_dc_isrc;k++){const DCCurrentSource_t *s=&ckt->dc_isrc[k];
    if(nm[s->t[0].node_id]>=0)I_vec[nm[s->t[0].node_id]]+=s->I_dc;
    if(nm[s->t[1].node_id]>=0)I_vec[nm[s->t[1].node_id]]-=s->I_dc;}
  free(nm);return n;}

/* Matrix Trace and Special Operations L3 */
double matrix_trace(const double *A,int n){if(!A||n<=0)return 0;double t=0;for(int i=0;i<n;i++)t+=A[i*n+i];return t;}
void matrix_scale(double *A,double scalar,int rows,int cols){if(!A)return;for(int i=0;i<rows*cols;i++)A[i]*=scalar;}
void matrix_negate(double *A,int rows,int cols){matrix_scale(A,-1.0,rows,cols);}

/* Special Matrices for Circuit Analysis L3 */
void incidence_matrix_build(const Circuit_t *ckt,double *A,int *rows,int *cols){
  if(!ckt||!A||!rows||!cols)return;*rows=ckt->n_nodes;*cols=ckt->n_branches;
  for(int i=0;i<(*rows)*(*cols);i++)A[i]=0.0;
  for(int k=0;k<ckt->n_branches;k++){const Branch_t *b=&ckt->branches[k];
    if(b->node_from>=0&&b->node_from<ckt->n_nodes)A[b->node_from*(*cols)+k]=1.0;
    if(b->node_to>=0&&b->node_to<ckt->n_nodes)A[b->node_to*(*cols)+k]=-1.0;}}

/* Reduced Incidence Matrix (no ground) L3 */
void reduced_incidence_matrix(const Circuit_t *ckt,double *A,int *rows,int *cols){
  if(!ckt||!A||!rows||!cols)return;*rows=ckt->n_nodes-1;*cols=ckt->n_branches;
  for(int i=0;i<(*rows)*(*cols);i++)A[i]=0.0;int r=0;
  for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)continue;
    for(int k=0;k<ckt->n_branches;k++){const Branch_t *b=&ckt->branches[k];
      if(b->node_from==i)A[r*(*cols)+k]=1.0;
      else if(b->node_to==i)A[r*(*cols)+k]=-1.0;}r++;}}
double matrix_symmetry_check(const double *A,int n){
  if(!A)return -1;double max_err=0;
  for(int i=0;i<n;i++)for(int j=0;j<i;j++){double err=fabs(A[i*n+j]-A[j*n+i]);if(err>max_err)max_err=err;}return max_err;}
int matrix_is_positive_definite(const double *A,int n){
  if(!A||n<=0)return 0;
  double *L=(double*)calloc((size_t)(n*n),sizeof(double));if(!L)return 0;
  for(int j=0;j<n;j++){double sum=0;for(int k=0;k<j;k++)sum+=L[j*n+k]*L[j*n+k];
    L[j*n+j]=sqrt(A[j*n+j]-sum);if(A[j*n+j]-sum<=0){free(L);return 0;}
    for(int i=j+1;i<n;i++){sum=0;for(int k=0;k<j;k++)sum+=L[i*n+k]*L[j*n+k];L[i*n+j]=(A[i*n+j]-sum)/L[j*n+j];}}
  free(L);return 1;}

/* Matrix Exponential via Pade Approximation L8 */
void matrix_exp_pade(const double *A,int n,double *expA){
  if(!A||!expA||n<=0)return;
  double *I=(double*)calloc((size_t)(n*n),sizeof(double));
  double *A2=(double*)calloc((size_t)(n*n),sizeof(double));
  double *num=(double*)calloc((size_t)(n*n),sizeof(double));
  double *den=(double*)calloc((size_t)(n*n),sizeof(double));
  if(!I||!A2||!num||!den){free(I);free(A2);free(num);free(den);return;}
  for(int i=0;i<n;i++)I[i*n+i]=1.0;
  /* Pade (1,1): exp(A) ~ (I + A/2) / (I - A/2) */
  for(int i=0;i<n*n;i++){num[i]=I[i]+0.5*A[i];den[i]=I[i]-0.5*A[i];}
  int *pivot=(int*)malloc((size_t)n*sizeof(int));
  for(int j=0;j<n;j++){
    double *rhs=(double*)calloc((size_t)n,sizeof(double));
    double *col=(double*)malloc((size_t)n*sizeof(double));
    for(int i=0;i<n;i++)rhs[i]=num[i*n+j];
    double *D=(double*)malloc((size_t)(n*n)*sizeof(double));
    for(int i=0;i<n*n;i++)D[i]=den[i];
    matrix_lu_decompose(D,n,pivot);matrix_lu_solve(D,pivot,rhs,col,n);
    for(int i=0;i<n;i++)expA[i*n+j]=col[i];free(rhs);free(col);free(D);}
  free(I);free(A2);free(num);free(den);free(pivot);}

/* Cholesky Decomposition for SPD matrices L5 */
int matrix_cholesky(const double *A,int n,double *L){
  if(!A||!L||n<=0)return -1;
  for(int i=0;i<n*n;i++)L[i]=0.0;
  for(int j=0;j<n;j++){double sum=0;for(int k=0;k<j;k++)sum+=L[j*n+k]*L[j*n+k];
    double diag=A[j*n+j]-sum;if(diag<=0)return -1;L[j*n+j]=sqrt(diag);
    for(int i=j+1;i<n;i++){sum=0;for(int k=0;k<j;k++)sum+=L[i*n+k]*L[j*n+k];L[i*n+j]=(A[i*n+j]-sum)/L[j*n+j];}}
  return 0;}

/* QR Decomposition via Gram-Schmidt L5 */
int matrix_qr_decompose(const double *A,int m,int n,double *Q,double *R){
  if(!A||!Q||!R||m<=0||n<=0)return -1;
  for(int i=0;i<m*n;i++){Q[i]=0;R[i]=0;}
  double *v=(double*)malloc((size_t)m*sizeof(double));if(!v)return -1;
  for(int k=0;k<n;k++){
    for(int i=0;i<m;i++)v[i]=A[i*n+k];
    for(int j=0;j<k;j++){double dot=0;for(int i=0;i<m;i++)dot+=Q[i*n+j]*A[i*n+k];R[j*n+k]=dot;
      for(int i=0;i<m;i++)v[i]-=dot*Q[i*n+j];}
    double norm=0;for(int i=0;i<m;i++)norm+=v[i]*v[i];R[k*n+k]=sqrt(norm);
    if(R[k*n+k]>1e-15){for(int i=0;i<m;i++)Q[i*n+k]=v[i]/R[k*n+k];}}
  free(v);return 0;}

/* Givens Rotation for sparse elimination L8 */
void givens_rotation(double a,double b,double *c,double *s,double *r){
  if(b==0){*c=1;*s=0;*r=a;return;}
  if(fabs(b)>fabs(a)){double tau=-a/b;*s=1.0/sqrt(1+tau*tau);*c=*s*tau;}
  else{double tau=-b/a;*c=1.0/sqrt(1+tau*tau);*s=*c*tau;}
  *r=(*c)*a-(*s)*b;}
void givens_apply_row(double *A,int n,int row_i,int row_j,double c,double s){
  for(int k=0;k<n;k++){double ai=A[row_i*n+k];double aj=A[row_j*n+k];
    A[row_i*n+k]=c*ai-s*aj;A[row_j*n+k]=s*ai+c*aj;}}

/* Householder Reflection L8 */
double householder_vector(const double *x,double *v,int n){
  if(!x||!v||n<=0)return 0;
  double norm=0;for(int i=0;i<n;i++)norm+=x[i]*x[i];norm=sqrt(norm);if(norm<1e-15)return 0;
  double sigma=0;for(int i=1;i<n;i++)sigma+=x[i]*x[i];v[0]=x[0]+(x[0]>=0?1:-1)*norm;
  for(int i=1;i<n;i++)v[i]=x[i];double beta=2.0/(sigma+v[0]*v[0]);return beta;}
void householder_apply(double *A,int m,int n,double *v,double beta,int col){
  if(!A||!v)return;
  for(int j=0;j<n;j++){double dot=0;for(int i=col;i<m;i++)dot+=v[i-col]*A[i*n+j];
    for(int i=col;i<m;i++)A[i*n+j]-=beta*dot*v[i-col];}}

/* Matrix Rank via Gaussian Elimination L3 */
int matrix_rank(const double *A,int m,int n){
  if(!A||m<=0||n<=0)return 0;
  int max_rc=(m<n)?m:n;double *B=(double*)malloc((size_t)(m*n)*sizeof(double));if(!B)return 0;
  for(int i=0;i<m*n;i++)B[i]=A[i];int rank=0;
  for(int col=0;col<n&&rank<m;col++){
    int pivot_row=-1;for(int row=rank;row<m;row++){if(fabs(B[row*n+col])>1e-12){pivot_row=row;break;}}
    if(pivot_row<0)continue;
    if(pivot_row!=rank){for(int j=col;j<n;j++){double t=B[rank*n+j];B[rank*n+j]=B[pivot_row*n+j];B[pivot_row*n+j]=t;}}
    double pivot=B[rank*n+col];for(int row=rank+1;row<m;row++){double factor=B[row*n+col]/pivot;for(int j=col;j<n;j++)B[row*n+j]-=factor*B[rank*n+j];}
    rank++;}
  free(B);return rank;}

/* Matrix Null Space Dimension */
int matrix_nullity(const double *A,int m,int n){int r=matrix_rank(A,m,n);return n-r;}

/* Least Squares Solve: min ||Ax-b||_2 L5 */
int matrix_least_squares(const double *A,int m,int n,const double *b,double *x){
  if(!A||!b||!x||m<=0||n<=0||m<n)return -1;
  double *AtA=(double*)calloc((size_t)(n*n),sizeof(double));
  double *Atb=(double*)calloc((size_t)n,sizeof(double));
  if(!AtA||!Atb){free(AtA);free(Atb);return -1;}
  for(int i=0;i<n;i++){for(int j=0;j<n;j++){double s=0;for(int k=0;k<m;k++)s+=A[k*n+i]*A[k*n+j];AtA[i*n+j]=s;}}
  for(int i=0;i<n;i++){double s=0;for(int k=0;k<m;k++)s+=A[k*n+i]*b[k];Atb[i]=s;}
  int ret=matrix_solve(AtA,Atb,x,n);free(AtA);free(Atb);return ret;}

/* Pseudoinverse via SVD approximation for well-conditioned matrices L8 */
int matrix_pseudoinverse(const double *A,int m,int n,double *Apinv){
  if(!A||!Apinv||m<=0||n<=0)return -1;
  if(m>=n){double *AtA=(double*)calloc((size_t)(n*n),sizeof(double));
    double *At=(double*)malloc((size_t)(n*m)*sizeof(double));if(!AtA||!At){free(AtA);free(At);return -1;}
    matrix_transpose(A,At,m,n);matrix_mul(At,A,AtA,n,m,n);
    int ret=matrix_inverse(AtA,Apinv,n);free(AtA);free(At);if(ret<0)return -1;
    double *temp=(double*)malloc((size_t)(n*m)*sizeof(double));matrix_mul(Apinv,At,temp,n,n,m);
    for(int i=0;i<n*m;i++)Apinv[i]=temp[i];free(temp);return 0;}
  else{return -1;}}

/* Power Iteration for dominant eigenvalue L8 */
double matrix_power_iteration(const double *A,int n,double *eigenvector,int max_iter,double tol){
  if(!A||!eigenvector||n<=0)return 0;
  for(int i=0;i<n;i++)eigenvector[i]=1.0/sqrt((double)n);double lambda=0;
  for(int iter=0;iter<max_iter;iter++){double *y=(double*)calloc((size_t)n,sizeof(double));
    for(int i=0;i<n;i++){for(int j=0;j<n;j++)y[i]+=A[i*n+j]*eigenvector[j];}
    double y_norm=0;for(int i=0;i<n;i++)y_norm+=y[i]*y[i];y_norm=sqrt(y_norm);
    if(y_norm<1e-15){free(y);return 0;}
    double lambda_new=0;for(int i=0;i<n;i++)eigenvector[i]=y[i]/y_norm;
    for(int i=0;i<n;i++){double Ax=0;for(int j=0;j<n;j++)Ax+=A[i*n+j]*eigenvector[j];lambda_new+=eigenvector[i]*Ax;}
    if(fabs(lambda_new-lambda)<tol){free(y);return lambda_new;}
    lambda=lambda_new;free(y);}
  return lambda;}

/* Singular Value Decomposition (2x2 analytical) L8 */
int svd_2x2(double a11,double a12,double a21,double a22,double *U,double *S,double *Vt){
  double AtA11=a11*a11+a21*a21;double AtA12=a11*a12+a21*a22;
  double AtA22=a12*a12+a22*a22;
  double trace=AtA11+AtA22;double det=AtA11*AtA22-AtA12*AtA12;
  double disc=sqrt(trace*trace-4*det);
  double lambda1=(trace+disc)/2;double lambda2=(trace-disc)/2;
  if(lambda1<0)lambda1=0;if(lambda2<0)lambda2=0;
  S[0]=sqrt(lambda1);S[1]=sqrt(fmax(lambda2,0));S[2]=0;S[3]=0;
  if(S[0]>1e-15){double v1=AtA12;double v2=lambda1-AtA11;double n=sqrt(v1*v1+v2*v2);if(n>1e-15){Vt[0]=v1/n;Vt[1]=v2/n;}else{Vt[0]=1;Vt[1]=0;}}
  else{Vt[0]=1;Vt[1]=0;}
  Vt[2]=-Vt[1];Vt[3]=Vt[0];
  if(S[0]>1e-15){U[0]=(a11*Vt[0]+a12*Vt[1])/S[0];U[2]=(a21*Vt[0]+a22*Vt[1])/S[0];}
  else{U[0]=1;U[2]=0;}
  if(S[1]>1e-15){U[1]=(a11*Vt[2]+a12*Vt[3])/S[1];U[3]=(a21*Vt[2]+a22*Vt[3])/S[1];}
  else{U[1]=-U[2];U[3]=U[0];}
  return 0;}

/* Matrix Exponential via Taylor Series L8 */
void matrix_exp_taylor(const double *A,int n,double *expA,int terms){
  if(!A||!expA||n<=0||terms<=0)return;
  double *term=(double*)calloc((size_t)(n*n),sizeof(double));
  double *power=(double*)calloc((size_t)(n*n),sizeof(double));
  for(int i=0;i<n*n;i++)expA[i]=0;for(int i=0;i<n;i++)expA[i*n+i]=1.0;
  for(int i=0;i<n*n;i++)power[i]=A[i];double fact=1.0;
  for(int k=1;k<=terms;k++){fact*=k;for(int i=0;i<n*n;i++){term[i]=power[i]/fact;expA[i]+=term[i];}
    double *temp=(double*)calloc((size_t)(n*n),sizeof(double));matrix_mul(power,A,temp,n,n,n);
    for(int i=0;i<n*n;i++)power[i]=temp[i];free(temp);}
  free(term);free(power);}

/* Companion Matrix for Polynomial L8 */
void companion_matrix(const double *coeff,int degree,double *C){
  if(!coeff||!C||degree<=1)return;
  for(int i=0;i<degree*degree;i++)C[i]=0;
  for(int i=1;i<degree;i++)C[i*degree+(i-1)]=1.0;
  for(int j=0;j<degree;j++)C[(degree-1)*degree+j]=-coeff[j]/coeff[degree];}

/* Diagonalization check for 2x2 */
int matrix_diagonalize_2x2(const double *A,double *P,double *D){
  if(!A||!P||!D)return 0;
  double lr1,li1,lr2,li2;int real_ev=matrix_eigen_2x2(A[0],A[1],A[2],A[3],&lr1,&li1,&lr2,&li2);
  if(!real_ev)return 0;
  D[0]=lr1;D[1]=0;D[2]=0;D[3]=lr2;
  P[0]=1;P[2]=1;
  if(fabs(A[2])>1e-15){P[1]=(lr1-A[0])/A[2];P[3]=(lr2-A[0])/A[2];}
  else if(fabs(A[1])>1e-15){P[1]=(lr1-A[3])/A[1];P[3]=(lr2-A[3])/A[1];}
  else{return 0;}return 1;}
