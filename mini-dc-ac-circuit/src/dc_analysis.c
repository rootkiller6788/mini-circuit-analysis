#include "dc_analysis.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "circuit_utils.h"
#include "ac_analysis.h"
#include "circuit_matrix.h"

/* Ohm Law 1827 L4 */
double ohm_law_voltage(double I,double R){return I*R;}
double ohm_law_current(double V,double R){if(R==0.0){errno=ERANGE;return (V>=0.0)?INFINITY:-INFINITY;}return V/R;}
double ohm_law_resistance(double V,double I){if(I==0.0)return INFINITY;return V/I;}

/* Joule Law 1841 L4 */
double joule_power_from_vi(double V,double I){return V*I;}
double joule_power_from_ir(double I,double R){return I*I*R;}
double joule_power_from_vr(double V,double R){if(R==0.0)return INFINITY;return (V*V)/R;}
double conductance_from_resistance(double R){if(R==0.0)return INFINITY;return 1.0/R;}
double resistance_from_conductance(double G){if(G==0.0)return INFINITY;return 1.0/G;}

/* KCL KVL 1845 L4 */
double kcl_residual(const double *c,int n){if(!c||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++)s+=c[i];return s;}
double kvl_residual(const double *v,int n){if(!v||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++)s+=v[i];return s;}

/* Voltage Divider L2 */
double voltage_divider(double Vin,double R1,double R2){double Rt=R1+R2;if(Rt==0.0)return 0.0;return Vin*R2/Rt;}
int voltage_divider_n(double Vin,const double *R,double *V,int n){
  if(!R||!V||n<=0)return -1;double Rt=0;
  for(int i=0;i<n;i++){if(R[i]<0)return -1;Rt+=R[i];}
  if(Rt==0.0){for(int i=0;i<n;i++)V[i]=0;return 0;}
  for(int i=0;i<n;i++)V[i]=Vin*R[i]/Rt;return 0;}

/* Current Divider L2 */
double current_divider_2r(double It,double Rother,double Rthis){double Rs=Rthis+Rother;if(Rs==0.0)return INFINITY;return It*Rother/Rs;}
int current_divider_n(double It,const double *R,double *Io,int n){
  if(!R||!Io||n<=0)return -1;
  double *G=(double*)malloc((size_t)n*sizeof(double));if(!G)return -1;
  double Gt=0;for(int i=0;i<n;i++){if(R[i]<0){free(G);return -1;}G[i]=(R[i]==0.0)?INFINITY:(1.0/R[i]);Gt+=G[i];}
  if(Gt==0.0||isinf(Gt)){for(int i=0;i<n;i++)Io[i]=0;}
  else{for(int i=0;i<n;i++)Io[i]=It*G[i]/Gt;}free(G);return 0;}

/* Equivalent R L C L2 */
double resistance_series(const double *R,int n){if(!R||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(R[i]<0)return -1.0;s+=R[i];}return s;}
double resistance_parallel_2(double R1,double R2){double s=R1+R2;if(s==0.0)return 0.0;return(R1*R2)/s;}
double resistance_parallel(const double *R,int n){if(!R||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(R[i]<0)return -1.0;if(R[i]==0.0)return 0.0;s+=1.0/R[i];}if(s==0.0)return INFINITY;return 1.0/s;}
double capacitance_series(const double *C,int n){if(!C||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(C[i]<0)return -1.0;if(C[i]==0.0)return 0.0;s+=1.0/C[i];}if(s==0.0)return INFINITY;return 1.0/s;}
double capacitance_parallel(const double *C,int n){if(!C||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(C[i]<0)return -1.0;s+=C[i];}return s;}
double inductance_series(const double *L,int n){if(!L||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(L[i]<0)return -1.0;s+=L[i];}return s;}
double inductance_parallel(const double *L,int n){if(!L||n<=0)return 0.0;double s=0;for(int i=0;i<n;i++){if(L[i]<0)return -1.0;if(L[i]==0.0)return 0.0;s+=1.0/L[i];}if(s==0.0)return INFINITY;return 1.0/s;}

/* Gaussian Elimination L5 O(n^3) */
static int gaussian_elimination(double *A,const double *b,double *x,int n){
  if(!A||!b||!x||n<=0)return -1;
  for(int i=0;i<n;i++)x[i]=b[i];
  for(int col=0;col<n;col++){
    int max_row=col;double max_val=fabs(A[col*n+col]);
    for(int row=col+1;row<n;row++){double val=fabs(A[row*n+col]);if(val>max_val){max_val=val;max_row=row;}}
    if(max_val<1e-15)return -1;
    if(max_row!=col){
      for(int j=col;j<n;j++){double tmp=A[col*n+j];A[col*n+j]=A[max_row*n+j];A[max_row*n+j]=tmp;}
      double tmp=x[col];x[col]=x[max_row];x[max_row]=tmp;}
    double pivot=A[col*n+col];
    for(int row=col+1;row<n;row++){double factor=A[row*n+col]/pivot;if(factor==0.0)continue;A[row*n+col]=0.0;for(int j=col+1;j<n;j++)A[row*n+j]-=factor*A[col*n+j];x[row]-=factor*x[col];}}
  for(int i=n-1;i>=0;i--){double sum=x[i];for(int j=i+1;j<n;j++)sum-=A[i*n+j]*x[j];x[i]=sum/A[i*n+i];}
  return 0;}

/* Nodal Analysis L5 */
int nodal_analysis_solve(int n,const double *G,const double *I,double *V){
  if(!G||!I||!V||n<=0)return -1;
  int size=n*n;double *A=(double*)malloc((size_t)size*sizeof(double));if(!A)return -1;
  for(int i=0;i<size;i++)A[i]=G[i];
  int result=gaussian_elimination(A,I,V,n);free(A);return result;}

/* Build Conductance Matrix L5 MNA */
int nodal_build_conductance(const Circuit_t *ckt,double *G,double *I){
  if(!ckt||!G||!I)return -1;if(ckt->n_nodes<2)return -1;
  int n=ckt->n_nodes-1;int n_size=n*n;
  for(int i=0;i<n_size;i++)G[i]=0.0;for(int i=0;i<n;i++)I[i]=0.0;
  int *node_to_mat=(int*)malloc((size_t)ckt->n_nodes*sizeof(int));
  int *mat_to_node=(int*)malloc((size_t)n*sizeof(int));
  if(!node_to_mat||!mat_to_node){free(node_to_mat);free(mat_to_node);return -1;}
  int mat_idx=0;
  for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)node_to_mat[i]=-1;else{node_to_mat[i]=mat_idx;mat_to_node[mat_idx]=i;mat_idx++;}}
  for(int k=0;k<ckt->n_resistors;k++){
    const Resistor_t *r=&ckt->resistors[k];int ni=r->t[0].node_id;int nj=r->t[1].node_id;
    if(r->R<=0.0)continue;double g=1.0/r->R;
    if(node_to_mat[ni]>=0)G[node_to_mat[ni]*n+node_to_mat[ni]]+=g;
    if(node_to_mat[nj]>=0)G[node_to_mat[nj]*n+node_to_mat[nj]]+=g;
    if(node_to_mat[ni]>=0&&node_to_mat[nj]>=0){G[node_to_mat[ni]*n+node_to_mat[nj]]-=g;G[node_to_mat[nj]*n+node_to_mat[ni]]-=g;}}
  for(int k=0;k<ckt->n_dc_isrc;k++){
    const DCCurrentSource_t *src=&ckt->dc_isrc[k];int ni=src->t[0].node_id;int nj=src->t[1].node_id;
    if(node_to_mat[ni]>=0)I[node_to_mat[ni]]+=src->I_dc;
    if(node_to_mat[nj]>=0)I[node_to_mat[nj]]-=src->I_dc;}
  for(int k=0;k<ckt->n_dc_vsrc;k++){
    const DCVoltageSource_t *src=&ckt->dc_vsrc[k];int np=src->t[0].node_id;int nm=src->t[1].node_id;
    if(nm==ckt->ground_node&&node_to_mat[np]>=0){int idx=node_to_mat[np];for(int j=0;j<n;j++)G[idx*n+j]=0.0;G[idx*n+idx]=1.0;I[idx]=src->V_dc;}
    else if(np==ckt->ground_node&&node_to_mat[nm]>=0){int idx=node_to_mat[nm];for(int j=0;j<n;j++)G[idx*n+j]=0.0;G[idx*n+idx]=1.0;I[idx]=-src->V_dc;}}
  free(node_to_mat);free(mat_to_node);return n;}

/* Mesh Analysis L5 */
int mesh_analysis_solve(int m,const double *R,const double *V,double *I_mesh){
  if(!R||!V||!I_mesh||m<=0)return -1;int size=m*m;
  double *A=(double*)malloc((size_t)size*sizeof(double));if(!A)return -1;
  for(int i=0;i<size;i++)A[i]=R[i];
  int result=gaussian_elimination(A,V,I_mesh,m);free(A);return result;}
int mesh_build_resistance(const Circuit_t *ckt,double *R_mat,double *V_vec){
  if(!ckt||!R_mat||!V_vec)return -1;if(ckt->n_loops<=0)return 0;int m=ckt->n_loops;int ms=m*m;
  for(int i=0;i<ms;i++)R_mat[i]=0.0;for(int i=0;i<m;i++)V_vec[i]=0.0;
  for(int i=0;i<m;i++){const Loop_t *li=&ckt->loops[i];
    for(int bi=0;bi<li->branch_count;bi++){int br_idx=li->branches[bi];const Branch_t *br=&ckt->branches[br_idx];
      if(br->elem_type==ELEM_RESISTOR)R_mat[i*m+i]+=br->value;
      for(int j=i+1;j<m;j++){const Loop_t *lj=&ckt->loops[j];
        for(int bj=0;bj<lj->branch_count;bj++){if(lj->branches[bj]==br_idx){
          double sign=(li->directions[bi]==lj->directions[bj])?-1.0:1.0;
          if(br->elem_type==ELEM_RESISTOR){R_mat[i*m+j]+=sign*br->value;R_mat[j*m+i]+=sign*br->value;}}}}
      if(br->elem_type==ELEM_DC_VOLTAGE_SRC)V_vec[i]+=li->directions[bi]*br->value;}}
  return m;}

/* Thevenin Norton L4 L6 */
TheveninEquiv_t thevenin_from_oc_sc(double Voc,double Isc){TheveninEquiv_t th;th.V_th=Voc;th.R_th=(Isc==0.0)?INFINITY:(Voc/Isc);return th;}
TheveninEquiv_t thevenin_from_voc_rth(double Voc,double Rth){TheveninEquiv_t th;th.V_th=Voc;th.R_th=Rth;return th;}
NortonEquiv_t norton_from_sc_oc(double Isc,double Voc){NortonEquiv_t no;no.I_n=Isc;no.R_n=(Isc==0.0)?INFINITY:(Voc/Isc);return no;}
NortonEquiv_t thevenin_to_norton(TheveninEquiv_t th){NortonEquiv_t no;no.I_n=(th.R_th==0.0)?INFINITY:(th.V_th/th.R_th);no.R_n=th.R_th;return no;}
TheveninEquiv_t norton_to_thevenin(NortonEquiv_t no){TheveninEquiv_t th;th.V_th=no.I_n*no.R_n;th.R_th=no.R_n;return th;}

/* Max Power Transfer L4 Jacobi 1840 */
double max_power_transfer_load_r(double Rth){return Rth;}
double max_power_transfer_pmax(double Vth,double Rth){if(Rth==0.0)return INFINITY;return(Vth*Vth)/(4.0*Rth);}
double load_power(double Vth,double Rth,double Rl){double Rt=Rth+Rl;if(Rt==0.0)return 0.0;return(Vth*Vth*Rl)/(Rt*Rt);}
double load_efficiency(double Rth,double Rl){double Rt=Rth+Rl;if(Rt==0.0)return 0.0;return Rl/Rt;}
double max_power_efficiency(void){return 0.5;}

/* Superposition L4 Helmholtz 1853 */
int superposition_solve(const Circuit_t *ckt,double *V_nodes,int n_nodes){
  if(!ckt||!V_nodes||n_nodes<=0)return -1;
  for(int i=0;i<n_nodes;i++)V_nodes[i]=0.0;
  int n_total_src=ckt->n_dc_vsrc+ckt->n_dc_isrc;if(n_total_src==0)return 0;
  for(int s=0;s<ckt->n_dc_vsrc;s++){
    int n=ckt->n_nodes-1;
    double *G=(double*)calloc((size_t)(n*n),sizeof(double));
    double *I=(double*)calloc((size_t)n,sizeof(double));
    double *Vp=(double*)calloc((size_t)n,sizeof(double));
    if(!G||!I||!Vp){free(G);free(I);free(Vp);return -1;}
    int *nm=(int*)malloc((size_t)ckt->n_nodes*sizeof(int));
    int mi=0;for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)nm[i]=-1;else nm[i]=mi++;}
    for(int k=0;k<ckt->n_resistors;k++){const Resistor_t *r=&ckt->resistors[k];
      int ni=r->t[0].node_id;int nj=r->t[1].node_id;if(r->R<=0.0)continue;double g=1.0/r->R;
      if(nm[ni]>=0)G[nm[ni]*n+nm[ni]]+=g;if(nm[nj]>=0)G[nm[nj]*n+nm[nj]]+=g;
      if(nm[ni]>=0&&nm[nj]>=0){G[nm[ni]*n+nm[nj]]-=g;G[nm[nj]*n+nm[ni]]-=g;}}
    const DCVoltageSource_t *src=&ckt->dc_vsrc[s];int np=src->t[0].node_id;int nmn=src->t[1].node_id;
    if(nmn==ckt->ground_node&&nm[np]>=0){int idx=nm[np];for(int j=0;j<n;j++)G[idx*n+j]=0.0;G[idx*n+idx]=1.0;I[idx]=src->V_dc;}
    else if(np==ckt->ground_node&&nm[nmn]>=0){int idx=nm[nmn];for(int j=0;j<n;j++)G[idx*n+j]=0.0;G[idx*n+idx]=1.0;I[idx]=-src->V_dc;}
    if(gaussian_elimination(G,I,Vp,n)==0){for(int i=0;i<ckt->n_nodes;i++){if(nm[i]>=0)V_nodes[i]+=Vp[nm[i]];}}
    free(nm);free(G);free(I);free(Vp);}
  for(int s=0;s<ckt->n_dc_isrc;s++){
    int n=ckt->n_nodes-1;
    double *G=(double*)calloc((size_t)(n*n),sizeof(double));
    double *I=(double*)calloc((size_t)n,sizeof(double));
    double *Vp=(double*)calloc((size_t)n,sizeof(double));
    if(!G||!I||!Vp){free(G);free(I);free(Vp);return -1;}
    int *nm=(int*)malloc((size_t)ckt->n_nodes*sizeof(int));
    int mi=0;for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)nm[i]=-1;else nm[i]=mi++;}
    for(int k=0;k<ckt->n_resistors;k++){const Resistor_t *r=&ckt->resistors[k];
      int ni=r->t[0].node_id;int nj=r->t[1].node_id;if(r->R<=0.0)continue;double g=1.0/r->R;
      if(nm[ni]>=0)G[nm[ni]*n+nm[ni]]+=g;if(nm[nj]>=0)G[nm[nj]*n+nm[nj]]+=g;
      if(nm[ni]>=0&&nm[nj]>=0){G[nm[ni]*n+nm[nj]]-=g;G[nm[nj]*n+nm[ni]]-=g;}}
    const DCCurrentSource_t *src=&ckt->dc_isrc[s];int ni=src->t[0].node_id;int nj=src->t[1].node_id;
    if(nm[ni]>=0)I[nm[ni]]+=src->I_dc;if(nm[nj]>=0)I[nm[nj]]-=src->I_dc;
    if(gaussian_elimination(G,I,Vp,n)==0){for(int i=0;i<ckt->n_nodes;i++){if(nm[i]>=0)V_nodes[i]+=Vp[nm[i]];}}
    free(nm);free(G);free(I);free(Vp);}
  if(ckt->ground_node>=0&&ckt->ground_node<n_nodes)V_nodes[ckt->ground_node]=0.0;
  return 0;}

/* Delta-Wye Transform L5 Kennelly 1899 */
void delta_to_wye(double Ra,double Rb,double Rc,double *R1,double *R2,double *R3){
  double sum=Ra+Rb+Rc;if(sum==0.0){if(R1)*R1=0;if(R2)*R2=0;if(R3)*R3=0;return;}
  if(R1)*R1=(Rb*Rc)/sum;if(R2)*R2=(Ra*Rc)/sum;if(R3)*R3=(Ra*Rb)/sum;}
void wye_to_delta(double R1,double R2,double R3,double *Ra,double *Rb,double *Rc){
  double num=R1*R2+R2*R3+R3*R1;
  if(Ra)*Ra=(R1==0.0)?INFINITY:(num/R1);if(Rb)*Rb=(R2==0.0)?INFINITY:(num/R2);if(Rc)*Rc=(R3==0.0)?INFINITY:(num/R3);}

/* Wheatstone Bridge L6 Wheatstone 1843 */
double wheatstone_bridge_vout(double Vin,double R1,double R2,double R3,double Rx){
  double Rls=R1+R2,Rrs=R3+Rx;double Va=(Rls==0.0)?0.0:Vin*R2/Rls;double Vb=(Rrs==0.0)?0.0:Vin*Rx/Rrs;return Va-Vb;}
double wheatstone_bridge_rx_balanced(double R1,double R2,double R3){if(R1==0.0)return INFINITY;return R3*R2/R1;}
int wheatstone_is_balanced(double R1,double R2,double R3,double Rx,double tol){
  double left=R1*Rx,right=R2*R3;double diff=fabs(left-right);double avg=(fabs(left)+fabs(right))/2.0;
  if(avg==0.0)return 1;return(diff/avg)<tol;}

/* DC Operating Point L5 MNA */
int dc_operating_point(const Circuit_t *ckt,double *node_voltages,double *branch_currents){
  if(!ckt||!node_voltages)return -1;int n=ckt->n_nodes-1;if(n<=0)return -1;
  double *G=(double*)calloc((size_t)(n*n),sizeof(double));
  double *Iv=(double*)calloc((size_t)n,sizeof(double));
  double *V=(double*)calloc((size_t)n,sizeof(double));
  if(!G||!Iv||!V){free(G);free(Iv);free(V);return -1;}
  int neq=nodal_build_conductance(ckt,G,Iv);int ret=nodal_analysis_solve(neq,G,Iv,V);
  if(ret==0){for(int i=0;i<ckt->n_nodes;i++){if(i==ckt->ground_node)node_voltages[i]=0.0;else{int mat_i=0;for(int j=0;j<i;j++){if(j!=ckt->ground_node)mat_i++;}if(mat_i<n)node_voltages[i]=V[mat_i];}}}
  if(branch_currents&&ret==0){for(int k=0;k<ckt->n_branches;k++){const Branch_t *b=&ckt->branches[k];double vd=node_voltages[b->node_from]-node_voltages[b->node_to];if(b->elem_type==ELEM_RESISTOR&&b->value>0.0)branch_currents[k]=vd/b->value;else branch_currents[k]=0.0;}}
  free(G);free(Iv);free(V);return ret;}

/* DC Analysis Additional Methods L5 */
int dc_sweep_solve(const Circuit_t *ckt,double *V_sweep,int n_steps,double V_start,double V_stop,double *V_out,int out_node){
  if(!ckt||!V_sweep||!V_out||n_steps<=0)return -1;
  for(int step=0;step<n_steps;step++){double V=V_start+(V_stop-V_start)*step/(n_steps-1);V_sweep[step]=V;
    Circuit_t ckt_copy=*ckt;
    int n=ckt->n_nodes-1;if(n>0){
      double *G=(double*)calloc((size_t)(n*n),sizeof(double));
      double *I=(double*)calloc((size_t)n,sizeof(double));
      double *Vv=(double*)calloc((size_t)n,sizeof(double));
      nodal_build_conductance(ckt,G,I);nodal_analysis_solve(n,G,I,Vv);
      if(out_node>=0&&out_node<ckt->n_nodes){if(out_node==ckt->ground_node)V_out[step]=0;else{int mi=0;for(int j=0;j<out_node;j++)if(j!=ckt->ground_node)mi++;V_out[step]=(mi<n)?Vv[mi]:0;}}
      free(G);free(I);free(Vv);}}return 0;}

/* Sensitivity Analysis L8 */
double dc_sensitivity_resistor(const Circuit_t *ckt,int resistor_idx,int meas_node){
  if(!ckt||resistor_idx<0||resistor_idx>=ckt->n_resistors)return 0;
  double R_orig=ckt->resistors[resistor_idx].R;double delta=0.01*R_orig;if(delta<1e-6)delta=1e-6;
  double *V_nom=(double*)calloc((size_t)ckt->n_nodes,sizeof(double));
  dc_operating_point(ckt,V_nom,NULL);double V_nom_meas=V_nom[meas_node];
  Circuit_t ckt_pert=*ckt;ckt_pert.resistors[resistor_idx].R=R_orig+delta;
  double *V_pert=(double*)calloc((size_t)ckt_pert.n_nodes,sizeof(double));
  dc_operating_point(&ckt_pert,V_pert,NULL);double V_pert_meas=V_pert[meas_node];
  free(V_nom);free(V_pert);return(V_pert_meas-V_nom_meas)/delta;}

/* Monte Carlo Tolerance Analysis L8 */
int monte_carlo_dc(const Circuit_t *ckt,int n_runs,double tolerances[],double *V_mean,double *V_stddev,int n_nodes){
  if(!ckt||n_runs<=0)return -1;
  for(int i=0;i<n_nodes;i++){V_mean[i]=0;V_stddev[i]=0;}
  double **V_samples=(double**)malloc((size_t)n_runs*sizeof(double*));
  for(int r=0;r<n_runs;r++){V_samples[r]=(double*)calloc((size_t)n_nodes,sizeof(double));
    Circuit_t ckt_mc=*ckt;
    for(int k=0;k<ckt->n_resistors&&k<n_nodes;k++){
      ckt_mc.resistors[k].R=ckt->resistors[k].R*(1.0+tolerances[k]*(random_uniform(0,1)-0.5)*2.0);}
    dc_operating_point(&ckt_mc,V_samples[r],NULL);
    for(int i=0;i<n_nodes;i++)V_mean[i]+=V_samples[r][i];}
  for(int i=0;i<n_nodes;i++)V_mean[i]/=n_runs;
  for(int r=0;r<n_runs;r++)for(int i=0;i<n_nodes;i++){double d=V_samples[r][i]-V_mean[i];V_stddev[i]+=d*d;}
  for(int i=0;i<n_nodes;i++)V_stddev[i]=sqrt(V_stddev[i]/n_runs);
  for(int r=0;r<n_runs;r++)free(V_samples[r]);free(V_samples);return 0;}

/* Nodal Admittance for AC sweep wrapper L5 */
int ac_sweep_solve(const Circuit_t *ckt,double f_start,double f_stop,int n_points,int out_node,double *magnitude,double *phase_deg){
  if(!ckt||!magnitude||!phase_deg||n_points<=0)return -1;
  for(int i=0;i<n_points;i++){double f=f_start+(f_stop-f_start)*i/(n_points-1);double w=2*M_PI*f;
    int n=ckt->n_nodes-1;if(n>0){
      double complex *Y=(double complex*)calloc((size_t)(n*n),sizeof(double complex));
      double complex *Is=(double complex*)calloc((size_t)n,sizeof(double complex));
      double complex *V=(double complex*)calloc((size_t)n,sizeof(double complex));
      ac_nodal_build_Y(ckt,w,Y,Is);ac_nodal_solve(n,Y,Is,V);
      double complex Vout=0;if(out_node!=ckt->ground_node){int mi=0;for(int j=0;j<out_node;j++)if(j!=ckt->ground_node)mi++;Vout=(mi<n)?V[mi]:0;}
      magnitude[i]=cabs(Vout);phase_deg[i]=carg(Vout)*180.0/M_PI;
      free(Y);free(Is);free(V);}}return 0;}

/* Nonlinear DC Analysis L8 */
double diode_ideal_current(double V_diode,double Is,double V_T){if(V_T<=0)return 0;return Is*(exp(V_diode/V_T)-1.0);}
double diode_ideal_voltage(double I_diode,double Is,double V_T){if(I_diode<=0||Is<=0||V_T<=0)return 0;return V_T*log(I_diode/Is+1.0);}
double diode_small_signal_resistance(double I_dc,double V_T){if(I_dc<=0||V_T<=0)return INFINITY;return V_T/I_dc;}
double bjt_collector_current(double Is,double V_BE,double V_T,double V_CE,double V_A){
  if(V_T<=0)return 0;return Is*exp(V_BE/V_T)*(1+V_CE/V_A);}
double bjt_transconductance(double I_C,double V_T){if(V_T<=0)return 0;return I_C/V_T;}
double bjt_base_emitter_resistance(double beta,double gm){if(gm<=0)return INFINITY;return beta/gm;}
double mosfet_drain_current_saturation(double k,double V_GS,double V_TH,double lambda,double V_DS){
  double Vov=V_GS-V_TH;if(Vov<=0)return 0;return 0.5*k*Vov*Vov*(1+lambda*V_DS);}
double mosfet_transconductance(double k,double V_GS,double V_TH){double Vov=V_GS-V_TH;if(Vov<=0)return 0;return k*Vov;}
double mosfet_output_resistance(double lambda,double I_D){if(lambda<=0||I_D<=0)return INFINITY;return 1.0/(lambda*I_D);}

/* Newton-Raphson DC Solver L8 */
int newton_raphson_dc(int n_eq,int max_iter,double tol,void (*func)(const double*,double*),
                        void (*jacobian)(const double*,double*),double *x){
  if(!func||!jacobian||!x||n_eq<=0)return -1;double *f=(double*)malloc((size_t)n_eq*sizeof(double));
  double *J=(double*)malloc((size_t)(n_eq*n_eq)*sizeof(double));double *dx=(double*)malloc((size_t)n_eq*sizeof(double));
  if(!f||!J||!dx){free(f);free(J);free(dx);return -1;}
  for(int iter=0;iter<max_iter;iter++){func(x,f);jacobian(x,J);
    for(int i=0;i<n_eq;i++)dx[i]=-f[i];
    int ret=matrix_solve(J,f,dx,n_eq);if(ret<0){free(f);free(J);free(dx);return -1;}
    double err=0;for(int i=0;i<n_eq;i++){x[i]+=dx[i];err+=fabs(dx[i]);}
    if(err<tol){free(f);free(J);free(dx);return iter+1;}}
  free(f);free(J);free(dx);return -2;}

/* Circuit Analysis Helper L6 */
int count_loops_planar(const Circuit_t *ckt){if(!ckt)return 0;return ckt->n_branches-ckt->n_nodes+1;}
int is_planar_estimate(const Circuit_t *ckt){if(!ckt)return 0;return ckt->n_branches<=3*ckt->n_nodes-6;}
double power_consumption_total(const Circuit_t *ckt,const double *V,const double *I){
  if(!ckt||!V||!I)return 0;double P=0;for(int i=0;i<ckt->n_branches;i++)P+=V[i]*I[i];return P;}

/* Bus and Interconnect Modeling L7 */
double voltage_drop_wire(double I,double R_per_m,double length_m){return I*R_per_m*length_m;}
double wire_resistance_from_awg_length(int awg,double length_m){return awg_to_resistance_per_meter(awg)*length_m;}
double ground_bounce_voltage(double L_ground_nH,double di_dt_A_per_ns){return L_ground_nH*di_dt_A_per_ns;}
double decoupling_capacitor_min(double I_transient,double dt,double dV_allowed){
  if(dV_allowed<=0||dt<=0)return 0;return I_transient*dt/dV_allowed;}
double simultaneous_switching_noise(double L_eff,double N_drivers,double I_per_driver,double trise){
  if(trise<=0)return INFINITY;return L_eff*N_drivers*I_per_driver/trise;}

/* Temperature Sensing Circuits L7 */
double ntc_thermistor_resistance(double R25,double B_value,double T_celsius){
  double T_K=T_celsius+273.15;return R25*exp(B_value*(1.0/T_K-1.0/298.15));}
double ntc_voltage_divider_output(double Vcc,double R_fixed,double R_ntc){return Vcc*R_ntc/(R_fixed+R_ntc);}
double pt100_resistance(double T_celsius){return 100.0*(1.0+3.9083e-3*T_celsius-5.775e-7*T_celsius*T_celsius);}
double thermocouple_voltage_k_type(double T_celsius){return T_celsius*0.041;}
double thermocouple_cold_junction_compensation(double V_measured,double T_cj){return V_measured+thermocouple_voltage_k_type(T_cj);}

/* Strain Gauge Bridge L7 */
double strain_gauge_resistance_change(double R_nominal,double GF,double strain){return R_nominal*GF*strain;}
double strain_gauge_bridge_output(double V_excitation,double GF,double strain){
  return V_excitation*GF*strain/4.0;}
double strain_gauge_quarter_bridge_completion(double R_gauge,double R_completion){return R_completion;}

/* Photodiode Circuit L7 */
double photodiode_transimpedance_gain(double R_feedback){return R_feedback;}
double photodiode_bandwidth(double R_feedback,double C_feedback){return 1.0/(2*M_PI*R_feedback*C_feedback);}
double photodiode_responsivity_A_W(double I_photo,double P_optical){if(P_optical<=0)return 0;return I_photo/P_optical;}
double photodiode_NEP(double noise_current,double responsivity){if(responsivity<=0)return 0;return noise_current/responsivity;}

/* Hall Effect Sensor L7 */
double hall_voltage(double I_bias,double B_tesla,double d_thickness,double n_carriers,double q_electron){
  if(d_thickness<=0||n_carriers<=0)return 0;return I_bias*B_tesla/(n_carriers*q_electron*d_thickness);}
double hall_sensitivity_V_per_T(double V_hall,double B){if(B<=0)return 0;return V_hall/B;}
double gmr_sensor_resistance_change(double R0,double GMR_ratio,double B,double B_sat){
  if(B_sat<=0)return R0;double ratio=(B/B_sat<1)?(B/B_sat):1;return R0*(1.0-GMR_ratio*ratio);}

/* Current Mirror Circuits L7 */
double current_mirror_output_current(double I_ref,double W2_L2_ratio,double W1_L1_ratio){if(W1_L1_ratio<=0)return 0;return I_ref*W2_L2_ratio/W1_L1_ratio;}
double current_mirror_output_resistance(double V_A,double I_out){if(I_out<=0)return INFINITY;return V_A/I_out;}
double current_mirror_matching_error(double dVth,double V_ov,double gm,double gds){return gm*dVth;}
double widlar_current_source_R(double V_T,double I_out,double I_ref){if(I_out<=0||I_ref<=0)return 0;return V_T/I_out*log(I_ref/I_out);}

/* Differential Pair L7 */
double diff_pair_gain(double gm,double R_load){return gm*R_load;}
double diff_pair_input_range(double Vcc,double V_T,double I_tail,double R_tail){return Vcc-4*V_T-I_tail*R_tail;}
double diff_pair_cmrr_db(double Adm,double Acm){if(Acm<=0)return INFINITY;return 20*log10(Adm/Acm);}
double diff_pair_offset_voltage(double dVth,double dWL_ratio,double V_ov){return dVth+0.5*V_ov*dWL_ratio;}

/* Bandgap Voltage Reference L7 */
double bandgap_ptat_voltage(double k,double T,double q,double n){return k*T/q*log(n);}
double bandgap_vref(double Vbe,double K,double delta_Vbe){return Vbe+K*delta_Vbe;}
double bandgap_curvature_correction(double V_T,double T,double T0){return V_T*(T/T0-1)*log(T/T0);}
