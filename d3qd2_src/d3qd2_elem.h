/*
 * d3qd2_elem.h
 *
 *  Created on: Dec 14, 2018
 *      Author: ohta
 */

#ifndef D3QD2_ELEM_H_
#define D3QD2_ELEM_H_

#include "d3qd2_const.h"

typedef struct sub_domain_data{
  int Ne;   // number of elements
  int *sid; // signed element id

  double complex **P,**dP;    // velocity potential boundary value, P[element id][element node id(0~3)]
  double complex ***pv,***dpv;// particle velocity boundary value, pv[element id][element node id][dimension id(0:pv_x,1:pv_y,2:pv_z)]
}SUBD;

typedef struct boundary_data{
  // periodicity data
  int ps;  // 1 : quasi-periodic boundary condition, 0 : non periodic
  QPD2 qd; // data for quasi-periodic green's function
  // gmsh data
  int Nn;      // number of nodes
  int Ne;      // number of elements 
  double **rn; // node data, rn[node_id][dimension id(0:x,1:y,2:z)]=node coordinate of each dimension
  int **ed;    // element data, ed[element_id][position_id(0,1,2,3)]=node id
  // element data
  int NN;      // number of nodes on element
  int **eni;   // eni[element id][position id(0,1,2,3)]=element node id
  int *md;     // md[element id] = main domain id
  int *sd;     // sd[element id] = sub domain id
  int *gd;     // gd[element id] = group id (gmsh elementary geometrical entity)
  // element constant
  double ***cr,***cw;   // geometric constant c[element id][dimension id(0:x,1:y,2:z)][coefficient id(0:a,1:b,2:c,3:d)]
  double ***ren,***wen; // element node r and w data, ren[element id][element node id(0~3)][dimension id]
  double complex **Pi,**dPi; // incident field data, Ui[element id][element node id][dimension id]
  double zt_44[4],et_44[4],wt_44[4],zt_49[9],et_49[9],wt_49[9]; // gaussian quadrature node and weight for bi-linear element
  double zt_34[4],et_34[4],wt_34[4],zt_37[7],et_37[7],wt_37[7]; // gaussian quadrature node and weight for linear triangular elemen
  double xli[GLN],wli[GLN]; // gauss-legendre node and weight
  double xhi[GHN],whi[GHN]; //
  // sub domain data
  SUBD *sb;
}BOUD;

typedef struct domain_data {
  char ped_fn[128],med_fn[128], msh_fn[128]; // periodicity_datafile_name, medium_datafile_name , mesh data file name
  double rv[3],th;                           // rv:vector defining rotation axis, theta:rotation angle ( for Rodrigues' rotation formula )
  double tv[3];                              // tv:translation vector

  int MN;             // number of mediums
  double *rho0,*c0;   // density of the medium, wave speed in the medium.
  
  double *k0,*K0;       // wave number k0, buld modulus K0 (K0=rho0*c0*c0)
  double complex *k1,*k2;// wave constant (k1=I*omega/K0, k2=I*omega*rho0)

  Apw pw;  // incident field data, aw_pw.h
  BOUD bd; // boundary data
} AQD2;

typedef struct tmp_data{
  double complex k; // wave number
  double zeta,eta;
  double zeta_t,eta_t;
  double r0,r1,th0,th1;
  double sth,cth;
  double rt[3],vt[3],vu[3],aW;
  double cr[3][4],cw[3][3];
  double *xi,*wi,*xhi,*whi; // gauss node and weight data
  int nid,gn,ghn;
  QPD2 *qpd;
}TDATA;

typedef struct coef_matrix{
  int type; // coefficient integral type
  int MN; // material number

  // coefficient matrix file name of each domain
  char **tgfn; // G matrix name of each domain
  char **thfn; // H matrix name of each domain

  // derivative coefficient data file name
  char **tdgfn;
  char **tdhfn;
  char **tdffn;

  size_t nn; // total element node number
  // compressed sparse row (CSR) storage file data
  size_t na; // matrix size na*na
  size_t nnz; // number of none zero element
  char *aval; // value
  char *aptr; // pointer
  char *aidx; // index
  // right hand side (rhs) vector
  char *b; //
}CMD; // for boundary value solver


// -- common function -----------------------------------------------------------------------------------------------------------
int check_element_type(int seid,BOUD *bd); //seid:signed element id, return ELT3 or ELT4 defined number
void r_bd(double *r,int seid,double zeta,double eta,BOUD *bd); // calc coordinate on element
void r_bd_node(double *r,int seid,int nid,BOUD *bd); // calc node coordinate on element, nid: element node id 0~3
void n_bd_node(double *n,int seid,int nid,BOUD *bd); // calc normal vector at element node, nid: element node id 0~3
void r_tz_bd(double *r,double *vtz,int seid,double zeta,double eta,BOUD *bd); // calc coordinate and zeta-direction unit tangential vector on element
void r_te_bd(double *r,double *vte,int seid,double zeta,double eta,BOUD *bd); // calc coordinate and  eta-direction unit tangential vector on element
void Tz_bd(double *vTz,int seid,double zeta,double eta,BOUD *bd); // Tz is zeta direction tangential vector, non unit vector
void Te_bd(double *vTe,int seid,double zeta,double eta,BOUD *bd); // Te is eta  direction tangential vector, non unit vector
void tz_te_bd_node(double *vtz,double *vte,int seid,int nid,BOUD *bd);


// coef_rt(), CC[0~3]:G, CC[4~7]:H, CC[8]:F, type=0:4p, type=1:9p(bi-linear) or 7p(linear-tri), type=2:GLN, type=3:GHN, type=4:DE
void coef_rt(double complex *CC,double *rt,int s,double complex k,int type,BOUD *bd);
void coef_bd(double complex *CC,double *rt,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd);
void coef_bd_node(double complex *CC,int t,int nid,int s,double complex k,int type,BOUD *bd);
// dcoef_rt(), dCC[0~3]:dG/dvt, dCC[4~7]:dH/dvt, dCC[8]:dF/dvt
void dcoef_rt(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,int type,BOUD *bd);
void dcoef_rt_grad(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,int type,BOUD *bd);
// dcoef_bd_tz() and dcoef_bd_te()
void dcoef_bd_tz(double complex *CC,double complex *dCC,double *rt,double *vtz,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd);
void dcoef_bd_te(double complex *CC,double complex *dCC,double *rt,double *vte,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd);
void dcoef_bd_t2(double complex *CC,double complex *dCdtz,double complex *dCdte,
                 double *rt,double *vtz,double *vte,int t,double zeta_t,double eta_t,int s,double complex k,int type, BOUD *bd);
void dcoef_bd_node_t2(double complex *CC,double complex *dCdtz,double complex *dCdte,
                      int t,int nid,double *vtz,double *vte,int s,double complex k,int type,BOUD *bd);

// quasi-peridic function
void q0_coef_rt(double complex *CC,double *rt,int s,double complex k,int type,BOUD *bd);
void q0_dcoef_rt(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,int type,BOUD *bd);
void q0_dcoef_rt_grad(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,int type,BOUD *bd);
void q0_coef_bd(double complex *CC,double *rt,int t,double zeta_t,double eta_t,int s,double complex k,int type,BOUD *bd);
void q0_coef_bd_node(double complex *CC,int t,int nid,int s,double complex k,int type,BOUD *bd);
void q0_dcoef_bd_node_t2(double complex *CC,double complex *dCdtz,double complex *dCdte,
                  int t,int nid,double *vtz,double *vte,int s,double complex k,int type,BOUD *bd);
void q0_dcoef_bd_node_t2_cc(double complex *CC,double complex *dCdtz,double complex *dCdte,
                  int t,int nid,double *vtz,double *vte,int s,double complex k,int type,BOUD *bd); // calc CC[0~8],dCdtz[8],dCdte[8] only

// ---- bi-linear element (4-node quadrangle) -----------------------------------------------------------------------------------
double bil_Mn(int n,double zeta,double eta); // shape function
double bil_Nn(int n,double zeta,double eta); // interpolate function
void bil_copy_elem_const_rw(double cr[][4],double cw[][3],int s,BOUD *bd); // copy element constant of element id s
void bil_convert_CC(double complex *CC); // convert coefficient data (normal vector direction)
void bil_convert_dCC(double complex *dCC); // convert derivative coefficient data (normal vector direction)
void bil_r_zeta_eta(double *r,double zeta,double eta,double cr[][4]); // calc coordinate at zeta,eta on element
void bil_w_zeta_eta(double *w,double zeta,double eta,double cw[][3]); // calc normal vector at zeta,eta on element
void bil_rw_zeta_eta(double *r,double *w,double zeta,double eta,double cr[][4],double cw[][3]);
void bil_t_zeta(double *t,double zeta,double eta,double cr[][4]); // calc zeta-direction tangential vector at zeta,eta on element
void bil_t_eta (double *t,double zeta,double eta,double cr[][4]); // calc  eta-direction tangential vector at zeta,eta on element
int bil_check_plane(double cr[][4],double cw[][3]); // ret 0: element is not plane, ret 1 : element is plane
int bil_check_on_plane(double *rt,double cr[][4],double cw[][3]); // ret 0: rt is not on element plane, ret 1: rt is on element plane
void bil_calc_node_r_th(double *r,double *th,double zeta_t,double eta_t);
// -- bil_coef.c --
// 4p:4 point rule, 9p:9 point rule, GL:GLN*GLN (defined b_const.h) point rule, DE:double exponential method
// CC[0~3]:G, CC[4~7]:H, CC[8]:F
void bil_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // 4-point rule
void bil_coef_9p(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // 9-point rule
void bil_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // GLN * GLN point rule ( defined b_const.h )
void bil_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // GHN * GHN point rule ( defined b_const.h )
void bil_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // Double exponential integration. for test
void bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // rt is on boundary
// -- bil_dcoef.c --
// dCC[0~3]:dG/dv, dCC[4~7]:dH/dv, dCC[8]:dF/dv, tz=T_zeta/|T_zeta|, te=T_eta/|T_eta|
void bil_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void bil_dcoef_9p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void bil_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void bil_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void bil_dcoef_DE0(double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd); // for test
void bil_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd); // for test
void bil_dcoef_bd_tz0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_zeta directional derivative
void bil_dcoef_bd_te0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_eta  directional derivative
void bil_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value, t_zeta directional derivative
void bil_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value, t_eta  directional derivative
void bil_dcoef_bd_tz_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_zeta directional derivative
void bil_dcoef_bd_te_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_eta directional derivative
void bil_dcoef_bd_tz_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_zeta direction derivative, for test
void bil_dcoef_bd_te_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_eta direction derivative, for test
// -- bil_dcoef_v2.c --
void bil_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void bil_dcoef_v2_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void bil_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void bil_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void bil_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void bil_dcoef_bd_vt2_PV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
void bil_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// --bil_dcoef_grad.c --
void bil_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void bil_dcoef_grad_9p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void bil_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void bil_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);


// ----  linear triangular element (3-node triangle) ----------------------------------------------------------------------------
double lit_Mn(int n,double zeta,double eta); // shape function
double lit_Nn(int n,double zeta,double eta); // interpolate function
void lit_copy_elem_const_rw(double cr[][4],double cw[][3],int s,BOUD *bd); // copy element constant of element id s
void lit_convert_CC(double complex *CC); // convert coefficient data (normal vector direction)
void lit_convert_dCC(double complex *dCC);// convert derivative coefficient data (normal vector direction)
void lit_r_zeta_eta(double *r,double zeta,double eta,double cr[][4]); // calc coordinate at zeta,eta on element
void lit_w_zeta_eta(double *w,double zeta,double eta,double cw[][3]); // calc normal vector at zeta,eta on element
void lit_rw_zeta_eta(double *r,double *w,double zeta,double eta,double cr[][4],double cw[][3]);
void lit_t_zeta(double *t,double zeta,double eta,double cr[][4]); // calc zeta-direction tangential vector at zeta,eta on element
void lit_t_eta (double *t,double zeta,double eta,double cr[][4]); // calc  eta-direction tangential vector at zeta,eta on element
int lit_check_on_plane(double *rt,double cr[][4],double cw[][3]); // ret 0: rt is not on element plane, rt 1: rt is on element plane
void lit_calc_node_r_th(double *r,double *th,double zeta_t,double eta_t);
// -- lit_coef.c --
// 4p:4 point rule, 7p:7 point rule, GL:3 part GLN*GLN (defined b_const.h) point rule, DE:double exponential method
void lit_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void lit_coef_7p(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void lit_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void lit_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void lit_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // for test
void lit_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // rt is on boundary
// -- lit_dcoef.c --
void lit_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void lit_dcoef_7p(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void lit_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void lit_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd);
void lit_dcoef_DE0(double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd); // for test
void lit_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *vt,int s,double complex k,BOUD *bd); // for test
void lit_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_zeta direction derivative, for test
void lit_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_eta direction derivative, for test
void lit_dcoef_bd_tz0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_zeta directional derivative
void lit_dcoef_bd_te0_PV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_eta directional derivative
void lit_dcoef_bd_tz_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_zeta direction derivative, for test
void lit_dcoef_bd_te_DV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // differential value of t_eta direction derivative, for test
void lit_dcoef_bd_tz_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_zeta directional derivative
void lit_dcoef_bd_te_PV(double complex *CC,double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // principal value, t_eta directional derivative
// -- lit_dcoef_v2.c --
void lit_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void lit_dcoef_v2_7p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void lit_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void lit_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void lit_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *vt0,double *vt1,int s,double complex k,BOUD *bd);
void lit_dcoef_bd_vt2_PV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
void lit_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// --lit_dcoef_grad.c --
void lit_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void lit_dcoef_grad_7p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void lit_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void lit_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);


// -- quasi-periodic --------------------------------------------------------
// -- bi-linear element
// q0_bil_coef.c
void q0_bil_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // 4-point rule
void q0_bil_coef_9p(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // 9-point rule
void q0_bil_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // GLN * GLN point rule ( defined b_const.h )
void q0_bil_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // GHN * GHN point rule ( defined b_const.h )
void q0_bil_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // Double exponential integration. for test
void q0_bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // rt is on boundary
// q0_bil_dcoef.c
void q0_bil_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_9p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd); // for test
// q0_bil_dcoef_v2.c
void q0_bil_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// q0_bil_dcoef_v2_cc.c
void q0_bil_dcoef_v2_cc_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_cc_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_cc_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_cc_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_v2_cc_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_bd_vt2_cc_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// q0_bil_dcoef_grad.c
void q0_bil_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_grad_9p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_bil_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);


// -- linear-triangular element
// q0_lit_coef.c
void q0_lit_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void q0_lit_coef_7p(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void q0_lit_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void q0_lit_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd);
void q0_lit_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd); // for test
void q0_lit_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// q0_lit_dcoef.c
void q0_lit_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_7p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd); // for test
// q0_lit_dcoef_v2.c
void q0_lit_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_7p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// q0_lit_dcoef_v2_cc.c
void q0_lit_dcoef_v2_cc_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_cc_7p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_cc_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_cc_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_v2_cc_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_bd_vt2_cc_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
// q0_lit_dcoef_grad.c
void q0_lit_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_grad_7p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);
void q0_lit_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd);


#endif
