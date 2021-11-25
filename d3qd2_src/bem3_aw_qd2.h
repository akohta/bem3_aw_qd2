/*
 * bem3_aw_qd2.h
 *
 *  Created on: Nov 22, 2021
 *      Author: ohta
 */
#ifndef BEM3_AW_QD2_H_
#define BEM3_AW_QD2_H_
 
#include "d3qd2_elem.h"

// -- bem3_aw_qd2.c --
void read_aqd2(int argc,char **argv,AQD2 *qd); // read datafile.
void print_aqd2(AQD2 *qd);                     // print data.
void initialize_aqd2(AQD2 *qd);                // memory allocation and initalize coefficients.
void finalize_aqd2(AQD2 *qd);                  // memory free.
int domain_id_aqd2(double *rt,AQD2 *qd);       // return domain id of point rt. return the main domain id if on boundary.
void domain_id_l_aqd2(int *l1,int *l2,double *rt,AQD2 *qd); // calc domain number l1 and l2 at point rt.
void dat_write_aqd2(char *fname,AQD2 *qd);     // write datafile
void dat_read_aqd2(char *fname,AQD2 *qd);      // read datafile 


// -- bem3_aw_qd2_solve_bieq.c --
void solve_bieq_aqd2(AQD2 *qd); // solve boundary integral equations


// -- bem3_aw_qd2_force.c --
void force_FN_aqd2(double *F,double *N,double *rc,int type,AQD2 *qd);
// outputs
// F : radiation force,  F[0]=F_x, F[1]=F_y, F[2]=F_z.
// N : radiation torque, N[0]=N_x, N[1]=N_y, N[2]=N_z.
// inputs
// rc : coordinate of rotation center, rc[0]=x, rc[1]=y, rc[2]=z.
// type :  setting of numerical integration, type=0:4-point GL, type!=0:9-point or 7-point GL.
// qd : pointer of AQD2 object.
 
void surface_area_aqd2(double *S,int did,int type,AQD2 *qd);
// output
// S : surface area of domain id "did".
// intputs
// did : domain id.
// type : setting of numerical integration, type=0:4-point GL, type!=0:9-point or 7-point GL.
// qd : pointer of AQD2 object.


// -- bem3_aw_qd1_field.c --
// velocity potential
int phi_s_aqd2(double complex *phi,double *rt,int type,AQD2 *qd); // scattered field or internal field
int phi_t_aqd2(double complex *phi,double *rt,int type,AQD2 *qd); // total field (add incident field to scattered field)
int phi_i_aqd2(double complex *phi,double *rt,int type,AQD2 *qd); // incidient field
// output
// phi : velocity potential.
// inputs
// rt : coordinate of calculation point, rt[0]=x, rt[1]=y, rt[2]=z.
// type : setting of numerical integration, type=0:4pGL,1:9pGL,2:GLN-point GL,3:GHN-point GL,4:DE. GLN and GHN are defined in d3qd2_const.h.
// qd : point of AQD2 object.
// return domain id of point rt.

// sound pressure
int p_s_aqd2(double complex *p,double *rt,int type,AQD2 *qd); // scattered field
int p_t_aqd2(double complex *p,double *rt,int type,AQD2 *qd); // total field
int p_i_aqd2(double complex *p,double *rt,int type,AQD2 *qd); // incident field
// output
// p : sound pressure.
// others are the same as phi_*_aqd2().

// sound pressure and particle velocity
int pv_s_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *qd); // scattered field
int pv_t_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *qd); // total field
int pv_i_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *qd); // incident field
// output
// p : sound pressure.
// pv: particle velocity, pv[0]=pv_x, pv[1]=pv_y, pv[2]=pv_z.
// others are the same as p_*_aqd2().

// boundary value of velocity potential
void phi_s_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
void phi_t_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
void phi_i_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
// output
// phi : velocity potential.
// inputs
// did : domain id.
// t  : element id. 
// zeta_t, eta_t :  parameter of the calculation point on the element (main domain coordinate), -1 < zeta_t < 1, -1 < eta_t < 1.
// type : setting of numerical integration, type=0:4pGL,1:9pGL,2:GLN-point GL,3:GHN-point GL,4:DE.
// qd : pointer of AQD2 object.

// boundary value of sound pressure and particle velocity 
void pv_s_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
void pv_t_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
void pv_i_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *qd);
// outputs
// p : sound pressure.
// pv: particle velocity, pv[0]=pv_x, pv[1]=pv_y, pv[2]=pv_z.
// others are the same as phi_*_bd_aqd2().


#endif
