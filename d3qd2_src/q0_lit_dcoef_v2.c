/*
 * q0_lit_dcoef_v2.c
 *
 *  Created on: Sep 18, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_lit_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,aV,*vV,rdV,rdv0,rdv1,Vdv0,Vdv1,Nn;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  for(i=0;i<9;i++){
    CC [i]=0.0;
    dC0[i]=0.0;
    dC1[i]=0.0;
  }
  vV=bd->wen[s][0];
  aV=vabs_d(vV);
  Vdv0=vdot_d(vV,v0);
  Vdv1=vdot_d(vV,v1);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    rdV=vdot_d(vR,vV)*i_aR;
    rdv0=vdot_d(vR,v0)*i_aR;
    rdv1=vdot_d(vR,v1)*i_aR;
    err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_v2_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i+0]=bd->wt_34[i]*qG;
    CC[i+4]=bd->wt_34[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2]);
    CC[8]+=bd->wt_34[i]*rdV*i_aR*i_aR;

    dC0[i+0]=bd->wt_34[i]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2]);
    dC0[i+4]=bd->wt_34[i]*( d2qG[0]*vV[0]*v0[0]+d2qG[1]*vV[1]*v0[1]+d2qG[2]*vV[2]*v0[2]
                           +d2qG[3]*(vV[0]*v0[1]+vV[1]*v0[0])
                           +d2qG[4]*(vV[1]*v0[2]+vV[2]*v0[1])
                           +d2qG[5]*(vV[2]*v0[0]+vV[0]*v0[2]));
    dC0[8 ]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdV-Vdv0 );

    dC1[i+0]=bd->wt_34[i]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2]);
    dC1[i+4]=bd->wt_34[i]*( d2qG[0]*vV[0]*v1[0]+d2qG[1]*vV[1]*v1[1]+d2qG[2]*vV[2]*v1[2]
                           +d2qG[3]*(vV[0]*v1[1]+vV[1]*v1[0])
                           +d2qG[4]*(vV[1]*v1[2]+vV[2]*v1[1])
                           +d2qG[5]*(vV[2]*v1[0]+vV[0]*v1[2]));
    dC1[8 ]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdV-Vdv1 );
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  rdv0=vdot_d(vR,v0)*i_aR;
  rdv1=vdot_d(vR,v1)*i_aR;
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
  if(err<0){
    printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_v2_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*qG*Nn;
    CC[i+4]+=bd->wt_34[3]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;

    dC0[i+0]+=bd->wt_34[3]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*Nn;
    dC0[i+4]+=bd->wt_34[3]*( d2qG[0]*vV[0]*v0[0]+d2qG[1]*vV[1]*v0[1]+d2qG[2]*vV[2]*v0[2]
                            +d2qG[3]*(vV[0]*v0[1]+vV[1]*v0[0])
                            +d2qG[4]*(vV[1]*v0[2]+vV[2]*v0[1])
                            +d2qG[5]*(vV[2]*v0[0]+vV[0]*v0[2]))*Nn;

    dC1[i+0]+=bd->wt_34[3]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*Nn;
    dC1[i+4]+=bd->wt_34[3]*( d2qG[0]*vV[0]*v1[0]+d2qG[1]*vV[1]*v1[1]+d2qG[2]*vV[2]*v1[2]
                            +d2qG[3]*(vV[0]*v1[1]+vV[1]*v1[0])
                            +d2qG[4]*(vV[1]*v1[2]+vV[2]*v1[1])
                            +d2qG[5]*(vV[2]*v1[0]+vV[0]*v1[2]))*Nn;
  }
  CC [8]+=bd->wt_34[3]*rdV*i_aR*i_aR;
  dC0[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdV-Vdv0 );
  dC1[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdV-Vdv1 );

  for(i=0;i<3;i++){
    CC[i+0]*=0.5*aV;
    CC[i+4]*=0.5;
    dC0[i+0]*=-0.5*aV;
    dC0[i+4]*=-0.5;
    dC1[i+0]*=-0.5*aV;
    dC1[i+4]*=-0.5;
  }
  CC[8]*=I_EP;
  dC0[8]*=I_EP;
  dC1[8]*=I_EP;
}

void q0_lit_dcoef_v2_7p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double cr[3][4],cw[3][3],zeta,eta,*vV,aV,vR[3],aR,i_aR,rdV,rdv0,rdv1,Vdv0,Vdv1,Nn;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2_7p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2_7p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  vV=bd->wen[s][0];
  aV=vabs_d(vV);
  Vdv0=vdot_d(vV,v0);
  Vdv1=vdot_d(vV,v1);

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      rdV=vdot_d(vR,vV)*i_aR;
      rdv0=vdot_d(vR,v0)*i_aR;
      rdv1=vdot_d(vR,v1)*i_aR;
      Nn=lit_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_v2_7p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC [ep  ]+=bd->wt_37[i]*qG*Nn;
      CC [ep+4]+=bd->wt_37[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;
      dC0[ep+0]+=bd->wt_37[i]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*Nn;
      dC0[ep+4]+=bd->wt_37[i]*( d2qG[0]*vV[0]*v0[0]+d2qG[1]*vV[1]*v0[1]+d2qG[2]*vV[2]*v0[2]
                               +d2qG[3]*(vV[0]*v0[1]+vV[1]*v0[0])
                               +d2qG[4]*(vV[1]*v0[2]+vV[2]*v0[1])
                               +d2qG[5]*(vV[2]*v0[0]+vV[0]*v0[2]))*Nn;
      dC1[ep+0]+=bd->wt_37[i]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*Nn;
      dC1[ep+4]+=bd->wt_37[i]*( d2qG[0]*vV[0]*v1[0]+d2qG[1]*vV[1]*v1[1]+d2qG[2]*vV[2]*v1[2]
                               +d2qG[3]*(vV[0]*v1[1]+vV[1]*v1[0])
                               +d2qG[4]*(vV[1]*v1[2]+vV[2]*v1[1])
                               +d2qG[5]*(vV[2]*v1[0]+vV[0]*v1[2]))*Nn;
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdV*i_aR*i_aR;
        dC0[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv0*rdV-Vdv0 );
        dC1[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv1*rdV-Vdv1 );
      }
    }
    CC[ep+0]*=0.5*aV;
    CC[ep+4]*=0.5;
    dC0[ep+0]*=-0.5*aV;
    dC0[ep+4]*=-0.5;
    dC1[ep+0]*=-0.5*aV;
    dC1[ep+4]*=-0.5;
    if(ep==0){
      CC [8]*=I_EP;
      dC0[8]*=I_EP;
      dC1[8]*=I_EP;
    }
  }
}

void q0_lit_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_v2_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[9],tpj[9],tpi[9];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.qpd=&(bd->qd);
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,td.cw);
  td.aW=vabs_d(vV);
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=v0[l];
    td.vu[l]=v1[l];
  }

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<9;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<9;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];

        for(l=0;l<9;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC[ep  ]+=tpi[0]*bd->wli[i];
      CC[ep+4]+=tpi[1]*bd->wli[i];
      dC0[ep  ]+=tpi[3]*bd->wli[i];
      dC0[ep+4]+=tpi[4]*bd->wli[i];
      dC1[ep  ]+=tpi[6]*bd->wli[i];
      dC1[ep+4]+=tpi[7]*bd->wli[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        dC0[8]+=tpi[5]*bd->wli[i];
        dC1[8]+=tpi[8]*bd->wli[i];
      }
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    dC0[ep  ]*=-1.0/24.0*td.aW;
    dC0[ep+4]*=-1.0/24.0;
    dC1[ep  ]*=-1.0/24.0*td.aW;
    dC1[ep+4]*=-1.0/24.0;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dC0[8]*=0.0625*SQ3*I_SPSQ3;
      dC1[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}

void q0_lit_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_v2_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[9],tpj[9],tpi[9];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.qpd=&(bd->qd);
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,td.cw);
  td.aW=vabs_d(vV);
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=v0[l];
    td.vu[l]=v1[l];
  }

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<9;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<9;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_v2_zeta_eta(ret,&td);
        for(l=0;l<9;l++) tpj[l]+=ret[l];

        for(l=0;l<9;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC [ep  ]+=tpi[0]*bd->whi[i];
      CC [ep+4]+=tpi[1]*bd->whi[i];
      dC0[ep  ]+=tpi[3]*bd->whi[i];
      dC0[ep+4]+=tpi[4]*bd->whi[i];
      dC1[ep  ]+=tpi[6]*bd->whi[i];
      dC1[ep+4]+=tpi[7]*bd->whi[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->whi[i];
        dC0[8]+=tpi[5]*bd->whi[i];
        dC1[8]+=tpi[8]*bd->whi[i];
      }
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    dC0[ep  ]*=-1.0/24.0*td.aW;
    dC0[ep+4]*=-1.0/24.0;
    dC1[ep  ]*=-1.0/24.0*td.aW;
    dC1[ep+4]*=-1.0/24.0;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dC0[8]*=0.0625*SQ3*I_SPSQ3;
      dC1[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}

void q0_lit_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  void q0_lit_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd); // q0_lit_dcoef.c

  q0_lit_coef_DE(CC,rt,s,k,bd);
  q0_lit_dcoef_DE0(dC0,rt,v0,s,k,bd);
  q0_lit_dcoef_DE0(dC1,rt,v1,s,k,bd);
}

void q0_lit_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  void q0_lit_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // q0_lit_coef.c
  void q0_lit_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
  void q0_lit_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);

  q0_lit_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  q0_lit_dcoef_bd_tz0_DV(dCdtz,zeta_t,eta_t,s,k,bd);
  q0_lit_dcoef_bd_te0_DV(dCdte,zeta_t,eta_t,s,k,bd);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_lit_sdqGHF_v2_zeta_eta(double complex *ret,TDATA *td)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],vV[3],aR,i_aR,Np,rdV,rdv0,rdv1,Vdv0,Vdv1;
  int l,err;

  vV[0]=td->cw[0][0];  vV[1]=td->cw[1][0];  vV[2]=td->cw[2][0];
  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);
  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  rdv0=vdot_d(vR,td->vt)*i_aR;
  rdv1=vdot_d(vR,td->vu)*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_dcoef_v2.c, q0_lit_sdqGHF_v2_zeta_eta(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  ret[0]=qG*Np;
  ret[1]=(dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*Np;

  ret[3]=(dqG[0]*td->vt[0]+dqG[1]*td->vt[1]+dqG[2]*td->vt[2])*Np;
  ret[4]=( d2qG[0]*vV[0]*td->vt[0]+d2qG[1]*vV[1]*td->vt[1]+d2qG[2]*vV[2]*td->vt[2]
          +d2qG[3]*(vV[0]*td->vt[1]+vV[1]*td->vt[0])
          +d2qG[4]*(vV[1]*td->vt[2]+vV[2]*td->vt[1])
          +d2qG[5]*(vV[2]*td->vt[0]+vV[0]*td->vt[2]))*Np;

  ret[6]=(dqG[0]*td->vu[0]+dqG[1]*td->vu[1]+dqG[2]*td->vu[2])*Np;
  ret[7]=( d2qG[0]*vV[0]*td->vu[0]+d2qG[1]*vV[1]*td->vu[1]+d2qG[2]*vV[2]*td->vu[2]
          +d2qG[3]*(vV[0]*td->vu[1]+vV[1]*td->vu[0])
          +d2qG[4]*(vV[1]*td->vu[2]+vV[2]*td->vu[1])
          +d2qG[5]*(vV[2]*td->vu[0]+vV[0]*td->vu[2]))*Np;

  if(td->nid==0){
    ret[2]=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
    Vdv0=vdot_d(vV,td->vt);
    Vdv1=vdot_d(vV,td->vu);
    ret[5]=i_aR*i_aR*i_aR*( 3.0*rdv0*rdV-Vdv0 );
    ret[8]=i_aR*i_aR*i_aR*( 3.0*rdv1*rdV-Vdv1 );
  }
  else {
    ret[2]=0.0;
    ret[5]=0.0;
    ret[8]=0.0;
  }
}

void q0_lit_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_bd_tz0_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_bd_tz0_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_t_zeta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;
  q0_bil_coef_bd(Cp,zeta_t+CDFH,eta_t,s,k,bd);
  q0_bil_coef_bd(Cm,zeta_t-CDFH,eta_t,s,k,bd);
  for(i=0;i<4;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }
  dCC[8]=i_aT*(Cp[8]-Cm[8])/(2.0*CDFH);
}

void q0_lit_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_bd_te0_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_v2.c, q0_lit_dcoef_bd_te0_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_t_eta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;
  q0_bil_coef_bd(Cp,zeta_t,eta_t+CDFH,s,k,bd);
  q0_bil_coef_bd(Cm,zeta_t,eta_t-CDFH,s,k,bd);
  for(i=0;i<4;i++){
    dCC[i+0]=i_aT*(Cp[i+0]-Cm[i+0])/(2.0*CDFH);
    dCC[i+4]=i_aT*(Cp[i+4]-Cm[i+4])/(2.0*CDFH);
  }
  dCC[8]=i_aT*(Cp[8]-Cm[8])/(2.0*CDFH);
}
