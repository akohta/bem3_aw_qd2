/*
 * q0_bil_dcoef_v2_cc.c
 *
 *  Created on: Sep 25, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_bil_dcoef_v2_cc_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double vR[3],aR,i_aR,aW,*vW,rdW,rdv0,rdv1,Wdv0,Wdv1;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  for(i=0;i<9;i++){
    CC[i]=0.0;
    dC0[i]=0.0;
    dC1[i]=0.0;
  }

  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    rdW=vdot_d(vR,vW)*i_aR;
    rdv0=vdot_d(vR,v0)*i_aR;
    rdv1=vdot_d(vR,v1)*i_aR;
    Wdv0=vdot_d(vW,v0);
    Wdv1=vdot_d(vW,v1);

    err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_v2_cc_4p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i  ]=qG*aW; // qG
    CC[i+4]=dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2]; // qH
    CC[8]+=bd->wt_44[i]*rdW*i_aR*i_aR; // F

    dC0[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0); // dF/dv0
    dC1[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1); // dF/dv1
  }
  CC [8]*=I_FP;
  dC0[8]*=I_FP;
  dC1[8]*=I_FP;
}

void q0_bil_dcoef_v2_cc_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,rdv0,rdv1,rdW,Wdv0,Wdv1;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc_9p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc_9p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<9;i++){
      zeta=bd->zt_49[i];
      eta =bd->et_49[i];
      bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      aW=vabs_d(vW);
      rdW=vdot_d(vR,vW)*i_aR;
      rdv0=vdot_d(vR,v0)*i_aR;
      rdv1=vdot_d(vR,v1)*i_aR;
      Wdv0=vdot_d(vW,v0);
      Wdv1=vdot_d(vW,v1);
      Np=bil_Nn(ep,zeta,eta);

      err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_v2_cc_9p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_49[i]*qG*Np*aW; // qG
      CC[ep+4]+=bd->wt_49[i]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np; // qH

      if(ep==0){
        CC[8]+=bd->wt_49[i]*rdW*i_aR*i_aR;
        dC0[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
        dC1[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
      }
    }

    if(ep==0){
      CC [8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_cc_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv0,rdv1,Wdv0,Wdv1,Np,tmpf,tmpdf0,tmpdf1;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){ // init
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg =0.0;
      tmph =0.0;
      tmpf =0.0;
      tmpdf0=0.0;
      tmpdf1=0.0;
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        rdv0=vdot_d(vR,v0)*i_aR;
        rdv1=vdot_d(vR,v1)*i_aR;
        Wdv0=vdot_d(vW,v0);
        Wdv1=vdot_d(vW,v1);
        Np=bil_Nn(ep,zeta,eta);

        err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_v2_cc_GL(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->wli[j]*qG*Np*aW;
        tmph +=bd->wli[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        if(ep==0){
          tmpf +=bd->wli[j]*rdW*i_aR*i_aR;
          tmpdf0+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC [ep  ]+=bd->wli[i]*tmpg;
      CC [ep+4]+=bd->wli[i]*tmph;
      if(ep==0){
        CC[8] +=bd->wli[i]*tmpf;
        dC0[8]+=bd->wli[i]*tmpdf0;
        dC1[8]+=bd->wli[i]*tmpdf1;
      }
    }
    if(ep==0){
      CC[8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_cc_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv0,rdv1,Wdv0,Wdv1,Np,tmpf,tmpdf0,tmpdf1;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){ // init
    CC[l]=0.0;
    dC0[l]=0.0;
    dC1[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg =0.0;
      tmph =0.0;
      tmpf =0.0;
      tmpdf0=0.0;
      tmpdf1=0.0;
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        rdv0=vdot_d(vR,v0)*i_aR;
        rdv1=vdot_d(vR,v1)*i_aR;
        Wdv0=vdot_d(vW,v0);
        Wdv1=vdot_d(vW,v1);
        Np=bil_Nn(ep,zeta,eta);

        err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_v2_cc_GH(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->whi[j]*qG*Np*aW;
        tmph +=bd->whi[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        if(ep==0){
          tmpf +=bd->whi[j]*rdW*i_aR*i_aR;
          tmpdf0+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      if(ep==0){
        CC[8] +=bd->whi[i]*tmpf;
        dC0[8]+=bd->whi[i]*tmpdf0;
        dC1[8]+=bd->whi[i]*tmpdf1;
      }
    }
    if(ep==0){
      CC[8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_cc_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  void q0_bil_dcoef_DE0_cc(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);

  q0_bil_coef_DE(CC,rt,s,k,bd);
  q0_bil_dcoef_DE0_cc(dC0,rt,v0,s,k,bd);
  q0_bil_dcoef_DE0_cc(dC1,rt,v1,s,k,bd);
}

void q0_bil_dcoef_bd_vt2_cc_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  void q0_bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // q0_bil_coef.c
  void q0_bil_dcoef_bd_tz0_cc_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
  void q0_bil_dcoef_bd_te0_cc_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);

  q0_bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  q0_bil_dcoef_bd_tz0_cc_DV(dCdtz,zeta_t,eta_t,s,k,bd);
  q0_bil_dcoef_bd_te0_cc_DV(dCdte,zeta_t,eta_t,s,k,bd);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_bil_dcoef_DE0_cc(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double bil_sdF_zeta(double zeta,void *tmp); // bil_dcoef.c

  TDATA td;
  double err;
  int i,l;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_DE0_cc(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_DE0_cc(), wave number error. k must be a real number. k=%g %+g. Exit...\n",creal(k),cimag(k));
    exit(0);
  }

  td.k=k;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=v[l];
  }

  for(i=0;i<9;i++) dCC[i]=0.0; //init

  // dF/dv
  dCC[8]=I_FP*deintd(bil_sdF_zeta,-1.0,1.0,&td,IEPS,&err);
  if(err<0){ printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_DE0_cc(), bil_sdF_zeta() DE integration error! err=%g. Exit...\n",err);  exit(0);  }
}

void q0_bil_dcoef_bd_tz0_cc_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double q0_bil_coef_bd_F(double zeta_t,double eta_t,int s,double complex k,BOUD *bd);

  double Fp,Fm,cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_bd_tz0_cc_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_bd_tz0_cc_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);
  bil_t_zeta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;
  Fp=q0_bil_coef_bd_F(zeta_t+CDFH,eta_t,s,k,bd);
  Fm=q0_bil_coef_bd_F(zeta_t-CDFH,eta_t,s,k,bd);
  dCC[8]=i_aT*(Fp-Fm)/(2.0*CDFH);
}

void q0_bil_dcoef_bd_te0_cc_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double q0_bil_coef_bd_F(double zeta_t,double eta_t,int s,double complex k,BOUD *bd);

  double Fp,Fm,cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_bd_te0_cc_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_dcoef_bd_te0_cc_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);
  bil_t_eta(vT,zeta_t,eta_t,cr);
  i_aT=1.0/vabs_d(vT);

  for(i=0;i<9;i++) dCC[i]=0.0;
  Fp=q0_bil_coef_bd_F(zeta_t,eta_t+CDFH,s,k,bd);
  Fm=q0_bil_coef_bd_F(zeta_t,eta_t-CDFH,s,k,bd);
  dCC[8]=i_aT*(Fp-Fm)/(2.0*CDFH);
}

double q0_bil_coef_bd_F(double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double bil_sF_theta(double theta,void *tmp); // bil_coef.c

  TDATA td;
  double r[5],th[5],err,bdC,tD,F;
  int m,j;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_coef_bd_F(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_cc.c, q0_bil_coef_bd_F(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.zeta_t=zeta_t;
  td.eta_t = eta_t;
  if(SW_BDR_BIEQ==0){
    td.gn=GLN;
    td.xi=bd->xli;
    td.wi=bd->wli;
  }
  else {
    td.gn=GHN;
    td.xi=bd->xhi;
    td.wi=bd->whi;
  }
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  bil_calc_node_r_th(r,th,zeta_t,eta_t);
  td.qpd=&(bd->qd);

  // F
  if(bil_check_plane(td.cr,td.cw)==1){ // element is plane
    F=0.0;
  }
  else {
    F=0.0;
    bdC=td.cr[0][1]*td.cw[0][2]+td.cr[1][1]*td.cw[1][2]+td.cr[2][1]*td.cw[2][2];
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];
      if(SW_BDT_BIEQ==0){
        tD=0.0;
        for(m=0;m<GHN;m++) tD+=( bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tD*=0.25*(th[j+1]-th[j]);
        F+=tD;
      }
      else {
        F+=deintd(bil_sF_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("q0_bil_coef.c, q0_bil_coef_bd(), bil_sF_theta() DE integration error. err=%f Exit...\n",err);  exit(0);  }
      }
    }
    F*=I_FP*bdC;
  }

  return F;
}
