/*
 * q0_bil_dcoef_v2.c
 *
 *  Created on: Sep 12, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_bil_dcoef_v2_4p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,aW,*vW,rdW,rdv0,rdv1,Wdv0,Wdv1;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  CC[8]=0.0;
  dC0[8]=0.0;
  dC1[8]=0.0;
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
    err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_v2_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i  ]=qG*aW; // qG
    CC[i+4]=dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2]; // qH
    CC[8]+=bd->wt_44[i]*rdW*i_aR*i_aR; // F

    dC0[i  ]=-(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*aW; // dqG/dv0
    dC1[i  ]=-(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*aW; // dqG/dv1
    dC0[i+4]=-(d2qG[0]*vW[0]*v0[0]+d2qG[1]*vW[1]*v0[1]+d2qG[2]*vW[2]*v0[2]
               +d2qG[3]*(vW[0]*v0[1]+vW[1]*v0[0])
               +d2qG[4]*(vW[1]*v0[2]+vW[2]*v0[1])
               +d2qG[5]*(vW[2]*v0[0]+vW[0]*v0[2])); // dqH/dv0
    dC1[i+4]=-(d2qG[0]*vW[0]*v1[0]+d2qG[1]*vW[1]*v1[1]+d2qG[2]*vW[2]*v1[2]
               +d2qG[3]*(vW[0]*v1[1]+vW[1]*v1[0])
               +d2qG[4]*(vW[1]*v1[2]+vW[2]*v1[1])
               +d2qG[5]*(vW[2]*v1[0]+vW[0]*v1[2])); // dqH/dv1
    dC0[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0); // dF/dv0
    dC1[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1); // dF/dv1
  }
  CC [8]*=I_FP;
  dC0[8]*=I_FP;
  dC1[8]*=I_FP;
}

void q0_bil_dcoef_v2_9p(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,rdv0,rdv1,rdW,Wdv0,Wdv1;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_9p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_9p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
      err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_v2_9p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_49[i]*qG*Np*aW; // qG
      CC[ep+4]+=bd->wt_49[i]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np; // qH

      dC0[ep  ]+=bd->wt_49[i]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*aW*Np; // dqG/dv0
      dC1[ep  ]+=bd->wt_49[i]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*aW*Np; // dqG/dv1
      dC0[ep+4]+=bd->wt_49[i]*( d2qG[0]*vW[0]*v0[0]+d2qG[1]*vW[1]*v0[1]+d2qG[2]*vW[2]*v0[2]
                               +d2qG[3]*(vW[0]*v0[1]+vW[1]*v0[0])
                               +d2qG[4]*(vW[1]*v0[2]+vW[2]*v0[1])
                               +d2qG[5]*(vW[2]*v0[0]+vW[0]*v0[2]))*Np; // dqH/dv0
      dC1[ep+4]+=bd->wt_49[i]*( d2qG[0]*vW[0]*v1[0]+d2qG[1]*vW[1]*v1[1]+d2qG[2]*vW[2]*v1[2]
                               +d2qG[3]*(vW[0]*v1[1]+vW[1]*v1[0])
                               +d2qG[4]*(vW[1]*v1[2]+vW[2]*v1[1])
                               +d2qG[5]*(vW[2]*v1[0]+vW[0]*v1[2]))*Np; // dqH/dv1
      if(ep==0){
        CC[8]+=bd->wt_49[i]*rdW*i_aR*i_aR;
        dC0[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
        dC1[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
      }
    }
    dC0[ep  ]*=-1.0;
    dC1[ep  ]*=-1.0;
    dC0[ep+4]*=-1.0;
    dC1[ep+4]*=-1.0;
    if(ep==0){
      CC [8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_GL(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6],tmpg,tmph,tmpdg0,tmpdg1,tmpdh0,tmpdh1;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv0,rdv1,Wdv0,Wdv1,Np,tmpf,tmpdf0,tmpdf1;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
      tmpdg0=0.0;
      tmpdg1=0.0;
      tmph =0.0;
      tmpdh0=0.0;
      tmpdh1=0.0;
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
        err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_v2_GL(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->wli[j]*qG*Np*aW;
        tmpdg0+=bd->wli[j]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*aW*Np;
        tmpdg1+=bd->wli[j]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*aW*Np;
        tmph  +=bd->wli[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        tmpdh0+=bd->wli[j]*( d2qG[0]*vW[0]*v0[0]+d2qG[1]*vW[1]*v0[1]+d2qG[2]*vW[2]*v0[2]
                            +d2qG[3]*(vW[0]*v0[1]+vW[1]*v0[0])
                            +d2qG[4]*(vW[1]*v0[2]+vW[2]*v0[1])
                            +d2qG[5]*(vW[2]*v0[0]+vW[0]*v0[2]))*Np;
        tmpdh1+=bd->wli[j]*( d2qG[0]*vW[0]*v1[0]+d2qG[1]*vW[1]*v1[1]+d2qG[2]*vW[2]*v1[2]
                            +d2qG[3]*(vW[0]*v1[1]+vW[1]*v1[0])
                            +d2qG[4]*(vW[1]*v1[2]+vW[2]*v1[1])
                            +d2qG[5]*(vW[2]*v1[0]+vW[0]*v1[2]))*Np;
        if(ep==0){
          tmpf +=bd->wli[j]*rdW*i_aR*i_aR;
          tmpdf0+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC [ep  ]+=bd->wli[i]*tmpg;
      CC [ep+4]+=bd->wli[i]*tmph;
      dC0[ep  ]+=bd->wli[i]*tmpdg0;
      dC1[ep  ]+=bd->wli[i]*tmpdg1;
      dC0[ep+4]+=bd->wli[i]*tmpdh0;
      dC1[ep+4]+=bd->wli[i]*tmpdh1;
      if(ep==0){
        CC[8] +=bd->wli[i]*tmpf;
        dC0[8]+=bd->wli[i]*tmpdf0;
        dC1[8]+=bd->wli[i]*tmpdf1;
      }
    }
    dC0[ep  ]*=-1.0;
    dC1[ep  ]*=-1.0;
    dC0[ep+4]*=-1.0;
    dC1[ep+4]*=-1.0;
    if(ep==0){
      CC[8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_GH(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6],tmpg,tmph,tmpdg0,tmpdg1,tmpdh0,tmpdh1;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv0,rdv1,Wdv0,Wdv1,Np,tmpf,tmpdf0,tmpdf1;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
      tmpdg0=0.0;
      tmpdg1=0.0;
      tmph =0.0;
      tmpdh0=0.0;
      tmpdh1=0.0;
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
        err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_v2_GH(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->whi[j]*qG*Np*aW;
        tmpdg0+=bd->whi[j]*(dqG[0]*v0[0]+dqG[1]*v0[1]+dqG[2]*v0[2])*aW*Np;
        tmpdg1+=bd->whi[j]*(dqG[0]*v1[0]+dqG[1]*v1[1]+dqG[2]*v1[2])*aW*Np;
        tmph  +=bd->whi[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        tmpdh0+=bd->whi[j]*( d2qG[0]*vW[0]*v0[0]+d2qG[1]*vW[1]*v0[1]+d2qG[2]*vW[2]*v0[2]
                            +d2qG[3]*(vW[0]*v0[1]+vW[1]*v0[0])
                            +d2qG[4]*(vW[1]*v0[2]+vW[2]*v0[1])
                            +d2qG[5]*(vW[2]*v0[0]+vW[0]*v0[2]))*Np;
        tmpdh1+=bd->whi[j]*( d2qG[0]*vW[0]*v1[0]+d2qG[1]*vW[1]*v1[1]+d2qG[2]*vW[2]*v1[2]
                            +d2qG[3]*(vW[0]*v1[1]+vW[1]*v1[0])
                            +d2qG[4]*(vW[1]*v1[2]+vW[2]*v1[1])
                            +d2qG[5]*(vW[2]*v1[0]+vW[0]*v1[2]))*Np;
        if(ep==0){
          tmpf +=bd->whi[j]*rdW*i_aR*i_aR;
          tmpdf0+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv0*rdW-Wdv0);
          tmpdf1+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv1*rdW-Wdv1);
        }
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      dC0[ep  ]+=bd->whi[i]*tmpdg0;
      dC1[ep  ]+=bd->whi[i]*tmpdg1;
      dC0[ep+4]+=bd->whi[i]*tmpdh0;
      dC1[ep+4]+=bd->whi[i]*tmpdh1;
      if(ep==0){
        CC[8] +=bd->whi[i]*tmpf;
        dC0[8]+=bd->whi[i]*tmpdf0;
        dC1[8]+=bd->whi[i]*tmpdf1;
      }
    }
    dC0[ep  ]*=-1.0;
    dC1[ep  ]*=-1.0;
    dC0[ep+4]*=-1.0;
    dC1[ep+4]*=-1.0;
    if(ep==0){
      CC[8]*=I_FP;
      dC0[8]*=I_FP;
      dC1[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_v2_DE(double complex *CC,double complex *dC0,double complex *dC1,double *rt,double *v0,double *v1,int s,double complex k,BOUD *bd)
{
  void q0_bil_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd); // q0_bil_dcoef.c

  q0_bil_coef_DE(CC,rt,s,k,bd);
  q0_bil_dcoef_DE0(dC0,rt,v0,s,k,bd);
  q0_bil_dcoef_DE0(dC1,rt,v1,s,k,bd);
}

void q0_bil_dcoef_bd_vt2_DV(double complex *CC,double complex *dCdtz,double complex *dCdte,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  void q0_bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd); // q0_bil_coef.c
  void q0_bil_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);
  void q0_bil_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd);

  q0_bil_coef_bd(CC,zeta_t,eta_t,s,k,bd);
  q0_bil_dcoef_bd_tz0_DV(dCdtz,zeta_t,eta_t,s,k,bd);
  q0_bil_dcoef_bd_te0_DV(dCdte,zeta_t,eta_t,s,k,bd);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_bil_dcoef_bd_tz0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_bd_tz0_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_bd_tz0_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);
  bil_t_zeta(vT,zeta_t,eta_t,cr);
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

void q0_bil_dcoef_bd_te0_DV(double complex *dCC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex Cp[9],Cm[9];
  double cr[3][4],cw[3][3],vT[3],i_aT;
  int i;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_bd_te0_DV() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_v2.c, q0_bil_dcoef_bd_te0_DV() wave number error. k must be same as open region wave number. k=%g %+g. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);
  bil_t_eta(vT,zeta_t,eta_t,cr);
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
