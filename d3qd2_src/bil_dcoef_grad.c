/*
 * bil_dcoef_grad.c
 *
 *  Created on: Sep 26, 2019
 *      Author: ohta
 */
#include "d3qd2_elem.h"

void bil_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb;
  double vR[3],aR,i_aR,aW,tD,*vW,rdW;
  int sig,i,l,j;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  CC[8]=0.0;
  for(j=0;j<3;j++)
    for(i=0;i<9;i++) dC[j][i]=0.0;

  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    ca=I*k*aR;
    cb=ca-1.0;
    if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
    else ce=cexp(ca);
    rdW=vdot_d(vR,vW)*i_aR;
    tD=rdW*i_aR*i_aR;

    CC[i  ]=I_FP*ce*i_aR*aW; // G
    CC[i+4]=I_FP*(ca-1.0)*ce*tD; // H
    CC[8]+=bd->wt_44[i]*tD; // F

    for(j=0;j<3;j++){
      dC[j][i  ]=-I_FP*cb*ce*i_aR*i_aR*vR[j]*i_aR*aW;
      dC[j][i+4]= I_FP*(ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*vR[j]*i_aR*rdW-cb*vW[j]));
      dC[j][8  ]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*vR[j]*i_aR*rdW-vW[j]); // dF/dv
    }
  }

  CC[8]*=I_FP;
  for(j=0;j<3;j++){
    dC[j][8]*=I_FP;
  }

  if(sig==1){
    bil_convert_CC(CC);
    for(j=0;j<3;j++) bil_convert_dCC(dC[j]);
  }

}

void bil_dcoef_grad_9p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,rdW;
  int sig,i,l,ep,j;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC [l]=0.0;
  for(j=0;j<3;j++)
    for(l=0;l<9;l++) dC[j][l]=0.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<9;i++){
      zeta=bd->zt_49[i];
      eta =bd->et_49[i];
      bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      aW=vabs_d(vW);
      ca=I*k*aR;
      cb=ca-1.0;
      if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
      else ce=cexp(ca);
      rdW=vdot_d(vR,vW)*i_aR;
      tD=rdW*i_aR*i_aR;
      Np=bil_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_49[i]*ce*i_aR*aW*Np;
      CC[ep+4]+=bd->wt_49[i]*(ca-1.0)*ce*tD*Np;

      for(j=0;j<3;j++){
        dC[j][ep  ]+=bd->wt_49[i]*cb*ce*i_aR*i_aR*vR[j]*i_aR*aW*Np;
        dC[j][ep+4]+=bd->wt_49[i]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*vR[j]*i_aR*rdW-cb*vW[j])*Np;
      }
      if(ep==0){
        CC [8]+=bd->wt_49[i]*tD;
        for(j=0;j<3;j++) dC[j][8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*vR[j]*i_aR*rdW-vW[j]);
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    for(j=0;j<3;j++){
      dC[j][ep  ]*=-I_FP;
      dC[j][ep+4]*= I_FP;
    }
    if(ep==0){
      CC [8]*=I_FP;
      for(j=0;j<3;j++) dC[j][8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    for(j=0;j<3;j++) bil_convert_dCC(dC[j]);
  }
}

void bil_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg[3],tmpdh[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf[3],rdW;
  int sig,i,j,l,ep,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC [l]=0.0;
  for(p=0;p<3;p++)
    for(l=0;l<9;l++) dC[p][l]=0.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf =0.0;
      for(p=0;p<3;p++){
        tmpdg[p]=0.0;
        tmpdh[p]=0.0;
        tmpdf[p]=0.0;
      }
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        cb=ca-1.0;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        rdW=vdot_d(vR,vW)*i_aR;
        tD=rdW*i_aR*i_aR;
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->wli[j]*ce*i_aR*aW*Np;
        tmph+=bd->wli[j]*(ca-1.0)*ce*tD*Np;

        for(p=0;p<3;p++){
          tmpdg[p]+=bd->wli[j]*cb*ce*i_aR*i_aR*vR[p]*i_aR*aW*Np;
          tmpdh[p]+=bd->wli[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*vR[p]*i_aR*rdW-cb*vW[p])*Np;
        }
        if(ep==0){
          tmpf +=bd->wli[j]*tD;
          for(p=0;p<3;p++) tmpdf[p]+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]);
        }
      }
      CC [ep  ]+=bd->wli[i]*tmpg;
      CC [ep+4]+=bd->wli[i]*tmph;
      for(p=0;p<3;p++){
        dC[p][ep  ]+=bd->wli[i]*tmpdg[p];
        dC[p][ep+4]+=bd->wli[i]*tmpdh[p];
      }
      if(ep==0){
        CC [8]+=bd->wli[i]*tmpf;
        for(p=0;p<3;p++) dC[p][8]+=bd->wli[i]*tmpdf[p];
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-I_FP;
      dC[p][ep+4]*= I_FP;
    }
    if(ep==0){
      CC [8]*=I_FP;
      for(p=0;p<3;p++) dC[p][8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    for(p=0;p<3;p++) bil_convert_dCC(dC[p]);
  }
}

void bil_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,ce,cb,tmpdg[3],tmpdh[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,tD,zeta,eta,cr[3][4],cw[3][3],Np,tmpf,tmpdf[3],rdW;
  int sig,i,j,l,ep,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC [l]=0.0;
  for(p=0;p<3;p++)
    for(l=0;l<9;l++) dC[p][l]=0.0;

  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf=0.0;
      for(p=0;p<3;p++){
        tmpdg[p]=0.0;
        tmpdh[p]=0.0;
        tmpdf[p]=0.0;
      }
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        ca=I*k*aR;
        cb=ca-1.0;
        if(cimag(k)==0.0) ce=cos(creal(k)*aR)+I*sin(creal(k)*aR);
        else ce=cexp(ca);
        rdW=vdot_d(vR,vW)*i_aR;
        tD=rdW*i_aR*i_aR;
        Np=bil_Nn(ep,zeta,eta);

        tmpg+=bd->whi[j]*ce*i_aR*aW*Np;
        tmph+=bd->whi[j]*(ca-1.0)*ce*tD*Np;

        for(p=0;p<3;p++){
          tmpdg[p]+=bd->whi[j]*cb*ce*i_aR*i_aR*vR[p]*i_aR*aW*Np;
          tmpdh[p]+=bd->whi[j]*ce*i_aR*i_aR*i_aR*((k*k*aR*aR+3.0*cb)*vR[p]*i_aR*rdW-cb*vW[p])*Np;
        }
        if(ep==0){
          tmpf +=bd->whi[j]*tD;
          for(p=0;p<3;p++) tmpdf[p]+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*vR[p]*i_aR*rdW-vW[p]);
        }
      }
      CC [ep  ]+=bd->whi[i]*tmpg;
      CC [ep+4]+=bd->whi[i]*tmph;
      for(p=0;p<3;p++){
        dC[p][ep  ]+=bd->whi[i]*tmpdg[p];
        dC[p][ep+4]+=bd->whi[i]*tmpdh[p];
      }
      if(ep==0){
        CC [8]+=bd->whi[i]*tmpf;
        for(p=0;p<3;p++) dC[p][8]+=bd->whi[i]*tmpdf[p];
      }
    }
    CC [ep  ]*= I_FP;
    CC [ep+4]*= I_FP;
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-I_FP;
      dC[p][ep+4]*= I_FP;
    }
    if(ep==0){
      CC [8]*=I_FP;
      for(p=0;p<3;p++) dC[p][8]*=I_FP;
    }
  }

  if(sig==1){
    bil_convert_CC(CC);
    for(p=0;p<3;p++) bil_convert_dCC(dC[p]);
  }
}
