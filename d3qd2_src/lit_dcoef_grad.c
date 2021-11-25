/*
 * lit_dcoef_grad.c
 *
 *  Created on: Sep 26, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void lit_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double vR[3],*vW,aR,i_aR,aW,rdW,Nn;
  int sig,i,l,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  for(i=0;i<9;i++) CC [i]=0.0;
  for(p=0;p<3;p++)
    for(i=0;i<9;i++) dC[p][i]=0.0;

  vW=bd->wen[s][0];
  aW=vabs_d(vW);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    kR=k*aR;
    ca=I*kR;
    cb=ca-1.0;
    if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
    else ce=cexp(ca);
    rdW=vdot_d(vR,vW)*i_aR;

    CC[i+0]=bd->wt_34[i]*ce*i_aR;
    CC[i+4]=bd->wt_34[i]*(ca-1.0)*ce*rdW*i_aR*i_aR;
    CC [8]+=bd->wt_34[i]*rdW*i_aR*i_aR; // F

    for(p=0;p<3;p++){
      dC[p][i+0]=bd->wt_34[i]*cb*ce*i_aR*i_aR*vR[p]*i_aR; // dG/dv
      dC[p][i+4]=bd->wt_34[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*vR[p]*i_aR*rdW-cb*vW[p] ); // dH/dv
      dC[p][8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdW-vW[p]); // dF/dv
    }
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  kR=k*aR;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  rdW=vdot_d(vR,vW)*i_aR;

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*ce*i_aR*Nn;
    CC[i+4]+=bd->wt_34[3]*(ca-1.0)*ce*rdW*i_aR*i_aR*Nn;
    for(p=0;p<3;p++){
      dC[p][i+0]+=bd->wt_34[3]*cb*ce*i_aR*i_aR*vR[p]*i_aR*Nn;
      dC[p][i+4]+=bd->wt_34[3]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*vR[p]*i_aR*rdW-cb*vW[p])*Nn;
    }
    if(i==0){
      CC [8]+=bd->wt_34[3]*rdW*i_aR*i_aR;
      for(p=0;p<3;p++) dC[p][8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdW-vW[p]);
    }
  }

  for(i=0;i<3;i++){
    CC[i+0]*=I_EP*aW;
    CC[i+4]*=I_EP;
    for(p=0;p<3;p++){
      dC[p][i+0]*=-I_EP*aW;
      dC[p][i+4]*= I_EP;
    }
    if(i==0){
      CC [8]*=I_EP;
      for(p=0;p<3;p++) dC[p][8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    for(p=0;p<3;p++) lit_convert_dCC(dC[p]);
  }
}

void lit_dcoef_grad_7p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex ca,cb,ce,kR;
  double cr[3][4],cw[3][3],zeta,eta,vW[3],aW,vR[3],aR,i_aR,rdW,Np;
  int sig,i,l,ep,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,cw);
  aW=vabs_d(vW);

  for(l=0;l<9;l++) CC [l]=0.0;
  for(p=0;p<3;p++)
    for(l=0;l<9;l++) dC[p][l]=0.0;

  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      kR=k*aR;
      ca=I*kR;
      cb=ca-1.0;
      if(cimag(k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
      else ce=cexp(ca);
      rdW=vdot_d(vR,vW)*i_aR;
      Np=lit_Nn(ep,zeta,eta);

      CC[ep  ]+=bd->wt_37[i]*ce*i_aR*Np;
      CC[ep+4]+=bd->wt_37[i]*(ca-1.0)*ce*rdW*i_aR*i_aR*Np;
      for(p=0;p<3;p++){
        dC[p][ep  ]+=bd->wt_37[i]*cb*ce*i_aR*i_aR*vR[p]*i_aR*Np;
        dC[p][ep+4]+=bd->wt_37[i]*ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*vR[p]*i_aR*rdW-cb*vW[p])*Np;
      }
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdW*i_aR*i_aR;
        for(p=0;p<3;p++) dC[p][8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdW-vW[p]);
      }
    }
    CC [ep  ]*= I_EP*aW;
    CC [ep+4]*= I_EP;
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-I_EP*aW;
      dC[p][ep+4]*= I_EP;
    }
    if(ep==0){
      CC [8]*=I_EP;
      for(p=0;p<3;p++) dC[p][8]*=I_EP;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    for(p=0;p<3;p++) lit_convert_dCC(dC[p]);
  }
}

void lit_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  void lit_sdGL_grad_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[12],tpj[12],tpi[12];
  double vW[3],a,b,wt;
  int sig,i,j,l,ep,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  td.aW=vabs_d(vW);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC [l]=0.0;
  for(p=0;p<3;p++)
    for(l=0;l<9;l++) dC[p][l]=0.0;

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<12;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<12;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];

        for(l=0;l<12;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC [ep  ]+=tpi[0]*bd->wli[i];
      CC [ep+4]+=tpi[1]*bd->wli[i];
      for(p=0;p<3;p++){
        dC[p][ep  ]+=tpi[3+p]*bd->wli[i];
        dC[p][ep+4]+=tpi[6+p]*bd->wli[i];
      }
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        for(p=0;p<3;p++) dC[p][8]+=tpi[9+p]*bd->wli[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
      dC[p][ep+4]*= 0.0625*SQ3*I_SPSQ3;
    }
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      for(p=0;p<3;p++) dC[p][8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    for(p=0;p<3;p++) lit_convert_dCC(dC[p]);
  }
}


void lit_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  void lit_sdGL_grad_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[12],tpj[12],tpi[12];
  double vW[3],a,b,wt;
  int sig,i,j,l,ep,p;

  if(s>0) sig=0;
  else {
    s=-s;
    sig=1;
  }

  td.k=k;
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vW,0.0,0.0,td.cw);
  td.aW=vabs_d(vW);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC [l]=0.0;
  for(p=0;p<3;p++)
    for(l=0;l<9;l++) dC[p][l]=0.0;

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<12;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<12;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        lit_sdGL_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];

        for(l=0;l<12;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC [ep  ]+=tpi[0]*bd->whi[i];
      CC [ep+4]+=tpi[1]*bd->whi[i];
      for(p=0;p<3;p++){
        dC[p][ep  ]+=tpi[3+p]*bd->whi[i];
        dC[p][ep+4]+=tpi[6+p]*bd->whi[i];
      }
      if(ep==0){
        CC [8]+=tpi[2]*bd->whi[i];
        for(p=0;p<3;p++) dC[p][8]+=tpi[9+p]*bd->whi[i];
      }
    }
    CC[ep  ]*=0.0625*SQ3*I_SPSQ3*td.aW;
    CC[ep+4]*=0.0625*SQ3*I_SPSQ3;
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-0.0625*SQ3*I_SPSQ3*td.aW;
      dC[p][ep+4]*= 0.0625*SQ3*I_SPSQ3;
    }
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      for(p=0;p<3;p++) dC[p][8]*=0.0625*SQ3*I_SPSQ3;
    }
  }

  if(sig==1){
    lit_convert_CC(CC);
    for(p=0;p<3;p++) lit_convert_dCC(dC[p]);
  }
}

//////////////////////////////////////////////////////////////////
void lit_sdGL_grad_zeta_eta(double complex *ret,TDATA *td)
{
  double complex ca,cb,ce,kR;
  double vR[3],aR,i_aR,rdW,Np;
  int l,p;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  kR=td->k*aR;
  ca=I*kR;
  cb=ca-1.0;
  if(cimag(td->k)==0.0) ce=cos(creal(kR))+I*sin(creal(kR));
  else ce=cexp(ca);
  rdW=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);

  ret[0]=ce*i_aR*Np; // G
  ret[1]=(ca-1.0)*ce*rdW*i_aR*i_aR*Np; // H
  ret[2]=rdW*i_aR*i_aR; // F
  for(p=0;p<3;p++){
    ret[3+p]=cb*ce*i_aR*i_aR*vR[p]*i_aR*Np; // dG/dv
    ret[6+p]=ce*i_aR*i_aR*i_aR*( (kR*kR+3.0*cb)*vR[p]*i_aR*rdW-cb*td->cw[p][0])*Np; // dH/dv
    if(td->nid==0) ret[9+p]=i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdW-td->cw[p][0]); // dF/dv
  }
}
