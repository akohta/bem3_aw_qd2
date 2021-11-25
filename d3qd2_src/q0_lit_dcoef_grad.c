/*
 * q0_lit_dcoef_grad.c
 *
 *  Created on: Sep 27, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_lit_dcoef_grad_4p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,aV,*vV,rdV,Nn,vp[3][3];
  int i,l,err,p;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_grad_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_grad_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  for(i=0;i<9;i++) CC [i]=0.0;
  for(p=0;p<3;p++){
    for(i=0;i<9;i++) dC[p][i]=0.0;
    for(i=0;i<3;i++) vp[p][i]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  vV=bd->wen[s][0];
  aV=vabs_d(vV);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    rdV=vdot_d(vR,vV)*i_aR;
    err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_lit_dcoef_grad.c, q0_lit_dcoef_grad_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i+0]=bd->wt_34[i]*qG;
    CC[i+4]=bd->wt_34[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2]);
    CC[8]+=bd->wt_34[i]*rdV*i_aR*i_aR;

    for(p=0;p<3;p++){
      dC[p][i+0]=bd->wt_34[i]*(dqG[0]*vp[p][0]+dqG[1]*vp[p][1]+dqG[2]*vp[p][2]);
      dC[p][i+4]=bd->wt_34[i]*( d2qG[0]*vV[0]*vp[p][0]+d2qG[1]*vV[1]*vp[p][1]+d2qG[2]*vV[2]*vp[p][2]
                               +d2qG[3]*(vV[0]*vp[p][1]+vV[1]*vp[p][0])
                               +d2qG[4]*(vV[1]*vp[p][2]+vV[2]*vp[p][1])
                               +d2qG[5]*(vV[2]*vp[p][0]+vV[0]*vp[p][2]) );
      dC[p][8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdV-vV[p]);
    }
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
  if(err<0){
    printf("q0_lit_dcoef_grad.c, q0_lit_dcoef_grad_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*qG*Nn;
    CC[i+4]+=bd->wt_34[3]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;

    for(p=0;p<3;p++){
      dC[p][i+0]+=bd->wt_34[3]*(dqG[0]*vp[p][0]+dqG[1]*vp[p][1]+dqG[2]*vp[p][2])*Nn;
      dC[p][i+4]+=bd->wt_34[3]*( d2qG[0]*vV[0]*vp[p][0]+d2qG[1]*vV[1]*vp[p][1]+d2qG[2]*vV[2]*vp[p][2]
                                +d2qG[3]*(vV[0]*vp[p][1]+vV[1]*vp[p][0])
                                +d2qG[4]*(vV[1]*vp[p][2]+vV[2]*vp[p][1])
                                +d2qG[5]*(vV[2]*vp[p][0]+vV[0]*vp[p][2]))*Nn;
    }
  }
  CC [8]+=bd->wt_34[3]*rdV*i_aR*i_aR;
  for(p=0;p<3;p++) dC[p][8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdV-vV[p]);

  for(i=0;i<3;i++){
    CC[i+0]*=0.5*aV;
    CC[i+4]*=0.5;
    for(p=0;p<3;p++){
      dC[p][i+0]*=-0.5*aV;
      dC[p][i+4]*=-0.5;
    }
  }
  CC[8]*=I_EP;
  for(p=0;p<3;p++) dC[p][8]*=I_EP;
}

void q0_lit_dcoef_grad_7p(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double cr[3][4],cw[3][3],zeta,eta,*vV,aV,vR[3],aR,i_aR,rdV,Nn,vp[3][3];
  int i,l,ep,err,p;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_grad_7p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_grad_7p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  vV=bd->wen[s][0];
  aV=vabs_d(vV);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(p=0;p<3;p++){
    for(l=0;l<9;l++) dC[p][l]=0.0;
    for(l=0;l<3;l++) vp[p][l]=0.0;
  }
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      rdV=vdot_d(vR,vV)*i_aR;
      Nn=lit_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_lit_dcoef_grad.c, q0_lit_dcoef_grad_7p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_37[i]*qG*Nn;
      CC[ep+4]+=bd->wt_37[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;
      for(p=0;p<3;p++){
        dC[p][ep+0]+=bd->wt_37[i]*(dqG[0]*vp[p][0]+dqG[1]*vp[p][1]+dqG[2]*vp[p][2])*Nn;
        dC[p][ep+4]+=bd->wt_37[i]*( d2qG[0]*vV[0]*vp[p][0]+d2qG[1]*vV[1]*vp[p][1]+d2qG[2]*vV[2]*vp[p][2]
                                   +d2qG[3]*(vV[0]*vp[p][1]+vV[1]*vp[p][0])
                                   +d2qG[4]*(vV[1]*vp[p][2]+vV[2]*vp[p][1])
                                   +d2qG[5]*(vV[2]*vp[p][0]+vV[0]*vp[p][2]))*Nn;
      }
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdV*i_aR*i_aR;
        for(p=0;p<3;p++) dC[p][8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdV-vV[p]);
      }
    }
    CC[ep+0]*=0.5*aV;
    CC[ep+4]*=0.5;
    for(p=0;p<3;p++){
      dC[p][ep+0]*=-0.5*aV;
      dC[p][ep+4]*=-0.5;
    }
    if(ep==0){
      CC [8]*=I_EP;
      for(p=0;p<3;p++) dC[p][8]*=I_EP;
    }
  }
}

void q0_lit_dcoef_grad_GL(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_grad_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[12],tpj[12],tpi[12];
  double vV[3],a,b,wt;
  int i,j,l,ep,p;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_grad_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_grad_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.qpd=&(bd->qd);
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,td.cw);
  td.aW=vabs_d(vV);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC[l]=0.0;
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
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];

        for(l=0;l<12;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC[ep  ]+=tpi[0]*bd->wli[i];
      CC[ep+4]+=tpi[1]*bd->wli[i];
      for(p=0;p<3;p++){
        dC[p][ep  ]+=tpi[3+p]*bd->wli[i];
        dC[p][ep+4]+=tpi[6+p]*bd->wli[i];
      }
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        for(p=0;p<3;p++) dC[p][8]+=tpi[9+p]*bd->wli[i];
      }
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-1.0/24.0*td.aW;
      dC[p][ep+4]*=-1.0/24.0;
    }
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      for(p=0;p<3;p++) dC[p][8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}

void q0_lit_dcoef_grad_GH(double complex *CC,double complex dC[][9],double *rt,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_grad_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[12],tpj[12],tpi[12];
  double vV[3],a,b,wt;
  int i,j,l,ep,p;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_grad_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_grad_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.qpd=&(bd->qd);
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,td.cw);
  td.aW=vabs_d(vV);
  for(l=0;l<3;l++) td.rt[l]=rt[l];

  for(l=0;l<9;l++) CC[l]=0.0;
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
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
        for(l=0;l<12;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_grad_zeta_eta(ret,&td);
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

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    for(p=0;p<3;p++){
      dC[p][ep  ]*=-1.0/24.0*td.aW;
      dC[p][ep+4]*=-1.0/24.0;
    }
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      for(p=0;p<3;p++) dC[p][8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}


//////////////////////////////////////////////////////////////////////////////////////////
void q0_lit_sdqGHF_grad_zeta_eta(double complex *ret,TDATA *td)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],vV[3],aR,i_aR,Np,rdV,vp[3][3];
  int l,err,p;

  for(p=0;p<3;p++)
    for(l=0;l<3;l++) vp[p][l]=0.0;
  vp[0][0]=1.0;
  vp[1][1]=1.0;
  vp[2][2]=1.0;

  vV[0]=td->cw[0][0];  vV[1]=td->cw[1][0];  vV[2]=td->cw[2][0];
  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);
  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_dcoef_grad.c, q0_lit_sdqGHF_grad_zeta_eta(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  ret[0]=qG*Np;
  ret[1]=(dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*Np;

  for(p=0;p<3;p++){
    ret[3+p]=(dqG[0]*vp[p][0]+dqG[1]*vp[p][1]+dqG[2]*vp[p][2])*Np;
    ret[6+p]=( d2qG[0]*vV[0]*vp[p][0]+d2qG[1]*vV[1]*vp[p][1]+d2qG[2]*vV[2]*vp[p][2]
              +d2qG[3]*(vV[0]*vp[p][1]+vV[1]*vp[p][0])
              +d2qG[4]*(vV[1]*vp[p][2]+vV[2]*vp[p][1])
              +d2qG[5]*(vV[2]*vp[p][0]+vV[0]*vp[p][2]))*Np;
  }

  if(td->nid==0){
    ret[2]=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
    for(p=0;p<3;p++) ret[9+p]=i_aR*i_aR*i_aR*( 3.0*vR[p]*i_aR*rdV-vV[p]);
  }
  else {
    ret[2]=0.0;
    for(p=0;p<3;p++) ret[9+p]=0.0;
  }
}
