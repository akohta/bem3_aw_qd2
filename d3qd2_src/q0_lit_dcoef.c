/*
 * q0_lit_dcoef.c
 *
 *  Created on: Sep 17, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_lit_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,aV,*vV,rdV,rdv,Vdv,Nn;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  for(i=0;i<9;i++){
    CC [i]=0.0;
    dCC[i]=0.0;
  }
  vV=bd->wen[s][0];
  aV=vabs_d(vV);
  Vdv=vdot_d(vV,v);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    rdV=vdot_d(vR,vV)*i_aR;
    rdv=vdot_d(vR,v)*i_aR;
    err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_lit_dcoef.c, q0_lit_dcoef_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i+0]=bd->wt_34[i]*qG;
    CC[i+4]=bd->wt_34[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2]);
    CC[8]+=bd->wt_34[i]*rdV*i_aR*i_aR;

    dCC[i+0]=bd->wt_34[i]*(dqG[0]*v[0]+dqG[1]*v[1]+dqG[2]*v[2]);
    dCC[i+4]=bd->wt_34[i]*( d2qG[0]*vV[0]*v[0]+d2qG[1]*vV[1]*v[1]+d2qG[2]*vV[2]*v[2]
                           +d2qG[3]*(vV[0]*v[1]+vV[1]*v[0])
                           +d2qG[4]*(vV[1]*v[2]+vV[2]*v[1])
                           +d2qG[5]*(vV[2]*v[0]+vV[0]*v[2]) );
    dCC[8]+=bd->wt_34[i]*i_aR*i_aR*i_aR*( 3.0*rdv*rdV-Vdv );
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  rdv=vdot_d(vR,v)*i_aR;
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
  if(err<0){
    printf("q0_lit_dcoef.c, q0_lit_dcoef_4p(), d3hm_qpgf_d2v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*qG*Nn;
    CC[i+4]+=bd->wt_34[3]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;

    dCC[i+0]+=bd->wt_34[3]*(dqG[0]*v[0]+dqG[1]*v[1]+dqG[2]*v[2])*Nn;
    dCC[i+4]+=bd->wt_34[3]*( d2qG[0]*vV[0]*v[0]+d2qG[1]*vV[1]*v[1]+d2qG[2]*vV[2]*v[2]
                            +d2qG[3]*(vV[0]*v[1]+vV[1]*v[0])
                            +d2qG[4]*(vV[1]*v[2]+vV[2]*v[1])
                            +d2qG[5]*(vV[2]*v[0]+vV[0]*v[2]))*Nn;
  }
  CC [8]+=bd->wt_34[3]*rdV*i_aR*i_aR;
  dCC[8]+=bd->wt_34[3]*i_aR*i_aR*i_aR*( 3.0*rdv*rdV-Vdv );

  for(i=0;i<3;i++){
    CC[i+0]*=0.5*aV;
    CC[i+4]*=0.5;
    dCC[i+0]*=-0.5*aV;
    dCC[i+4]*=-0.5;
  }
  CC[8]*=I_EP;
  dCC[8]*=I_EP;
}

void q0_lit_dcoef_7p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double cr[3][4],cw[3][3],zeta,eta,*vV,aV,vR[3],aR,i_aR,rdV,rdv,Vdv,Nn;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_7p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_7p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  vV=bd->wen[s][0];
  aV=vabs_d(vV);
  Vdv=vdot_d(vV,v);

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dCC[l]=0.0;
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
      rdv=vdot_d(vR,v)*i_aR;
      Nn=lit_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_lit_dcoef.c, q0_lit_dcoef_7p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC [ep  ]+=bd->wt_37[i]*qG*Nn;
      CC [ep+4]+=bd->wt_37[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;
      dCC[ep+0]+=bd->wt_37[i]*(dqG[0]* v[0]+dqG[1]* v[1]+dqG[2]* v[2])*Nn;
      dCC[ep+4]+=bd->wt_37[i]*( d2qG[0]*vV[0]*v[0]+d2qG[1]*vV[1]*v[1]+d2qG[2]*vV[2]*v[2]
                               +d2qG[3]*(vV[0]*v[1]+vV[1]*v[0])
                               +d2qG[4]*(vV[1]*v[2]+vV[2]*v[1])
                               +d2qG[5]*(vV[2]*v[0]+vV[0]*v[2]))*Nn;
      if(ep==0){
        CC [8]+=bd->wt_37[i]*rdV*i_aR*i_aR;
        dCC[8]+=bd->wt_37[i]*i_aR*i_aR*i_aR*( 3.0*rdv*rdV-Vdv );
      }
    }
    CC[ep+0]*=0.5*aV;
    CC[ep+4]*=0.5;
    dCC[ep+0]*=-0.5*aV;
    dCC[ep+4]*=-0.5;
    if(ep==0){
      CC [8]*=I_EP;
      dCC[8]*=I_EP;
    }
  }
}

void q0_lit_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[6],tpj[6],tpi[6];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
    td.vt[l]=v[l];
  }

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dCC[l]=0.0;
  }

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<6;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<6;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];

        for(l=0;l<6;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC[ep  ]+=tpi[0]*bd->wli[i];
      CC[ep+4]+=tpi[1]*bd->wli[i];
      dCC[ep  ]+=tpi[3]*bd->wli[i];
      dCC[ep+4]+=tpi[4]*bd->wli[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->wli[i];
        dCC[8]+=tpi[5]*bd->wli[i];
      }
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    dCC[ep  ]*=-1.0/24.0*td.aW;
    dCC[ep+4]*=-1.0/24.0;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dCC[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}

void q0_lit_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  void q0_lit_sdqGHF_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[6],tpj[6],tpi[6];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
    td.vt[l]=v[l];
  }

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dCC[l]=0.0;
  }

  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<6;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<6;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sdqGHF_zeta_eta(ret,&td);
        for(l=0;l<6;l++) tpj[l]+=ret[l];

        for(l=0;l<6;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC [ep  ]+=tpi[0]*bd->whi[i];
      CC [ep+4]+=tpi[1]*bd->whi[i];
      dCC[ep  ]+=tpi[3]*bd->whi[i];
      dCC[ep+4]+=tpi[4]*bd->whi[i];
      if(ep==0){
        CC [8]+=tpi[2]*bd->whi[i];
        dCC[8]+=tpi[5]*bd->whi[i];
      }
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    dCC[ep  ]*=-1.0/24.0*td.aW;
    dCC[ep+4]*=-1.0/24.0;
    if(ep==0){
      CC [8]*=0.0625*SQ3*I_SPSQ3;
      dCC[8]*=0.0625*SQ3*I_SPSQ3;
    }
  }
}

void q0_lit_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  void q0_lit_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);

  q0_lit_coef_DE(CC,rt,s,k,bd);
  q0_lit_dcoef_DE0(dCC,rt,v,s,k,bd);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_lit_sdqGHF_zeta_eta(double complex *ret,TDATA *td)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],vV[3],aR,i_aR,Np,rdV,rdv,Vdv;
  int l,err;

  vV[0]=td->cw[0][0];  vV[1]=td->cw[1][0];  vV[2]=td->cw[2][0];
  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);
  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  rdV=vdot_d(vR,vV)*i_aR;
  rdv=vdot_d(vR,td->vt)*i_aR;
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_dcoef.c, q0_lit_sdqGHF_zeta_eta(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  ret[0]=qG*Np;
  ret[1]=(dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*Np;

  ret[3]=(dqG[0]*td->vt[0]+dqG[1]*td->vt[1]+dqG[2]*td->vt[2])*Np;
  ret[4]=( d2qG[0]*vV[0]*td->vt[0]+d2qG[1]*vV[1]*td->vt[1]+d2qG[2]*vV[2]*td->vt[2]
          +d2qG[3]*(vV[0]*td->vt[1]+vV[1]*td->vt[0])
          +d2qG[4]*(vV[1]*td->vt[2]+vV[2]*td->vt[1])
          +d2qG[5]*(vV[2]*td->vt[0]+vV[0]*td->vt[2]) )*Np;

  if(td->nid==0){
    ret[2]=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
    Vdv=vdot_d(vV,td->vt);
    ret[5]=i_aR*i_aR*i_aR*( 3.0*rdv*rdV-Vdv );
  }
  else {
    ret[2]=0.0;
    ret[5]=0.0;
  }
}

void q0_lit_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex q0_lit_sdqG_alpha(double a,void *tmp);
  double complex q0_lit_sdqH_alpha(double a,void *tmp);
  double lit_sdF_alpha(double alpha,void *tmp); // lit_dcoef.c

  TDATA td;
  double errd,vV[3];
  int i,l;

  // check parameter
  if(s<0){
    printf("q0_lit_dcoef.c, q0_lit_dcoef_DE0() signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_dcoef.c, q0_lit_dcoef_DE0() wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
    td.vt[l]=v[l];
  }

  for(i=0;i<9;i++) dCC[i]=0.0;

  // dqG/dv, dqH/dv
  for(i=0;i<3;i++){
    td.nid=i;
    dCC[i  ]=-1.0/24.0*td.aW*deintz(q0_lit_sdqG_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0){ printf("q0_lit_dcoef.c, q0_lit_dcoef_DE0(), q0_lit_sdqG_alpha() DE integration error! nid=%d. err=%g. Exit...\n",td.nid,errd);  exit(0);  }
    dCC[i+4]=-1.0/24.0*deintz(q0_lit_sdqH_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0){ printf("q0_lit_dcoef.c, q0_lit_dcoef_DE0(), q0_lit_sdqH_alpha() DE integration error! nid=%d. err=%g. Exit...\n",td.nid,errd);  exit(0);  }
  }
  // dF/dv
  if(lit_check_on_plane(td.rt,td.cr,td.cw)==0){ // rt is not on element plane
    dCC[8]=0.0625*SQ3*I_SPSQ3*deintd(lit_sdF_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0.0){ printf("q0_lit_dcoef.c, q0_lit_dcoef_DE0(), lit_sdF_alpha() DE integration error! err=%g Exit...\n",errd);  exit(0);  }
  }
  else dCC[8]=0.0;
}

double complex q0_lit_sdqG_alpha(double a,void *tmp)
{
  double complex q0_lit_sdqG_beta(double b,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double errd;

  td->zeta_t=a;
  res= deintz(q0_lit_sdqG_beta,-1.0,1.0,td,IEPS,&errd);
  if(errd<0.0){ printf("q0_lit_dcoef.c, q0_lit_sdqG_alpha(), q0_lit_sdqG_beta() DE integration error. err=%g. Exit...\n",errd);  exit(0);  }
  return res;
}

double complex q0_lit_sdqG_beta(double b,void *tmp)
{
  double complex q0_lit_sdqG_zeta_eta(TDATA *td);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=q0_lit_sdqG_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=q0_lit_sdqG_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=q0_lit_sdqG_zeta_eta(td);

  return ret*wt;
}

double complex q0_lit_sdqG_zeta_eta(TDATA *td)
{
  double complex dqG[3];
  double vW[3],vR[3];
  int i,err;

  lit_rw_zeta_eta(vR,vW,td->zeta,td->eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_dcoef.c, q0_lit_sdqG_zeta_eta(), d3hm_qpgf_d2_dqG() error. err=%d. Exit...\n",err);
    printf("r=(%15.14e, %15.14e, %15.14e)\n",vR[0],vR[1],vR[2]);
    exit(0);
  }

  return (dqG[0]*td->vt[0]+dqG[1]*td->vt[1]+dqG[2]*td->vt[2])*lit_Nn(td->nid,td->zeta,td->eta);
}

double complex q0_lit_sdqH_alpha(double a,void *tmp)
{
  double complex q0_lit_sdqH_beta(double b,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double errd;

  td->zeta_t=a;
  res= deintz(q0_lit_sdqH_beta,-1.0,1.0,td,IEPS,&errd);
  if(errd<0.0){ printf("q0_lit_dcoef.c, q0_lit_sdqH_alpha(), q0_lit_sdqH_beta() DE integration error. err=%g. Exit...\n",errd);  exit(0);  }
  return res;
}

double complex q0_lit_sdqH_beta(double b,void *tmp)
{
  double complex q0_lit_sdqH_zeta_eta(TDATA *td);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=q0_lit_sdqH_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=q0_lit_sdqH_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=q0_lit_sdqH_zeta_eta(td);

  return ret*wt;
}

double complex q0_lit_sdqH_zeta_eta(TDATA *td)
{
  double complex qG,dqG[3],d2qG[6];
  double vV[3],vR[3];
  int i,err;

  lit_rw_zeta_eta(vR,vV,td->zeta,td->eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];

  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_dcoef.c, q0_lit_sdqH_zeta_eta(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return ( d2qG[0]*vV[0]*td->vt[0]
          +d2qG[1]*vV[1]*td->vt[1]
          +d2qG[2]*vV[2]*td->vt[2]
          +d2qG[3]*(vV[0]*td->vt[1]+vV[1]*td->vt[0])
           +d2qG[4]*(vV[1]*td->vt[2]+vV[2]*td->vt[1])
          +d2qG[5]*(vV[2]*td->vt[0]+vV[0]*td->vt[2]) )*lit_Nn(td->nid,td->zeta,td->eta);
}
