/*
 * q0_lit_coef.c
 *
 *  Created on: Sep 13, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_lit_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double vR[3],*vV,aR,i_aR,aV,tD,Nn;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n",creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  for(i=0;i<9;i++) CC[i]=0.0;
  vV=bd->wen[s][0];
  aV=vabs_d(vV);

  for(i=0;i<3;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    tD=(vR[0]*vV[0]+vR[1]*vV[1]+vR[2]*vV[2])*i_aR*i_aR*i_aR;
    err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_lit_coef.c, q0_lit_coef_4p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i+0]=bd->wt_34[i]*qG;
    CC[i+4]=bd->wt_34[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2]);
    CC[8]+=bd->wt_34[i]*tD;
  }

  for(l=0;l<3;l++) vR[l]=bd->ren[s][3][l]-rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  tD=(vR[0]*vV[0]+vR[1]*vV[1]+vR[2]*vV[2])*i_aR*i_aR*i_aR;
  err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_coef_4p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
    exit(0);
  }
  for(i=0;i<3;i++){
    Nn=lit_Nn(i,bd->zt_34[3],bd->et_34[3]);
    CC[i+0]+=bd->wt_34[3]*qG*Nn;
    CC[i+4]+=bd->wt_34[3]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Nn;
  }
  CC[8]+=bd->wt_34[3]*tD;

  for(i=0;i<3;i++){
    CC[i+0]*=0.5*aV;
    CC[i+4]*=0.5;
  }
  CC[8]*=I_EP;
}

void q0_lit_coef_7p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double cr[3][4],cw[3][3],zeta,eta,vV[3],aV,vR[3],aR,i_aR,tD,Np;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_7p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_7p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  lit_copy_elem_const_rw(cr,cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,cw);
  aV=vabs_d(vV);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<3;ep++){
    for(i=0;i<7;i++){
      zeta=bd->zt_37[i];
      eta =bd->et_37[i];
      lit_r_zeta_eta(vR,zeta,eta,cr);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      tD=(vR[0]*vV[0]+vR[1]*vV[1]+vR[2]*vV[2])*i_aR*i_aR*i_aR;
      Np=lit_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_lit_coef.c, q0_lit_coef_7p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_37[i]*qG*Np;
      CC[ep+4]+=bd->wt_37[i]*(dqG[0]*vV[0]+dqG[1]*vV[1]+dqG[2]*vV[2])*Np;
      if(ep==0) CC[8]+=bd->wt_37[i]*tD;
    }
    CC[ep+0]*=0.5*aV;
    CC[ep+4]*=0.5;
    if(ep==0) CC[8]*=I_EP;
  }
}

void q0_lit_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  void q0_lit_sqGHF_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[3],tpj[3],tpi[3];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GLN;i++){
      for(l=0;l<3;l++) tpi[l]=0.0;
      a=bd->xli[i];
      for(j=0;j<GLN;j++){
        for(l=0;l<3;l++) tpj[l]=0.0;
        b=bd->xli[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];

        for(l=0;l<3;l++) tpi[l]+=tpj[l]*wt*bd->wli[j];
      }
      CC[ep  ]+=tpi[0]*bd->wli[i];
      CC[ep+4]+=tpi[1]*bd->wli[i];
      if(ep==0) CC[8]+=tpi[2]*bd->wli[i];
    }

    CC[ep  ]*=1.0/(24.0)*td.aW;
    CC[ep+4]*=1.0/(24.0);
    if(ep==0) CC[8]*=0.0625*SQ3*I_SPSQ3;
  }
}

void q0_lit_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  void q0_lit_sqGHF_zeta_eta(double complex *ret,TDATA *td);

  TDATA td;
  double complex ret[3],tpj[3],tpi[3];
  double vV[3],a,b,wt;
  int i,j,l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    for(i=0;i<GHN;i++){
      for(l=0;l<3;l++) tpi[l]=0.0;
      a=bd->xhi[i];
      for(j=0;j<GHN;j++){
        for(l=0;l<3;l++) tpj[l]=0.0;
        b=bd->xhi[j];
        wt=1.0+0.25*a+0.25*b;
        // part 1
        td.zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
        td.eta =0.125*SQ3*(-a+b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 2
        td.zeta=0.0625*(-3.0+a-5.0*b-a*b);
        td.eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];
        // part 3
        td.zeta=0.0625*(-3.0-5.0*a+b-a*b);
        td.eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
        q0_lit_sqGHF_zeta_eta(ret,&td);
        for(l=0;l<3;l++) tpj[l]+=ret[l];

        for(l=0;l<3;l++) tpi[l]+=tpj[l]*wt*bd->whi[j];
      }
      CC[ep  ]+=tpi[0]*bd->whi[i];
      CC[ep+4]+=tpi[1]*bd->whi[i];
      if(ep==0) CC[8]+=tpi[2]*bd->whi[i];
    }

    CC[ep  ]*=1.0/24.0*td.aW;
    CC[ep+4]*=1.0/24.0;
    if(ep==0) CC[8]*=0.0625*SQ3*I_SPSQ3;
  }
}

void q0_lit_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex q0_lit_sqG_alpha(double a,void *tmp);
  double complex q0_lit_sqH_alpha(double a,void *tmp);
  double lit_sF_alpha(double alpha,void *tmp); // lit_coef.c

  TDATA td;
  double vV[3],errd;
  int l,ep;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_DE(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_DE(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
  // G,H
  for(ep=0;ep<3;ep++){
    td.nid=ep;
    CC[ep+0]=1.0/24.0*td.aW*deintz(q0_lit_sqG_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0.0){ printf("q0_lit_coef.c, q0_lit_coef_DE(), q0_lit_sqG_alpha() DE integration error. err=%g Exit...\n",errd);  exit(1);  }
    CC[ep+4]=1.0/24.0*deintz(q0_lit_sqH_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0.0){ printf("q0_lit_coef.c, q0_lit_coef_DE(), q0_lit_sqH_alpha() DE integration error. err=%g Exit...\n",errd);  exit(1);  }
  }
  // F
  if(lit_check_on_plane(td.rt,td.cr,td.cw)==0){ // rt is not on element plane
    CC[8]=0.0625*SQ3*I_SPSQ3*deintd(lit_sF_alpha,-1.0,1.0,&td,IEPS,&errd);
    if(errd<0.0){ printf("q0_lit_coef.c, q0_lit_coef_DE(), lit_sF_alpha() DE integration error. err=%g Exit...\n",errd);  exit(1);  }
  }
  else CC[8]=0.0;

}

void q0_lit_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double complex q0_lit_sqG_theta(double theta,void *tmp);
  double complex q0_lit_sqH_theta(double theta,void *tmp);

  TDATA td;
  double complex tG,tH;
  double r[4],th[4],errd,vV[3],aV;
  int i,j,l,m;

  // check parameter
  if(s<0){
    printf("q0_lit_coef_bd(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_lit_coef_bd(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  td.qpd=&(bd->qd);
  td.zeta_t=zeta_t;
  td.eta_t =eta_t;
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
  lit_copy_elem_const_rw(td.cr,td.cw,s,bd);
  lit_w_zeta_eta(vV,0.0,0.0,td.cw);
  aV=vabs_d(vV);
  lit_calc_node_r_th(r,th,td.zeta_t,td.eta_t);

  for(l=0;l<9;l++) CC[l]=0.0;

  for(i=0;i<3;i++){
    td.nid=i;

    for(j=0;j<3;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      // G,H
      if(SW_BDT_BIEQ==0){
        tG=0.0;
        for(m=0;m<GHN;m++) tG+=( q0_lit_sqG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +q0_lit_sqG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tG*=0.25*(th[j+1]-th[j]);
        CC[i]+=tG;

        tH=0.0;
        for(m=0;m<GHN;m++) tH+=( q0_lit_sqH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +q0_lit_sqH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tH*=0.25*(th[j+1]-th[j]);
        CC[i+4]+=tH;
      }
      else {
        CC[i]+=deintz(q0_lit_sqG_theta,th[j],th[j+1],&td,IEPS,&errd);
        if(errd<0.0){ printf("q0_lit_coef_bd(), q0_lit_sqG_theta() DE integration error! err=%f Exit...\n",errd);  exit(0);  }

        CC[i+4]+=deintz(q0_lit_sqH_theta,th[j],th[j+1],&td,IEPS,&errd);
        if(errd<0.0){ printf("q0_lit_coef_bd(), q0_lit_sqH_theta() DE integration error! err=%f Exit...\n",errd);  exit(0);  }
      }
    }
    CC[i+0]*=2.0/(3.0*SQ3)*aV;
    CC[i+4]*=2.0/(3.0*SQ3);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_lit_sqGHF_zeta_eta(double complex *ret,TDATA *td)
{
  double complex qG,dqG[3];
  double vR[3],aR,i_aR,tD,Np;
  int l,err;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  aR=vabs_d(vR);
  i_aR=1.0/aR;
  if(td->nid==0)tD=(vR[0]*td->cw[0][0]+vR[1]*td->cw[1][0]+vR[2]*td->cw[2][0])*i_aR*i_aR*i_aR;
  else tD=0.0;
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_sqGHF_zeta_eta(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
    exit(0);
  }

  ret[0]=qG*Np;
  ret[1]=(dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*Np;
  ret[2]=tD;
}

double complex q0_lit_sqG_alpha(double a,void *tmp)
{
  double complex q0_lit_sqG_beta(double b,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=a;
  res= deintz(q0_lit_sqG_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_lit_coef.c, q0_lit_sqG_alpha(), q0_lit_sqG_beta() DE integration error. err=%g. Exit...\n",err);  exit(0);  }
  return res;
}

double complex q0_lit_sqG_beta(double b,void *tmp)
{
  double complex q0_lit_sqG_zeta_eta(TDATA *td);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0.0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=q0_lit_sqG_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=q0_lit_sqG_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=q0_lit_sqG_zeta_eta(td);

  return ret*wt;
}

double complex q0_lit_sqG_zeta_eta(TDATA *td)
{
  double complex qG;
  double vR[3],Np;
  int l,err;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_qG(&qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_sqG_zeta_eta(), d3hm_qpgf_d2_qG() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return qG*Np;
}

double complex q0_lit_sqH_alpha(double a,void *tmp)
{
  double complex q0_lit_sqH_beta(double b,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta_t=a;
  res= deintz(q0_lit_sqH_beta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_lit_coef.c, q0_lit_sqH_alpha(), q0_lit_sqH_beta() DE integration error. err=%g. Exit...\n",err);  exit(0);  }
  return res;
}

double complex q0_lit_sqH_beta(double b,void *tmp)
{
  double complex q0_lit_sqH_zeta_eta(TDATA *td);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double a,wt;

  a=td->zeta_t;
  wt=1.0+0.25*(a+b);

  ret=0.0;
  // point 1
  td->zeta=0.125*(3.0+2.0*a+2.0*b+a*b);
  td->eta =0.125*SQ3*(-a+b);
  ret+=q0_lit_sqH_zeta_eta(td);
  // part 2
  td->zeta=0.0625*(-3.0+a-5.0*b-a*b);
  td->eta =0.0625*SQ3*(3.0+3.0*a+b+a*b);
  ret+=q0_lit_sqH_zeta_eta(td);
  // part 3
  td->zeta=0.0625*(-3.0-5.0*a+b-a*b);
  td->eta =0.0625*SQ3*(-3.0-a-3.0*b-a*b);
  ret+=q0_lit_sqH_zeta_eta(td);

  return ret*wt;
}

double complex q0_lit_sqH_zeta_eta(TDATA *td)
{
  double complex dqG[3];
  double vR[3],Np;
  int l,err;

  lit_r_zeta_eta(vR,td->zeta,td->eta,td->cr);

  for(l=0;l<3;l++) vR[l]-=td->rt[l];
  Np=lit_Nn(td->nid,td->zeta,td->eta);
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_sqH_zeta_eta(), d3hm_qpgf_d2_dqG() error. err=%d. Exit...\n",err);
    printf("r=(% 15.14e, % 15.14e, % 15.14e)\n",vR[0],vR[1],vR[2]);
    exit(0);
  }

  return (dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*Np;
}

double complex q0_lit_sqG_theta(double theta,void *tmp)
{
  double complex q0_lit_sqG_r(double r,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("q0_lit_coef.c, q0_lit_sqG_theta() integral region error! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=q0_lit_sqG_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex q0_lit_sqG_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex qG;
  double vR[3],zeta,eta;
  int i,err;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) vR[i]=r*(td->cr[i][1]*td->cth+td->cr[i][2]*td->sth);
  err=d3hm_qpgf_d2_qG(&qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_sqG_r(), d3hm_qpgf_d2_qG() error!\n");
    printf("err=%d. r=(%g, %g, %g). Exit...\n",err,vR[0],vR[1],vR[2]);
    exit(0);
  }

  return qG*lit_Nn(td->nid,zeta,eta)*r;
}

double complex q0_lit_sqH_theta(double theta,void *tmp)
{
  double complex q0_lit_sqH_r(double r,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("q0_lit_coef.c, q0_lit_sqH_theta() integral region error! rm=%f must be positive. Exit...\n",rm);  exit(1);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=q0_lit_sqH_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex q0_lit_sqH_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex dqG[3];
  double vR[3],zeta,eta;
  int i,err;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++) vR[i]=r*(td->cr[i][1]*td->cth+td->cr[i][2]*td->sth);
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_lit_coef.c, q0_lit_sqH_r(), d3hm_qpgf_d2_dqG() error!\n");
    printf("err=%d. r=(%g, %g, %g). Exit...\n",err,vR[0],vR[1],vR[2]);
    exit(0);
  }

  return (dqG[0]*td->cw[0][0]+dqG[1]*td->cw[1][0]+dqG[2]*td->cw[2][0])*lit_Nn(td->nid,zeta,eta)*r;
}
