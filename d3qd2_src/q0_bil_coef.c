/*
 * q0_bil_coef.c
 *
 *  Created on: Sep 11, 2019
 *      Author: ohta
 */
#include "d3qd2_elem.h"


void q0_bil_coef_4p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double vR[3],aR,i_aR,aW,tD,*vW;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  CC[8]=0.0;
  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    tD=(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
    err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_bil_coef.c, q0_bil_coef_4p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i  ]=qG*aW; // qG
    CC[i+4]=dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2]; // qH
    CC[8]+=bd->wt_44[i]*tD; // F
  }
  CC[8]*=I_FP;
}

void q0_bil_coef_9p(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3];
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_9p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_9p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<9;i++){
      zeta=bd->zt_49[i];
      eta =bd->et_49[i];
      bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

      for(l=0;l<3;l++) vR[l]-=rt[l];
      aR=vabs_d(vR);
      i_aR=1.0/aR;
      aW=vabs_d(vW);
      Np=bil_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_bil_coef.c, q0_bil_coef_9p(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_49[i]*qG*Np*aW; // qG
      CC[ep+4]+=bd->wt_49[i]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
      if(ep==0) CC[8]+=bd->wt_49[i]*(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
    }
    if(ep==0) CC[8]*=I_FP;
  }
}

void q0_bil_coef_GL(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,tmpf;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf=0.0;
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        Np=bil_Nn(ep,zeta,eta);
        err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_coef.c, q0_bil_coef_GL(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg+=bd->wli[j]*qG*Np*aW;
        tmph+=bd->wli[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        if(ep==0) tmpf+=bd->wli[j]*(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
      }
      CC[ep  ]+=bd->wli[i]*tmpg;
      CC[ep+4]+=bd->wli[i]*tmph;
      if(ep==0) CC[8]+=bd->wli[i]*tmpf;
    }
    if(ep==0) CC[8]*=I_FP;
  }
}

void q0_bil_coef_GH(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],tmpg,tmph;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,tmpf;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++) CC[l]=0.0;
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg=0.0;
      tmph=0.0;
      tmpf=0.0;
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        Np=bil_Nn(ep,zeta,eta);
        err=d3hm_qpgf_d2_v1(&qG,dqG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_coef.c, q0_bil_coef_GH(), d3hm_qpgf_d2_v1() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg+=bd->whi[j]*qG*Np*aW;
        tmph+=bd->whi[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        if(ep==0) tmpf+=bd->whi[j]*(vR[0]*vW[0]+vR[1]*vW[1]+vR[2]*vW[2])*i_aR*i_aR*i_aR;
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      if(ep==0) CC[8]+=bd->whi[i]*tmpf;
    }
    if(ep==0) CC[8]*=I_FP;
  }
}

void q0_bil_coef_DE(double complex *CC,double *rt,int s,double complex k,BOUD *bd)
{
  double bil_sF_zeta(double zeta,void *tmp); // bil_coef.c
  double complex q0_bil_sqG_zeta(double zeta,void *tmp);
  double complex q0_bil_sqH_zeta(double zeta,void *tmp);

  TDATA td;
  double err;
  int i,l;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_DE(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_DE(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  td.k=k;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  for(l=0;l<3;l++) td.rt[l]=rt[l];
  td.qpd=&(bd->qd);

  for(i=0;i<9;i++) CC[i]=0.0; //init
  // qG,qH
  for(i=0;i<4;i++){
    td.nid=i;
    CC[i  ]=deintz(q0_bil_sqG_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("q0_bil_coef.c, q0_bil_coef_DE(), q0_bil_sqG_zeta() DE integration error. nid=%d. Exit...\n",td.nid);  exit(0);  }
    CC[i+4]=deintz(q0_bil_sqH_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("q0_bil_coef.c, q0_bil_coef_DE(), q0_bil_sqH_zeta() DE integration error. nid=%d. Exit...\n",td.nid);  exit(0);  }
  }

  // F
  if(bil_check_on_plane(td.rt,td.cr,td.cw)==1){
    CC[8]=0.0; // F=0
  }
  else {
    CC[8]=I_FP*deintd(bil_sF_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("q0_bil_coef.c, q0_bil_coef_DE(), bil_sF_zeta() DE integration error. Exit...\n");  exit(1);  }
  }
}

void q0_bil_coef_bd(double complex *CC,double zeta_t,double eta_t,int s,double complex k,BOUD *bd)
{
  double bil_sF_theta(double theta,void *tmp); // bil_coef.c
  double complex q0_bil_sqG_theta(double theta,void *tmp);
  double complex q0_bil_sqH_theta(double theta,void *tmp);

  TDATA td;
  double complex tG,tH;
  double r[5],th[5],err,bdC,tD;
  int m,i,j;

  // check parameter
  if(s<0){
    printf("q0_bil_coef_bd(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_coef_bd(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
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
  for(i=0;i<9;i++) CC[i]=0.0; // init

  // qG,qH
  for(i=0;i<4;i++){
    td.nid=i;
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];

      if(SW_BDT_BIEQ==0){
        tG=0.0;
        for(m=0;m<GHN;m++) tG+=( q0_bil_sqG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +q0_bil_sqG_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tG*=0.25*(th[j+1]-th[j]);
        CC[i]+=tG;

        tH=0.0;
        for(m=0;m<GHN;m++) tH+=( q0_bil_sqH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +q0_bil_sqH_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tH*=0.25*(th[j+1]-th[j]);
        CC[i+4]+=tH;
      }
      else {
        CC[i  ]+=deintz(q0_bil_sqG_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("q0_bil_coef.c, q0_bil_coef_bd(), q0_bil_sqG_theta() DE integration error! err=%f Exit...\n",err);  exit(0);  }
        CC[i+4]+=deintz(q0_bil_sqH_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("q0_bil_coef.c, q0_bil_coef_bd(), q0_bil_sqH_theta() DE integration error! err=%f Exit...\n",err);  exit(0);  }
      }
    }
  }

  // F
  if(bil_check_plane(td.cr,td.cw)==1){ // element is plane
    CC[8]=0.0;
  }
  else {
    bdC=td.cr[0][1]*td.cw[0][2]+td.cr[1][1]*td.cw[1][2]+td.cr[2][1]*td.cw[2][2];
    for(j=0;j<4;j++){
      td.r0=r[j];      td.r1=r[j+1];
      td.th0=th[j];    td.th1=th[j+1];
      if(SW_BDT_BIEQ==0){
        tD=0.0;
        for(m=0;m<GHN;m++) tD+=( bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(th[j+1]+3.0*th[j]),&td)
                                +bil_sF_theta(0.25*(th[j+1]-th[j])*bd->xhi[m]+0.25*(3.0*th[j+1]+th[j]),&td) )*bd->whi[m];
        tD*=0.25*(th[j+1]-th[j]);
        CC[8]+=tD;
      }
      else {
        CC[8]+=deintd(bil_sF_theta,th[j],th[j+1],&td,IEPS,&err);
        if(err<0.0){ printf("q0_bil_coef.c, q0_bil_coef_bd(), bil_sF_theta() DE integration error. err=%f Exit...\n",err);  exit(0);  }
      }
    }
    CC[8]*=I_FP*bdC;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double complex q0_bil_sqG_zeta(double zeta,void *tmp)
{
  double complex q0_bil_sqG_eta(double eta,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res=deintz(q0_bil_sqG_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_bil_coef.c, q0_bil_sqG_zeta(), q0_bil_sqG_eta() DE integration error. Exit...\n");  exit(1);  }

  return res;
}

double complex q0_bil_sqG_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex qG;
  double vW[3],vR[3],aW;
  int i,err;

  bil_rw_zeta_eta(vR,vW,td->zeta,eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];
  aW=vabs_d(vW);
  err=d3hm_qpgf_d2_qG(&qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_coef.c, q0_bil_sqG_eta(), d3hm_qpgf_d2_qG() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return qG*bil_Nn(td->nid,td->zeta,eta)*aW;
}

double complex q0_bil_sqH_zeta(double zeta,void *tmp)
{
  double complex q0_bil_sqH_eta(double eta,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res=deintz(q0_bil_sqH_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_bil_coef.c, q0_bil_sqH_zeta(), q0_bil_sqH_eta() DE integration error. Exit...\n");  exit(1);  }

  return res;
}

double complex q0_bil_sqH_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex dqG[3];
  double vW[3],vR[3];
  int i,err;

  bil_rw_zeta_eta(vR,vW,td->zeta,eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_coef.c, q0_bil_dqH_eta(), d3hm_qpgf_d2_dqG() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return (dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*bil_Nn(td->nid,td->zeta,eta);
}

double complex q0_bil_sqG_theta(double theta,void *tmp)
{
  double complex q0_bil_sqG_r(double r,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("q0_bil_coef.c, q0_bil_sqG_theta() integral region error! rm=%f must be positive. Exit...\n",rm);  exit(0);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=q0_bil_sqG_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex q0_bil_sqG_r(double r,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex qG;
  double vR[3],vW[3],zeta,eta,aW;
  int i,err;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    vR[i]=r*(  td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
              +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth));
    vW[i]=td->cw[i][0]+td->cw[i][1]*zeta+td->cw[i][2]*eta;
  }
  aW=vabs_d(vW);
  err=d3hm_qpgf_d2_qG(&qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_coef.c, q0_bil_sqG_r(), d3hm_qpgf_d2_qG() error!\n");
    printf("err=%d. r=(%g, %g, %g). Exit...\n",err,vR[0],vR[1],vR[2]);
    exit(0);
  }

  return qG*bil_Nn(td->nid,zeta,eta)*aW*r;
}

double complex q0_bil_sqH_theta(double theta,void *tmp)
{
  double complex q0_bil_sqH_r(double r,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex ret;
  double rm;
  int i;

  rm=td->r0*td->r1*sin(td->th1-td->th0)/(td->r0*sin(theta-td->th0)-td->r1*sin(theta-td->th1));
  if(rm<0.0){ printf("q0_bil_coef.c, q0_bil_sqH_theta() integral region error! rm=%f must be positive. Exit...\n",rm);  exit(0);  }
  td->cth=cos(theta);
  td->sth=sin(theta);

  ret=0.0;
  for(i=0;i<td->gn;i++) ret+=q0_bil_sqH_r(0.5*rm*(td->xi[i]+1.0),td)*td->wi[i];
  ret*=0.5*rm;

  return ret;
}

double complex q0_bil_sqH_r(double r,void *tmp)
{
  TDATA *td=(TDATA*)tmp;
  double complex dqG[3];
  double vR[3],vW[3],zeta,eta;
  int i,err;

  zeta=r*td->cth+td->zeta_t;
  eta =r*td->sth+td->eta_t;
  for(i=0;i<3;i++){
    vR[i]=r*(  td->cr[i][1]*td->cth+td->cr[i][2]*td->sth
              +td->cr[i][3]*(r*td->cth*td->sth+td->zeta_t*td->sth+td->eta_t*td->cth));
    vW[i]=td->cw[i][0]+td->cw[i][1]*zeta+td->cw[i][2]*eta;
  }
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_coef.c, q0_bil_sqH_r(), d3hm_qpgf_d2_dqG() error!\n");
    printf("err=%d. r=(%g, %g, %g). Exit...\n",err,vR[0],vR[1],vR[2]);
    exit(0);
  }

  return (dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*bil_Nn(td->nid,zeta,eta)*r;
}

