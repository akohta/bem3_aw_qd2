/*
 * q0_bil_dcoef.c
 *
 *  Created on: Sep 12, 2019
 *      Author: ohta
 */

#include "d3qd2_elem.h"

void q0_bil_dcoef_4p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,aW,*vW,rdW,rdv,Wdv;
  int i,l,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_4p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_4p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  CC[8]=0.0;
  dCC[8]=0.0;
  for(i=0;i<4;i++){
    for(l=0;l<3;l++) vR[l]=bd->ren[s][i][l]-rt[l];
    vW=bd->wen[s][i];

    aR=vabs_d(vR);
    i_aR=1.0/aR;
    aW=vabs_d(vW);
    rdW=vdot_d(vR,vW)*i_aR;
    rdv=vdot_d(vR,v)*i_aR;
    Wdv=vdot_d(vW,v);
    err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
    if(err<0){
      printf("q0_bil_dcoef.c, q0_bil_dcoef_4p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
      exit(0);
    }

    CC[i  ]=qG*aW; // qG
    CC[i+4]=dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2]; // qH
    CC[8]+=bd->wt_44[i]*rdW*i_aR*i_aR; // F

    dCC[i  ]=-(dqG[0]*v[0]+dqG[1]*v[1]+dqG[2]*v[2])*aW; // dqG/dv
    dCC[i+4]=-(d2qG[0]*vW[0]*v[0]+d2qG[1]*vW[1]*v[1]+d2qG[2]*vW[2]*v[2]
               +d2qG[3]*(vW[0]*v[1]+vW[1]*v[0])+d2qG[4]*(vW[1]*v[2]+vW[2]*v[1])+d2qG[5]*(vW[2]*v[0]+vW[0]*v[2]));
    dCC[8]+=bd->wt_44[i]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv); // dF/dv
  }
  CC[8]*=I_FP;
  dCC[8]*=I_FP;
}

void q0_bil_dcoef_9p(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6];
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],Np,rdv,rdW,Wdv;
  int i,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_9p(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_9p(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n",creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){
    CC[l]=0.0;
    dCC[l]=0.0;
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
      rdv=vdot_d(vR,v)*i_aR;
      Wdv=vdot_d(vW,v);
      Np=bil_Nn(ep,zeta,eta);
      err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
      if(err<0){
        printf("q0_bil_dcoef.c, q0_bil_dcoef_9p(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
        exit(0);
      }

      CC[ep  ]+=bd->wt_49[i]*qG*Np*aW; // qG
      CC[ep+4]+=bd->wt_49[i]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np; // qH

      dCC[ep  ]+=bd->wt_49[i]*(dqG[0]*v[0]+dqG[1]*v[1]+dqG[2]*v[2])*aW*Np; // dqG/dv
      dCC[ep+4]+=bd->wt_49[i]*(d2qG[0]*vW[0]*v[0]+d2qG[1]*vW[1]*v[1]+d2qG[2]*vW[2]*v[2]
                                +d2qG[3]*(vW[0]*v[1]+vW[1]*v[0])+d2qG[4]*(vW[1]*v[2]+vW[2]*v[1])+d2qG[5]*(vW[2]*v[0]+vW[0]*v[2]))*Np; // dqH/dv
      if(ep==0){
        CC[8]+=bd->wt_49[i]*rdW*i_aR*i_aR;
        dCC[8]+=bd->wt_49[i]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
      }
    }
    dCC[ep  ]*=-1.0;
    dCC[ep+4]*=-1.0;
    if(ep==0){
      CC[8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_GL(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6],tmpg,tmph,tmpdg,tmpdh;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv,Wdv,Np,tmpf,tmpdf;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_GL(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_GL(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
        ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){ // init
    CC[l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GLN;i++){
      tmpg =0.0;
      tmpdg=0.0;
      tmph =0.0;
      tmpdh=0.0;
      tmpf =0.0;
      tmpdf=0.0;
      for(j=0;j<GLN;j++){
        zeta=bd->xli[i];
        eta =bd->xli[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        rdv=vdot_d(vR,v)*i_aR;
        Wdv=vdot_d(vW,v);
        Np=bil_Nn(ep,zeta,eta);
        err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef.c, q0_bil_dcoef_GL(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->wli[j]*qG*Np*aW;
        tmpdg+=bd->wli[j]*(dqG[0]* v[0]+dqG[1]* v[1]+dqG[2]* v[2])*aW*Np;
        tmph +=bd->wli[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        tmpdh+=bd->wli[j]*(d2qG[0]*vW[0]*v[0]+d2qG[1]*vW[1]*v[1]+d2qG[2]*vW[2]*v[2]
                          +d2qG[3]*(vW[0]*v[1]+vW[1]*v[0])+d2qG[4]*(vW[1]*v[2]+vW[2]*v[1])+d2qG[5]*(vW[2]*v[0]+vW[0]*v[2]))*Np;
        if(ep==0){
          tmpf +=bd->wli[j]*rdW*i_aR*i_aR;
          tmpdf+=bd->wli[j]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
        }
      }
      CC[ep  ]+=bd->wli[i]*tmpg;
      CC[ep+4]+=bd->wli[i]*tmph;
      dCC[ep  ]+=bd->wli[i]*tmpdg;
      dCC[ep+4]+=bd->wli[i]*tmpdh;
      if(ep==0){
        CC[8] +=bd->wli[i]*tmpf;
        dCC[8]+=bd->wli[i]*tmpdf;
      }
    }
    dCC[ep  ]*=-1.0;
    dCC[ep+4]*=-1.0;
    if(ep==0){
      CC[8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_GH(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double complex qG,dqG[3],d2qG[6],tmpg,tmph,tmpdg,tmpdh;
  double vR[3],aR,i_aR,vW[3],aW,zeta,eta,cr[3][4],cw[3][3],rdW,rdv,Wdv,Np,tmpf,tmpdf;
  int i,j,l,ep,err;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_GH(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0 && creal(k)==bd->qd.k){
    printf("q0_bil_dcoef_GH(), wave number error. k must be same as open region wave number. k=%g %+gI. qd.k=%g. Exit...\n"
          ,creal(k),cimag(k),bd->qd.k);
    exit(0);
  }

  bil_copy_elem_const_rw(cr,cw,s,bd);

  for(l=0;l<9;l++){ // init
    CC[l]=0.0;
    dCC[l]=0.0;
  }
  for(ep=0;ep<4;ep++){
    for(i=0;i<GHN;i++){
      tmpg =0.0;
      tmpdg=0.0;
      tmph =0.0;
      tmpdh=0.0;
      tmpf =0.0;
      tmpdf=0.0;
      for(j=0;j<GHN;j++){
        zeta=bd->xhi[i];
        eta =bd->xhi[j];
        bil_rw_zeta_eta(vR,vW,zeta,eta,cr,cw);

        for(l=0;l<3;l++) vR[l]-=rt[l];
        aR=sqrt(vR[0]*vR[0]+vR[1]*vR[1]+vR[2]*vR[2]);
        i_aR=1.0/aR;
        aW=sqrt(vW[0]*vW[0]+vW[1]*vW[1]+vW[2]*vW[2]);
        rdW=vdot_d(vR,vW)*i_aR;
        rdv=vdot_d(vR,v)*i_aR;
        Wdv=vdot_d(vW,v);
        Np=bil_Nn(ep,zeta,eta);
        err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,&(bd->qd));
        if(err<0){
          printf("q0_bil_dcoef.c, q0_bil_dcoef_GH(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
          exit(0);
        }

        tmpg +=bd->whi[j]*qG*Np*aW;
        tmpdg+=bd->whi[j]*(dqG[0]* v[0]+dqG[1]* v[1]+dqG[2]* v[2])*aW*Np;
        tmph +=bd->whi[j]*(dqG[0]*vW[0]+dqG[1]*vW[1]+dqG[2]*vW[2])*Np;
        tmpdh+=bd->whi[j]*(d2qG[0]*vW[0]*v[0]+d2qG[1]*vW[1]*v[1]+d2qG[2]*vW[2]*v[2]
                          +d2qG[3]*(vW[0]*v[1]+vW[1]*v[0])+d2qG[4]*(vW[1]*v[2]+vW[2]*v[1])+d2qG[5]*(vW[2]*v[0]+vW[0]*v[2]))*Np;
        if(ep==0){
          tmpf +=bd->whi[j]*rdW*i_aR*i_aR;
          tmpdf+=bd->whi[j]*i_aR*i_aR*i_aR*(3.0*rdv*rdW-Wdv);
        }
      }
      CC[ep  ]+=bd->whi[i]*tmpg;
      CC[ep+4]+=bd->whi[i]*tmph;
      dCC[ep  ]+=bd->whi[i]*tmpdg;
      dCC[ep+4]+=bd->whi[i]*tmpdh;
      if(ep==0){
        CC[8] +=bd->whi[i]*tmpf;
        dCC[8]+=bd->whi[i]*tmpdf;
      }
    }
    dCC[ep  ]*=-1.0;
    dCC[ep+4]*=-1.0;
    if(ep==0){
      CC[8]*=I_FP;
      dCC[8]*=I_FP;
    }
  }
}

void q0_bil_dcoef_DE(double complex *CC,double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  void q0_bil_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd);

  q0_bil_coef_DE(CC,rt,s,k,bd);
  q0_bil_dcoef_DE0(dCC,rt,v,s,k,bd);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void q0_bil_dcoef_DE0(double complex *dCC,double *rt,double *v,int s,double complex k,BOUD *bd)
{
  double bil_sdF_zeta(double zeta,void *tmp); // bil_dcoef.c
  double complex q0_bil_sdqG_zeta(double zeta,void *tmp);
  double complex q0_bil_sdqH_zeta(double zeta,void *tmp);

  TDATA td;
  double err;
  int i,l;

  // check parameter
  if(s<0){
    printf("q0_bil_dcoef_DE0(), signed element id error. s must be positive number. s=%d. Exit...\n",s);
    exit(0);
  }
  if(cimag(k)!=0.0){
    printf("q0_bil_dcoef_DE0(), wave number error. k must be a real number. k=%g %+g. Exit...\n",creal(k),cimag(k));
    exit(0);
  }

  td.k=k;
  bil_copy_elem_const_rw(td.cr,td.cw,s,bd);
  for(l=0;l<3;l++){
    td.rt[l]=rt[l];
    td.vt[l]=v[l];
  }
  td.qpd=&(bd->qd);

  for(i=0;i<9;i++) dCC[i]=0.0; //init
  // dqG/dv, dqH/dv
  for(i=0;i<4;i++){
    td.nid=i;
    dCC[i  ]=-deintz(q0_bil_sdqG_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("q0_bil_dcoef.c, q0_bil_dcoef_DE0(), q0_bil_sdqG_zeta() DE integration error! nid=%d. err=%g. Exit...\n",td.nid,err);  exit(0);  }
    dCC[i+4]=-deintz(q0_bil_sdqH_zeta,-1.0,1.0,&td,IEPS,&err);
    if(err<0){ printf("q0_bil_dcoef.c, q0_bil_dcoef_DE0(), q0_bil_sdqH_zeta() DE integration error! nid=%d. err=%g. Exit...\n",td.nid,err);  exit(0);  }
  }

  // dF/dv
  dCC[8]=I_FP*deintd(bil_sdF_zeta,-1.0,1.0,&td,IEPS,&err);
  if(err<0){ printf("q0_bil_dcoef.c, q0_bil_dcoef_DE0(), bil_sdF_zeta() DE integration error! err=%g. Exit...\n",err);  exit(0);  }
}

double complex q0_bil_sdqG_zeta(double zeta,void *tmp)
{
  double complex q0_bil_sdqG_eta(double eta,void *tmp);

  TDATA *td=(TDATA*)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res= deintz(q0_bil_sdqG_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_bil_dcoef.c, q0_bil_sdqG_zeta(), q0_bil_sdqG_eta() DE integration error! err=%g. Exit...\n",err);  exit(0);  }
  return res;
}

double complex q0_bil_sdqG_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex dqG[3];
  double vW[3],vR[3],aW;
  int i,err;

  bil_rw_zeta_eta(vR,vW,td->zeta,eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];
  aW=vabs_d(vW);
  err=d3hm_qpgf_d2_dqG(dqG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_dcoef.c, q0_bil_sdqG_eta(), d3hm_qpgf_d2_dqG() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return (dqG[0]*td->vt[0]+dqG[1]*td->vt[1]+dqG[2]*td->vt[2])*bil_Nn(td->nid,td->zeta,eta)*aW;
}

double complex q0_bil_sdqH_zeta(double zeta,void *tmp)
{
  double complex q0_bil_sdqH_eta(double eta,void *tmp);

  TDATA *td=(TDATA *)tmp;
  double complex res;
  double err;

  td->zeta=zeta;
  res=deintz(q0_bil_sdqH_eta,-1.0,1.0,td,IEPS,&err);
  if(err<0.0){ printf("q0_bil_dcoef.c, q0_bil_sdqH_zeta(), q0_bil_sdqH_eta() DE integration error! err=%g. Exit...\n",err);  exit(0);  }
  return res;
}

double complex q0_bil_sdqH_eta(double eta,void *tmp)
{
  TDATA *td=(TDATA *)tmp;
  double complex qG,dqG[3],d2qG[6];
  double vW[3],vR[3];
  int i,err;

  bil_rw_zeta_eta(vR,vW,td->zeta,eta,td->cr,td->cw);
  for(i=0;i<3;i++) vR[i]-=td->rt[i];
  err=d3hm_qpgf_d2_v2(&qG,dqG,d2qG,vR,MEPS,td->qpd);
  if(err<0){
    printf("q0_bil_dcoef.c, q0_bil_sdqH_eta(), d3hm_qpgf_d2_v2() error. err=%d. Exit...\n",err);
    exit(0);
  }

  return ( d2qG[0]*vW[0]*td->vt[0]
          +d2qG[1]*vW[1]*td->vt[1]
          +d2qG[2]*vW[2]*td->vt[2]
          +d2qG[3]*(vW[0]*td->vt[1]+vW[1]*td->vt[0])
          +d2qG[4]*(vW[1]*td->vt[2]+vW[2]*td->vt[1])
          +d2qG[5]*(vW[2]*td->vt[0]+vW[0]*td->vt[2]) )*bil_Nn(td->nid,td->zeta,eta);
}

