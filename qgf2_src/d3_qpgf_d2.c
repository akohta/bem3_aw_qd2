/*
 * d3_qpgf_d2_v2.c
 *
 *  Created on: Oct 14, 2019
 *      Author: ohta
 */

#include "d3_qpgf_d2.h"

int d3hm_qpgf_d2_qG(double complex *qG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_qG_ew(double complex *qG,double *r,double eps,QPD2 *qd);
  int d3hm_qpgf_d2_qG_fd(double complex *qG,double *r,double eps,QPD2 *qd);
  
  double bN,aNx,aNy,k2,idd,aNpx,aNmx,aNpy,aNmy,ad12,ad22,cpd,bNp,bNm;
  int n1,n2,err;

  k2=qd->k*qd->k;
  idd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  cpd=2.0*M_PI*idd;
  ad12=qd->vd1[0]*qd->vd1[0]+qd->vd1[1]*qd->vd1[1];
  ad22=qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1];

  n1=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );
  n2=(int)round( (qd->vk[0]*qd->vd2[0]+qd->vk[1]*qd->vd2[1])/(2.0*M_PI) );

  aNx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
  aNy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  bN=sqrt(k2-aNx*aNx-aNy*aNy);

  if(ad22<ad12){ // |g_1| < |g_2|
    aNpx=-qd->vk[0]+(double)(n1+S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)(n1+S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)(n1-S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)(n1-S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  }
  else { // |g_1| >= |g_2|
    aNpx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2+S_LIMIT)*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2+S_LIMIT)*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2-S_LIMIT)*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2-S_LIMIT)*cpd*qd->vd1[0];
  }
  bNp=k2-aNpx*aNpx-aNpy*aNpy;
  bNm=k2-aNmx*aNmx-aNmy*aNmy;

  if(bNp>0.0 || bNm>0.0){ // Ewald method
    err=d3hm_qpgf_d2_qG_ew(qG,r,eps,qd);
    if(err<0) err-=20;
    return err;
  }
  else {
    bNp=sqrt(-bNp);
    bNm=sqrt(-bNm);
    if( bN/bNp*exp(-bNp*fabs(r[2]))<eps && bN/bNm*exp(-bNm*fabs(r[2]))<eps ){ // Fourier domain sum
      err=d3hm_qpgf_d2_qG_fd(qG,r,eps,qd);
      if(err<0) err-=10;
      return err;
    }
    else { // Ewald method
      err=d3hm_qpgf_d2_qG_ew(qG,r,eps,qd);
      if(err<0) err-=20;
      return err;
    }
  }
}

int d3hm_qpgf_d2_dqG(double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_dqG_ew(double complex *dqG,double *r,double eps,QPD2 *qd);
  int d3hm_qpgf_d2_dqG_fd(double complex *dqG,double *r,double eps,QPD2 *qd);
  
  double bN,aNx,aNy,k2,idd,aNpx,aNmx,aNpy,aNmy,ad12,ad22,cpd,bNp,bNm;
  int n1,n2,err;

  k2=qd->k*qd->k;
  idd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  cpd=2.0*M_PI*idd;
  ad12=qd->vd1[0]*qd->vd1[0]+qd->vd1[1]*qd->vd1[1];
  ad22=qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1];

  n1=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );
  n2=(int)round( (qd->vk[0]*qd->vd2[0]+qd->vk[1]*qd->vd2[1])/(2.0*M_PI) );

  aNx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
  aNy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  bN=sqrt(k2-aNx*aNx-aNy*aNy);

  if(ad22<ad12){ // |g_1| < |g_2|
    aNpx=-qd->vk[0]+(double)(n1+S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)(n1+S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)(n1-S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)(n1-S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  }
  else { // |g_1| >= |g_2|
    aNpx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2+S_LIMIT)*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2+S_LIMIT)*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2-S_LIMIT)*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2-S_LIMIT)*cpd*qd->vd1[0];
  }
  bNp=k2-aNpx*aNpx-aNpy*aNpy;
  bNm=k2-aNmx*aNmx-aNmy*aNmy;

  if(bNp>0.0 || bNm>0.0){ // Ewald method
    err=d3hm_qpgf_d2_dqG_ew(dqG,r,eps,qd);
    if(err<0) err-=20;
    return err;
  }
  else {
    bNp=sqrt(-bNp);
    bNm=sqrt(-bNm);
    if( bN/bNp*exp(-bNp*fabs(r[2]))<eps && bN/bNm*exp(-bNm*fabs(r[2]))<eps ){ // Fourier domain sum
      err=d3hm_qpgf_d2_dqG_fd(dqG,r,eps,qd);
      if(err<0) err-=10;
      return err;
    }
    else { // Ewald method
      err=d3hm_qpgf_d2_dqG_ew(dqG,r,eps,qd);
      if(err<0) err-=20;
      return err;
    }
  }
}

int d3hm_qpgf_d2_v1(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
  double bN,aNx,aNy,k2,idd,aNpx,aNmx,aNpy,aNmy,ad12,ad22,cpd,bNp,bNm;
  int n1,n2,err;

  k2=qd->k*qd->k;
  idd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  cpd=2.0*M_PI*idd;
  ad12=qd->vd1[0]*qd->vd1[0]+qd->vd1[1]*qd->vd1[1];
  ad22=qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1];

  n1=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );
  n2=(int)round( (qd->vk[0]*qd->vd2[0]+qd->vk[1]*qd->vd2[1])/(2.0*M_PI) );

  aNx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
  aNy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  bN=sqrt(k2-aNx*aNx-aNy*aNy);

  if(ad22<ad12){ // |g_1| < |g_2|
    aNpx=-qd->vk[0]+(double)(n1+S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)(n1+S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)(n1-S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)(n1-S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  }
  else { // |g_1| >= |g_2|
    aNpx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2+S_LIMIT)*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2+S_LIMIT)*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2-S_LIMIT)*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2-S_LIMIT)*cpd*qd->vd1[0];
  }
  bNp=k2-aNpx*aNpx-aNpy*aNpy;
  bNm=k2-aNmx*aNmx-aNmy*aNmy;

  if(bNp>0.0 || bNm>0.0){ // Ewald method
    err=d3hm_qpgf_d2_v1_ew(qG,dqG,r,eps,qd);
    if(err<0) err-=20;
    return err;
  }
  else {
    bNp=sqrt(-bNp);
    bNm=sqrt(-bNm);
    if( bN/bNp*exp(-bNp*fabs(r[2]))<eps && bN/bNm*exp(-bNm*fabs(r[2]))<eps ){ // Fourier domain sum
      err=d3hm_qpgf_d2_v1_fd(qG,dqG,r,eps,qd);
      if(err<0) err-=10;
      return err;
    }
    else { // Ewald method
      err=d3hm_qpgf_d2_v1_ew(qG,dqG,r,eps,qd);
      if(err<0) err-=20;
      return err;
    }
  }
}

int d3hm_qpgf_d2_v2(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd)
{
  double bN,aNx,aNy,k2,idd,aNpx,aNmx,aNpy,aNmy,ad12,ad22,cpd,bNp,bNm;
  int n1,n2,err;

  k2=qd->k*qd->k;
  idd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  cpd=2.0*M_PI*idd;
  ad12=qd->vd1[0]*qd->vd1[0]+qd->vd1[1]*qd->vd1[1];
  ad22=qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1];

  n1=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );
  n2=(int)round( (qd->vk[0]*qd->vd2[0]+qd->vk[1]*qd->vd2[1])/(2.0*M_PI) );

  aNx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
  aNy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  bN=sqrt(k2-aNx*aNx-aNy*aNy);

  if(ad22<ad12){ // |g_1| < |g_2|
    aNpx=-qd->vk[0]+(double)(n1+S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)(n1+S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)(n1-S_LIMIT)*cpd*qd->vd2[1]-(double)n2*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)(n1-S_LIMIT)*cpd*qd->vd2[0]+(double)n2*cpd*qd->vd1[0];
  }
  else { // |g_1| >= |g_2|
    aNpx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2+S_LIMIT)*cpd*qd->vd1[1];
    aNpy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2+S_LIMIT)*cpd*qd->vd1[0];
    aNmx=-qd->vk[0]+(double)n1*cpd*qd->vd2[1]-(double)(n2-S_LIMIT)*cpd*qd->vd1[1];
    aNmy=-qd->vk[1]-(double)n1*cpd*qd->vd2[0]+(double)(n2-S_LIMIT)*cpd*qd->vd1[0];
  }
  bNp=k2-aNpx*aNpx-aNpy*aNpy;
  bNm=k2-aNmx*aNmx-aNmy*aNmy;

  if(bNp>0.0 || bNm>0.0){ // Ewald method
    err=d3hm_qpgf_d2_v2_ew(qG,dqG,d2qG,r,eps,qd);
    if(err<0) err-=20;
    return err;
  }
  else {
    bNp=sqrt(-bNp);
    bNm=sqrt(-bNm);
    if( bN/bNp*exp(-bNp*fabs(r[2]))<eps && bN/bNm*exp(-bNm*fabs(r[2]))<eps ){ // Fourier domain sum
      err=d3hm_qpgf_d2_v2_fd(qG,dqG,d2qG,r,eps,qd);
      if(err<0) err-=10;
      return err;
    }
    else { // Ewald method
      err=d3hm_qpgf_d2_v2_ew(qG,dqG,d2qG,r,eps,qd);
      if(err<0) err-=20;
      return err;
    }
  }
}

int d3hm_qpgf_d2_qG_fd(double complex *qG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_qG_fd_sd(double complex *qG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2_qG_fd_sd(qG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_dqG_fd(double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_dqG_fd_sd(double complex *dqG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2_dqG_fd_sd(dqG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    dqG[0]*=cf;
    dqG[1]*=cf;
    dqG[2]*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_v1_fd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2v1_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2v1_fd_sd(qG,dqG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
    dqG[0]*=cf;
    dqG[1]*=cf;
    dqG[2]*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_v2_fd(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2v2_fd_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err,i;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2v2_fd_sd(qG,dqG,d2qG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
    for(i=0;i<3;i++) dqG[i]*=cf;
    for(i=0;i<6;i++) d2qG[i]*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_qG_ew(double complex *qG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_qG_ew_sd(double complex *qG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2_qG_ew_sd(qG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_dqG_ew(double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2_dqG_ew_sd(double complex *dqG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2_dqG_ew_sd(dqG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    dqG[0]*=cf;
    dqG[1]*=cf;
    dqG[2]*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_v1_ew(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2v1_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2v1_ew_sd(qG,dqG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
    dqG[0]*=cf;
    dqG[1]*=cf;
    dqG[2]*=cf;
  }

  return err;
}

int d3hm_qpgf_d2_v2_ew(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd)
{
  int d3hm_qpgf_d2v2_ew_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd);

  double complex cf;
  double rt[3],i_detd,arg;
  int l1,l2,err,i;

  i_detd=1.0/(qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1]);
  l1=(int)floor( (r[0]*qd->vd2[1]-r[1]*qd->vd2[0])*i_detd+0.5 );
  l2=(int)floor(-(r[0]*qd->vd1[1]-r[1]*qd->vd1[0])*i_detd+0.5 );

  rt[0]=r[0]-(double)l1*qd->vd1[0]-(double)l2*qd->vd2[0];
  rt[1]=r[1]-(double)l1*qd->vd1[1]-(double)l2*qd->vd2[1];
  rt[2]=r[2];

  err=d3hm_qpgf_d2v2_ew_sd(qG,dqG,d2qG,rt,eps,qd);

  if(l1!=0 || l2!=0){
    arg=-(qd->vk[0]*((double)l1*qd->vd1[0]+(double)l2*qd->vd2[0])+qd->vk[1]*((double)l1*qd->vd1[1]+(double)l2*qd->vd2[1]));
    cf=cos(arg)+I*sin(arg);

    *qG*=cf;
    for(i=0;i<3;i++) dqG[i]*=cf;
    for(i=0;i<6;i++) d2qG[i]*=cf;
  }

  return err;
}

//////////////////////////////////////////////////////////////////////
int d3hm_sf_comp_z(double complex *v1,double complex *v2,double eps)
{ // return 1 : same, return 0 : different.
  if(    2.0*fabs(creal(*v1)-creal(*v2))<=eps*(fabs(creal(*v1))+fabs(creal(*v2)))
      && 2.0*fabs(cimag(*v1)-cimag(*v2))<=eps*(fabs(cimag(*v1))+fabs(cimag(*v2))) ) return 1;
  else return 0;
}

int d3hm_sf_diff_z(double complex *v1,double complex *v2,double eps)
{ // return 1 : |real(v2*eps)|=>|real(v1)| && |imag(v2*veps)|=>|imag(v1)|, return 0 : others
  if( fabs(creal(*v1))<=fabs(creal(*v2))*EFS*eps && fabs(cimag(*v1))<=fabs(cimag(*v2))*EFS*eps) return 1;
  else return 0;
}

int d3hm_qpgf_d2_qG_fd_sd(double complex *qG,double *r,double eps,QPD2 *qd)
{
  int fd_sum_qG_n2(double complex *tqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err);

  double complex oqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,err,cc,nc;

  if(r[2]==0.0) return -4; // z==0.0

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -3; // determinant d == 0
  idd=1.0/detd;

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  nc=(int)round((qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI));

  cc=0;
  *qG=0.0;

  n1=nc;
  cc+=fd_sum_qG_n2(qG,n1,vg1,vg2,r,eps,qd,&err);  if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+FD_LIMIT;n1++){
    oqG=*qG;
    cc+=fd_sum_qG_n2(qG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;

    if( d3hm_sf_comp_z(qG,&oqG,eps) ) break;
  }
  if(n1==nc+FD_LIMIT) return -1; // summation limit

  // n1 < nc
  for(n1=nc-1;n1>nc-FD_LIMIT;n1--){
    oqG=*qG;
    cc+=fd_sum_qG_n2(qG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;

    if( d3hm_sf_comp_z(qG,&oqG,eps) ) break;
  }
  if(n1==nc-FD_LIMIT) return -1; // summation limit

  *qG*=0.5*I*fabs(idd);

  return cc;
}

int fd_sum_qG_n2(double complex *tqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err)
{
  int fd_qG_cFn(double complex *Fn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd);

  double complex Fn;
  int n2,cc,nc;

  nc=(int)round( (qd->vk[0]*vg2[0]+qd->vk[1]*vg2[1]-(double)n1*(vg1[0]*vg2[0]+vg1[1]*vg2[1]))/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  cc=0;
  n2=nc;
  *err=fd_qG_cFn(&Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
  *tqG+=Fn;
  cc+=1;

  // n2 > nc
  for(n2=nc+1;n2<nc+FD_LIMIT;n2++){
    *err=fd_qG_cFn(&Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&Fn,tqG,eps) ) break;
    }
  }
  if(n2==nc+FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-FD_LIMIT;n2--){
    *err=fd_qG_cFn(&Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&Fn,tqG,eps) ) break;
    }
  }
  if(n2==nc-FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  return cc;
}

int fd_qG_cFn(double complex *Fn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd)
{
  double an[2],bn2,bn,arg;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];

  if(bn2>0.0){
    bn=sqrt(bn2);
    arg=an[0]*r[0]+an[1]*r[1]+bn*fabs(r[2]);
    *Fn=(cos(arg)+I*sin(arg))/bn;
    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    arg=an[0]*r[0]+an[1]*r[1];
    *Fn=-I*(cos(arg)+I*sin(arg))*exp(-bn*fabs(r[2]))/bn;
    return 1;
  }
  else {
    return -2; // bn==0.0
  }
}

int d3hm_qpgf_d2_dqG_fd_sd(double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int fd_sum_dqG_n2(double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err);

  double complex odqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,err,cc,nc;

  if(r[2]==0.0) return -4; // z==0.0

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -3; // determinant d == 0
  idd=1.0/detd;

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  nc=(int)round((qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI));

  cc=0;
  dqG[0]=0.0;
  dqG[1]=0.0;
  dqG[2]=0.0;

  n1=nc;
  cc+=fd_sum_dqG_n2(dqG,n1,vg1,vg2,r,eps,qd,&err);  if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+FD_LIMIT;n1++){
    odqG=dqG[2];
    cc+=fd_sum_dqG_n2(dqG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&dqG[2],&odqG,eps) ) break;
  }
  if(n1==nc+FD_LIMIT) return -1; // summation limit

  // n1 < nc
  for(n1=nc-1;n1>nc-FD_LIMIT;n1--){
    odqG=dqG[2];
    cc+=fd_sum_dqG_n2(dqG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&dqG[2],&odqG,eps) ) break;
  }
  if(n1==nc-FD_LIMIT) return -1; // summation limit

  dqG[0]*=-0.5*fabs(idd);
  dqG[1]*=-0.5*fabs(idd);
  dqG[2]*=-0.5*fabs(idd)*r[2]/fabs(r[2]);

  return cc;
}

int fd_sum_dqG_n2(double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err)
{
  int fd_dqG_cFn(double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd);

  double complex dFn[3];
  int n2,cc,nc;

  nc=(int)round( (qd->vk[0]*vg2[0]+qd->vk[1]*vg2[1]-(double)n1*(vg1[0]*vg2[0]+vg1[1]*vg2[1]))/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  cc=0;
  n2=nc;
  *err=fd_dqG_cFn(dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
  tdqG[0]+=dFn[0];
  tdqG[1]+=dFn[1];
  tdqG[2]+=dFn[2];
  cc+=1;

  // n2 > nc
  for(n2=nc+1;n2<nc+FD_LIMIT;n2++){
    *err=fd_dqG_cFn(dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
  }
  if(n2==nc+FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-FD_LIMIT;n2--){
    *err=fd_dqG_cFn(dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
  }
  if(n2==nc-FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  return cc;
}

int fd_dqG_cFn(double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd)
{
  double complex Fn;
  double an[2],bn2,bn,arg;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];

  if(bn2>0.0){
    bn=sqrt(bn2);
    arg=an[0]*r[0]+an[1]*r[1]+bn*fabs(r[2]);
    dFn[2]=cos(arg)+I*sin(arg);
    Fn=dFn[2]/bn;
    dFn[0]=Fn*an[0];
    dFn[1]=Fn*an[1];
    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    arg=an[0]*r[0]+an[1]*r[1];
    dFn[2]=(cos(arg)+I*sin(arg))*exp(-bn*fabs(r[2]));
    Fn=-I*dFn[2]/bn;
    dFn[0]=Fn*an[0];
    dFn[1]=Fn*an[1];
    return 1;
  }
  else {
    return -2; // bn==0.0
  }
}

int d3hm_qpgf_d2v1_fd_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int fd_sum_v1_n2(double complex *tqG,double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err);

  double complex odqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,err,cc,nc;

  if(r[2]==0.0) return -4; // z==0.0

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -3; // determinant d == 0
  idd=1.0/detd;

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  nc=(int)round((qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI));

  cc=0;
  *qG=0.0;
  dqG[0]=0.0;
  dqG[1]=0.0;
  dqG[2]=0.0;

  n1=nc;
  cc+=fd_sum_v1_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);  if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+FD_LIMIT;n1++){
    odqG=dqG[2];
    cc+=fd_sum_v1_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&dqG[2],&odqG,eps) ) break;
  }
  if(n1==nc+FD_LIMIT) return -1; // summation limit

  // n1 < nc
  for(n1=nc-1;n1>nc-FD_LIMIT;n1--){
    odqG=dqG[2];
    cc+=fd_sum_v1_n2(qG,dqG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&dqG[2],&odqG,eps) ) break;
  }
  if(n1==nc-FD_LIMIT) return -1; // summation limit

  *qG*=0.5*I*fabs(idd);
  dqG[0]*=-0.5*fabs(idd);
  dqG[1]*=-0.5*fabs(idd);
  dqG[2]*=-0.5*fabs(idd)*r[2]/fabs(r[2]);

  return cc;
}

int fd_sum_v1_n2(double complex *tqG,double complex *tdqG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err)
{
  int fdv1_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd);

  double complex Fn,dFn[3];
  int n2,cc,nc;

  nc=(int)round( (qd->vk[0]*vg2[0]+qd->vk[1]*vg2[1]-(double)n1*(vg1[0]*vg2[0]+vg1[1]*vg2[1]))/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  cc=0;
  n2=nc;
  *err=fdv1_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
  *tqG+=Fn;
  tdqG[0]+=dFn[0];
  tdqG[1]+=dFn[1];
  tdqG[2]+=dFn[2];
  cc+=1;

  // n2 > nc
  for(n2=nc+1;n2<nc+FD_LIMIT;n2++){
    *err=fdv1_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
  }
  if(n2==nc+FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-FD_LIMIT;n2--){
    *err=fdv1_cFn(&Fn,dFn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
  }
  if(n2==nc-FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  return cc;
}

int fdv1_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd)
{
  double an[2],bn2,bn,arg;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];

  if(bn2>0.0){
    bn=sqrt(bn2);
    arg=an[0]*r[0]+an[1]*r[1]+bn*fabs(r[2]);
    dFn[2]=cos(arg)+I*sin(arg);
    *Fn=dFn[2]/bn;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    arg=an[0]*r[0]+an[1]*r[1];
    dFn[2]=(cos(arg)+I*sin(arg))*exp(-bn*fabs(r[2]));
    *Fn=-I*dFn[2]/bn;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    return 1;
  }
  else {
    return -2; // bn==0.0
  }
}

int d3hm_qpgf_d2v2_fd_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd)
{
  int fd_sum_v2_n2(double complex *tqG,double complex *tdqG,double complex *td2qG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err);

  double complex od2qG;
  double vg1[2],vg2[2],detd,idd,sgnz;
  int n1,err,cc,nc,i;

  if(r[2]==0.0) return -4; // z==0.0

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -3; // determinant d == 0
  idd=1.0/detd;

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  nc=(int)round((qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI));

  cc=0;
  *qG=0.0;
  for(i=0;i<3;i++) dqG[i]=0.0;
  for(i=0;i<6;i++) d2qG[i]=0.0;

  n1=nc;
  cc+=fd_sum_v2_n2(qG,dqG,d2qG,n1,vg1,vg2,r,eps,qd,&err);  if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+FD_LIMIT;n1++){
    od2qG=d2qG[2];
    cc+=fd_sum_v2_n2(qG,dqG,d2qG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&d2qG[2],&od2qG,eps) ) break;
  }
  if(n1==nc+FD_LIMIT) return -1; // summation limit

  // n1 < nc
  for(n1=nc-1;n1>nc-FD_LIMIT;n1--){
    od2qG=d2qG[2];
    cc+=fd_sum_v2_n2(qG,dqG,d2qG,n1,vg1,vg2,r,eps,qd,&err);    if(err<0) return err;
    if( d3hm_sf_comp_z(&d2qG[2],&od2qG,eps) ) break;
  }
  if(n1==nc-FD_LIMIT) return -1; // summation limit

  sgnz=r[2]/fabs(r[2]);

  *qG*=0.5*I*fabs(idd);
  dqG[0]*=-0.5*fabs(idd);
  dqG[1]*=-0.5*fabs(idd);
  dqG[2]*=-0.5*fabs(idd)*sgnz;
  d2qG[0]*=-0.5*I*fabs(idd);
  d2qG[1]*=-0.5*I*fabs(idd);
  d2qG[2]*=-0.5*I*fabs(idd);
  d2qG[3]*=-0.5*I*fabs(idd);
  d2qG[4]*=-0.5*I*fabs(idd)*sgnz;
  d2qG[5]*=-0.5*I*fabs(idd)*sgnz;

  return cc;
}

int fd_sum_v2_n2(double complex *tqG,double complex *tdqG,double complex *td2qG,int n1,double *vg1,double *vg2,double *r,double eps,QPD2 *qd,int *err)
{
  int fdv2_cFn(double complex *Fn,double complex *dFn,double complex *d2Fn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd);

  double complex Fn,dFn[3],d2Fn[6];
  int n2,cc,nc,i;

  nc=(int)round( (qd->vk[0]*vg2[0]+qd->vk[1]*vg2[1]-(double)n1*(vg1[0]*vg2[0]+vg1[1]*vg2[1]))/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  cc=0;
  n2=nc;
  *err=fdv2_cFn(&Fn,dFn,d2Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
  *tqG+=Fn;
  for(i=0;i<3;i++) tdqG[i]+=dFn[i];
  for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];
  cc+=1;

  // n2 > nc
  for(n2=nc+1;n2<nc+FD_LIMIT;n2++){
    *err=fdv2_cFn(&Fn,dFn,d2Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    for(i=0;i<3;i++) tdqG[i]+=dFn[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&d2Fn[2],&td2qG[2],eps) ) break;
    }
  }
  if(n2==nc+FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-FD_LIMIT;n2--){
    *err=fdv2_cFn(&Fn,dFn,d2Fn,n1,n2,vg1,vg2,r,qd); if(*err<0) return 0; // bn==0
    *tqG+=Fn;
    for(i=0;i<3;i++) tdqG[i]+=dFn[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];
    cc+=1;

    if(*err>0){
      if( d3hm_sf_diff_z(&d2Fn[2],&td2qG[2],eps) ) break;
    }
  }
  if(n2==nc-FD_LIMIT){
    *err=-1; // summation limit
    return 0;
  }

  return cc;
}

int fdv2_cFn(double complex *Fn,double complex *dFn,double complex *d2Fn,int n1,int n2,double *vg1,double *vg2,double *r,QPD2 *qd)
{
  double an[2],bn2,bn,arg;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];

  if(bn2>0.0){
    bn=sqrt(bn2);
    arg=an[0]*r[0]+an[1]*r[1]+bn*fabs(r[2]);
    dFn[2]=cos(arg)+I*sin(arg);
    *Fn=dFn[2]/bn;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    d2Fn[0]=an[0]*dFn[0];
    d2Fn[1]=an[1]*dFn[1];
    d2Fn[2]=bn2**Fn;
    d2Fn[3]=an[0]*dFn[1];
    d2Fn[4]=an[1]*dFn[2];
    d2Fn[5]=an[0]*dFn[2];
    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2);
    arg=an[0]*r[0]+an[1]*r[1];
    dFn[2]=(cos(arg)+I*sin(arg))*exp(-bn*fabs(r[2]));
    *Fn=-I*dFn[2]/bn;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    d2Fn[0]=an[0]*dFn[0];
    d2Fn[1]=an[1]*dFn[1];
    d2Fn[2]=bn2**Fn;
    d2Fn[3]=an[0]*dFn[1];
    d2Fn[4]=an[1]*dFn[2];
    d2Fn[5]=an[0]*dFn[2];
    return 1;
  }
  else {
    return -2; // bn==0.0
  }
}

int d3hm_qpgf_d2_qG_ew_sd(double complex *qG,double *r,double eps,QPD2 *qd)
{
  int ew_qG1(double complex *qG1,double *r,double eps,double veps,QPD2 *qd);
  int ew_qG2(double complex *qG2,double *r,double eps,double veps,QPD2 *qd);

  double complex qG1,qG2;
  double detd,veps,adp,adm,r2_min;
  int c1,c2,l1,l2;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0

  adp=sqrt(pow(qd->vd1[0]+qd->vd2[0],2)+pow(qd->vd1[1]+qd->vd2[1],2));
  adm=sqrt(pow(qd->vd1[0]-qd->vd2[0],2)+pow(qd->vd1[1]-qd->vd2[1],2));
  veps=sqrt(fabs(detd)/(4.0*M_PI)*adp/adm);

  l1=0;
  l2=0;
  r2_min=pow(r[0]+(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0],2)
      +  pow(r[1]+(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1],2)+r[2]*r[2];

  if( qd->k*qd->k*veps*veps-r2_min/(4.0*veps*veps)>EW_EXPM )
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r2_min+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  c1=ew_qG1(&qG1,r,eps,veps,qd);
  *qG=qG1;
  if(c1<0) return c1;

  // qG2
  c2=ew_qG2(&qG2,r,eps,veps,qd);
  *qG+=qG2;
  if(c2<0) return c2;

  return c1+c2;
}

int ew_qG1(double complex *qG1,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG1_sum_l2(double complex *tqG,int l1,double *r,double eps,double veps,QPD2 *qd);

  double complex oqG;
  int l1,cc,lc;

  lc=0;

  cc=0;
  *qG1=0.0;

  l1=lc;
  cc+=ew_qG1_sum_l2(qG1,l1,r,eps,veps,qd);  if(cc<0) return cc;

  // l1 > lc
  for(l1=lc+1;l1<lc+EW_LIMIT;l1++){
    oqG=*qG1;
    cc+=ew_qG1_sum_l2(qG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(qG1,&oqG,eps) ) break;
  }
  if(l1==lc+EW_LIMIT){
    return -1;
  }

  // l1 < lc
  for(l1=lc-1;l1>lc-EW_LIMIT;l1--){
    oqG=*qG1;

    cc+=ew_qG1_sum_l2(qG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(qG1,&oqG,eps) ) break;
  }
  if(l1==lc-EW_LIMIT){
    return -1;
  }

  *qG1/=4.0*M_PI;

  return cc;
}

int ew_qG1_sum_l2(double complex *tqG,int l1,double *r,double eps,double veps,QPD2 *qd)
{
  void ew_qG1_cFl(double complex *Fl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd);

  double complex Fl;
  int cc,l2,lc;

  lc=-(int)round( ((r[0]+(double)l1*qd->vd1[0])*qd->vd2[0]+(r[1]+(double)l1*qd->vd1[1])*qd->vd2[1])
                    /(qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1]) );

  cc=0;
  l2=lc;
  cc+=1;
  ew_qG1_cFl(&Fl,l1,l2,r,eps,veps,qd);
  *tqG+=Fl;

  // l2 > lc
  for(l2=lc+1;l2<lc+EW_LIMIT;l2++){
    cc+=1;
    ew_qG1_cFl(&Fl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;

    if( d3hm_sf_diff_z(&Fl,tqG,eps) ) break;
  }
  if(l2==lc+EW_LIMIT){
    return -1;
  }

  // l2 < 0
  for(l2=lc-1;l2>lc-EW_LIMIT;l2--){
    cc+=1;
    ew_qG1_cFl(&Fl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;

    if( d3hm_sf_diff_z(&Fl,tqG,eps) ) break;
  }
  if(l2==lc-EW_LIMIT){
    return -1;
  }

  return cc;
}

void ew_qG1_cFl(double complex *Fl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd)
{
  double complex w,ce;
  double rl,i_rl,Rl[2],ke,i_ve,rl_2e,arg,ef;

  Rl[0]=(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0];
  Rl[1]=(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1];
  rl=sqrt((r[0]+Rl[0])*(r[0]+Rl[0])+(r[1]+Rl[1])*(r[1]+Rl[1])+r[2]*r[2]);
  i_rl=1.0/rl;

  ke=qd->k*veps;
  i_ve=1.0/veps;
  rl_2e=rl*0.5*i_ve;

  arg=qd->vk[0]*Rl[0]+qd->vk[1]*Rl[1];
  ce=cos(arg)+I*sin(arg);
  ef=exp(ke*ke-rl_2e*rl_2e);
  w=Faddeeva_w(ke+I*rl_2e,eps);

  *Fl=ce*i_rl*ef*creal(w);
}

int ew_qG2(double complex *qG2,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG2_sum_n2(double complex *tqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err);

  double complex oqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,cc,err,nc;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0
  idd=1.0/detd;

  nc=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  cc=0;
  *qG2=0.0;

  n1=nc;
  cc+=ew_qG2_sum_n2(qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+EW_LIMIT;n1++){
    oqG=*qG2;
    cc+=ew_qG2_sum_n2(qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0) if( d3hm_sf_comp_z(qG2,&oqG,eps) ) break;
  }
  if(n1==nc+EW_LIMIT){
    return -2;
  }

  // n1 < nc
  for(n1=nc-1;n1>nc-EW_LIMIT;n1--){
    oqG=*qG2;
    cc+=ew_qG2_sum_n2(qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0) if( d3hm_sf_comp_z(qG2,&oqG,eps) ) break;
  }
  if(n1==nc-EW_LIMIT){
    return -2;
  }

  *qG2*=0.25*I*fabs(idd);

  return cc;
}

int ew_qG2_sum_n2(double complex *tqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err)
{
  int ew_qG2_cFn(double complex *Fn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd);

  double complex Fn;
  int cc,n2,nc,flg;

  nc=(int)round(-((-qd->vk[0]+(double)n1*vg1[0])*vg2[0]+(-qd->vk[1]+(double)n1*vg1[1])*vg2[1])/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  flg=0;
  cc=0;
  n2=nc;
  cc+=1;
  *err=ew_qG2_cFn(&Fn,n1,n2,r,vg1,vg2,eps,veps,qd);  if(*err<0) return 0;
  *tqG+=Fn;

  // n2 > nc
  for(n2=nc+1;n2<nc+EW_LIMIT;n2++){
    cc+=1;
    *err=ew_qG2_cFn(&Fn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;

    if(*err>0){
      if( d3hm_sf_diff_z(&Fn,tqG,eps) ) break;
    }
    else flg++;
  }
  if(n2==nc+EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-EW_LIMIT;n2--){
    cc+=1;
    *err=ew_qG2_cFn(&Fn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;

    if(*err>0){
      if( d3hm_sf_diff_z(&Fn,tqG,eps) ) break;
    }
    else flg++;
  }
  if(n2==nc-EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  if(flg==0) *err=1;
  else *err=0;

  return cc;
}

int ew_qG2_cFn(double complex *Fn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd)
{
  double complex ce,w,ff;
  double an[2],bn2,bn,arg,azve,bnaz,ebz,ef,f1,f2;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];
  arg=an[0]*r[0]+an[1]*r[1];
  ce=cos(arg)+I*sin(arg);

  if(bn2>0.0){
    bn=sqrt(bn2);
    azve=r[2]/(2.0*veps);
    arg=bn*r[2];
    ef=exp(bn2*veps*veps-azve*azve);
    ff=cos(arg)+I*sin(arg);
    w=Faddeeva_w(bn*veps+I*azve,eps);

    *Fn=2.0*ce/bn*(I*ef*cimag(w)+ff);

    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2); // |bn|
    azve=r[2]/(2.0*veps);
    bnaz=bn*r[2];
    ebz=exp(bnaz);
    f1=ebz*erfc(bn*veps+azve);
    f2=1.0/ebz*erfc(bn*veps-azve);

    *Fn=ce/(I*bn)*(f1+f2);

    return 1;
  }
  else { // bn=0.0;
    return -3;
  }
}

int d3hm_qpgf_d2_dqG_ew_sd(double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int ew_dqG1(double complex *dqG1,double *r,double eps,double veps,QPD2 *qd);
  int ew_dqG2(double complex *dqG2,double *r,double eps,double veps,QPD2 *qd);

  double complex dqG1[3],dqG2[3];
  double detd,veps,adp,adm,r2_min;
  int c1,c2,l1,l2;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0

  adp=sqrt(pow(qd->vd1[0]+qd->vd2[0],2)+pow(qd->vd1[1]+qd->vd2[1],2));
  adm=sqrt(pow(qd->vd1[0]-qd->vd2[0],2)+pow(qd->vd1[1]-qd->vd2[1],2));
  veps=sqrt(fabs(detd)/(4.0*M_PI)*adp/adm);

  l1=0;
  l2=0;
  r2_min=pow(r[0]+(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0],2)
      +  pow(r[1]+(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1],2)+r[2]*r[2];

  if( qd->k*qd->k*veps*veps-r2_min/(4.0*veps*veps)>EW_EXPM )
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r2_min+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  c1=ew_dqG1(dqG1,r,eps,veps,qd);
  dqG[0]=dqG1[0];
  dqG[1]=dqG1[1];
  dqG[2]=dqG1[2];
  if(c1<0) return c1;

  // qG2
  c2=ew_dqG2(dqG2,r,eps,veps,qd);
  dqG[0]+=dqG2[0];
  dqG[1]+=dqG2[1];
  dqG[2]+=dqG2[2];
  if(c2<0) return c2;

  return c1+c2;
}

int ew_dqG1(double complex *dqG1,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_dqG1_sum_l2(double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd);

  double complex odqG[3];
  int l1,cc,lc;

  lc=0;

  cc=0;
  dqG1[0]=0.0;
  dqG1[1]=0.0;
  dqG1[2]=0.0;

  l1=lc;
  cc+=ew_dqG1_sum_l2(dqG1,l1,r,eps,veps,qd);  if(cc<0) return cc;

  // l1 > lc
  for(l1=lc+1;l1<lc+EW_LIMIT;l1++){
    odqG[0]=dqG1[0];
    odqG[1]=dqG1[1];
    odqG[2]=dqG1[2];
    cc+=ew_dqG1_sum_l2(dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(&dqG1[0],&odqG[0],eps) )
      if( d3hm_sf_comp_z(&dqG1[1],&odqG[1],eps) ) break;
  }
  if(l1==lc+EW_LIMIT){
    return -1;
  }

  // l1 < lc
  for(l1=lc-1;l1>lc-EW_LIMIT;l1--){
    odqG[0]=dqG1[0];
    odqG[1]=dqG1[1];
    odqG[2]=dqG1[2];
    cc+=ew_dqG1_sum_l2(dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(&dqG1[0],&odqG[0],eps) )
      if( d3hm_sf_comp_z(&dqG1[1],&odqG[1],eps) ) break;
  }
  if(l1==lc-EW_LIMIT){
    return -1;
  }

  dqG1[0]/=-4.0*M_PI;
  dqG1[1]/=-4.0*M_PI;
  dqG1[2]*=-r[2]/(4.0*M_PI);

  return cc;
}

int ew_dqG1_sum_l2(double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd)
{
  void ew_dqG1_cFl(double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd);

  double complex dFl[3];
  int cc,l2,lc;

  lc=-(int)round( ((r[0]+(double)l1*qd->vd1[0])*qd->vd2[0]+(r[1]+(double)l1*qd->vd1[1])*qd->vd2[1])
                    /(qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1]) );

  cc=0;
  l2=lc;
  cc+=1;
  ew_dqG1_cFl(dFl,l1,l2,r,eps,veps,qd);
  tdqG[0]+=dFl[0];
  tdqG[1]+=dFl[1];
  tdqG[2]+=dFl[2];

  // l2 > lc
  for(l2=lc+1;l2<lc+EW_LIMIT;l2++){
    cc+=1;
    ew_dqG1_cFl(dFl,l1,l2,r,eps,veps,qd);
    tdqG[0]+=dFl[0];
    tdqG[1]+=dFl[1];
    tdqG[2]+=dFl[2];

    if( d3hm_sf_diff_z(&dFl[0],&tdqG[0],eps) )
      if( d3hm_sf_diff_z(&dFl[1],&tdqG[1],eps) ) break;
  }
  if(l2==lc+EW_LIMIT){
    return -1;
  }

  // l2 < 0
  for(l2=lc-1;l2>lc-EW_LIMIT;l2--){
    cc+=1;
    ew_dqG1_cFl(dFl,l1,l2,r,eps,veps,qd);
    tdqG[0]+=dFl[0];
    tdqG[1]+=dFl[1];
    tdqG[2]+=dFl[2];

    if( d3hm_sf_diff_z(&dFl[0],&tdqG[0],eps) )
      if( d3hm_sf_diff_z(&dFl[1],&tdqG[1],eps) ) break;
  }
  if(l2==lc-EW_LIMIT){
    return -1;
  }

  return cc;
}

void ew_dqG1_cFl(double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd)
{
  double complex w,ce,cf;
  double rl,i_rl,Rl[2],ke,i_ve,rl_2e,arg,ef;

  Rl[0]=(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0];
  Rl[1]=(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1];
  rl=sqrt((r[0]+Rl[0])*(r[0]+Rl[0])+(r[1]+Rl[1])*(r[1]+Rl[1])+r[2]*r[2]);
  i_rl=1.0/rl;

  ke=qd->k*veps;
  i_ve=1.0/veps;
  rl_2e=rl*0.5*i_ve;

  arg=qd->vk[0]*Rl[0]+qd->vk[1]*Rl[1];
  ce=cos(arg)+I*sin(arg);
  ef=exp(ke*ke-rl_2e*rl_2e);
  w=Faddeeva_w(ke+I*rl_2e,eps);
  cf=ce*i_rl*i_rl*ef*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));

  dFl[0]=(r[0]+Rl[0])*cf;
  dFl[1]=(r[1]+Rl[1])*cf;
  dFl[2]=cf;
}

int ew_dqG2(double complex *dqG2,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_dqG2_sum_n2(double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err);

  double complex odqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,cc,err,nc;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0
  idd=1.0/detd;

  nc=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  cc=0;
  dqG2[0]=0.0;
  dqG2[1]=0.0;
  dqG2[2]=0.0;

  n1=nc;
  cc+=ew_dqG2_sum_n2(dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+EW_LIMIT;n1++){
    odqG=dqG2[2];
    cc+=ew_dqG2_sum_n2(dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&dqG2[2],&odqG,eps) ) break;
    }
  }
  if(n1==nc+EW_LIMIT){
    return -2;
  }

  // n1 < nc
  for(n1=nc-1;n1>nc-EW_LIMIT;n1--){
    odqG=dqG2[2];
    cc+=ew_dqG2_sum_n2(dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&dqG2[2],&odqG,eps) ) break;
    }
  }
  if(n1==nc-EW_LIMIT){
    return -2;
  }

  dqG2[0]*=-0.25*fabs(idd);
  dqG2[1]*=-0.25*fabs(idd);
  dqG2[2]*= 0.25*fabs(idd);

  return cc;
}

int ew_dqG2_sum_n2(double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err)
{
  int ew_dqG2_cFn(double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd);

  double complex dFn[3];
  int cc,n2,nc,flg;

  nc=(int)round(-((-qd->vk[0]+(double)n1*vg1[0])*vg2[0]+(-qd->vk[1]+(double)n1*vg1[1])*vg2[1])/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  flg=0;
  cc=0;
  n2=nc;
  cc+=1;
  *err=ew_dqG2_cFn(dFn,n1,n2,r,vg1,vg2,eps,veps,qd);  if(*err<0) return 0;
  tdqG[0]+=dFn[0];
  tdqG[1]+=dFn[1];
  tdqG[2]+=dFn[2];

  // n2 > nc
  for(n2=nc+1;n2<nc+EW_LIMIT;n2++){
    cc+=1;
    *err=ew_dqG2_cFn(dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc+EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-EW_LIMIT;n2--){
    cc+=1;
    *err=ew_dqG2_cFn(dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc-EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  if(flg==0) *err=1;
  else *err=0;

  return cc;
}

int ew_dqG2_cFn(double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd)
{
  double complex ce,w,ff,Fn;
  double an[2],bn2,bn,arg,azve,bnaz,ebz,ef,f1,f2;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];
  arg=an[0]*r[0]+an[1]*r[1];
  ce=cos(arg)+I*sin(arg);

  if(bn2>0.0){
    bn=sqrt(bn2);
    azve=r[2]/(2.0*veps);
    arg=bn*r[2];
    ef=exp(bn2*veps*veps-azve*azve);
    ff=cos(arg)+I*sin(arg);
    w=Faddeeva_w(bn*veps+I*azve,eps);

    Fn=2.0*ce/bn*(I*ef*cimag(w)+ff);
    dFn[0]=Fn*an[0];
    dFn[1]=Fn*an[1];
    dFn[2]=2.0*ce*(ef*creal(w)-ff);

    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2); // |bn|
    azve=r[2]/(2.0*veps);
    bnaz=bn*r[2];
    ebz=exp(bnaz);
    f1=ebz*erfc(bn*veps+azve);
    f2=1.0/ebz*erfc(bn*veps-azve);

    Fn=ce/(I*bn)*(f1+f2);
    dFn[0]=Fn*an[0];
    dFn[1]=Fn*an[1];
    dFn[2]=ce*(f1-f2);

    return 1;
  }
  else { // bn=0.0;
    return -3;
  }
}

int d3hm_qpgf_d2v1_ew_sd(double complex *qG,double complex *dqG,double *r,double eps,QPD2 *qd)
{
  int ew_qG1_v1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPD2 *qd);
  int ew_qG2_v1(double complex *qG2,double complex *dqG2,double *r,double eps,double veps,QPD2 *qd);

  double complex qG1,dqG1[3],qG2,dqG2[3];
  double detd,veps,adp,adm,r2_min;
  int c1,c2,l1,l2;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0

  adp=sqrt(pow(qd->vd1[0]+qd->vd2[0],2)+pow(qd->vd1[1]+qd->vd2[1],2));
  adm=sqrt(pow(qd->vd1[0]-qd->vd2[0],2)+pow(qd->vd1[1]-qd->vd2[1],2));
  veps=sqrt(fabs(detd)/(4.0*M_PI)*adp/adm);

  l1=0;
  l2=0;
  r2_min=pow(r[0]+(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0],2)
      +  pow(r[1]+(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1],2)+r[2]*r[2];

  if( qd->k*qd->k*veps*veps-r2_min/(4.0*veps*veps)>EW_EXPM )
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r2_min+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  c1=ew_qG1_v1(&qG1,dqG1,r,eps,veps,qd);
  *qG=qG1;
  dqG[0]=dqG1[0];
  dqG[1]=dqG1[1];
  dqG[2]=dqG1[2];
  if(c1<0) return c1;

  // qG2
  c2=ew_qG2_v1(&qG2,dqG2,r,eps,veps,qd);
  *qG+=qG2;
  dqG[0]+=dqG2[0];
  dqG[1]+=dqG2[1];
  dqG[2]+=dqG2[2];
  if(c2<0) return c2;

  return c1+c2;
}

int ew_qG1_v1(double complex *qG1,double complex *dqG1,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG1_v1_sum_l2(double complex *tqG,double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd);

  double complex oqG;
  int l1,cc,lc;

  lc=0;

  cc=0;
  *qG1=0.0;
  dqG1[0]=0.0;
  dqG1[1]=0.0;
  dqG1[2]=0.0;

  l1=lc;
  cc+=ew_qG1_v1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd);  if(cc<0) return cc;

  // l1 > lc
  for(l1=lc+1;l1<lc+EW_LIMIT;l1++){
    oqG=*qG1;
    cc+=ew_qG1_v1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(qG1,&oqG,eps) ) break;
  }
  if(l1==lc+EW_LIMIT){
    return -1;
  }

  // l1 < lc
  for(l1=lc-1;l1>lc-EW_LIMIT;l1--){
    oqG=*qG1;
    cc+=ew_qG1_v1_sum_l2(qG1,dqG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(qG1,&oqG,eps) ) break;
  }
  if(l1==lc-EW_LIMIT){
    return -1;
  }

  *qG1/=4.0*M_PI;
  dqG1[0]/=-4.0*M_PI;
  dqG1[1]/=-4.0*M_PI;
  dqG1[2]*=-r[2]/(4.0*M_PI);

  return cc;
}

int ew_qG1_v1_sum_l2(double complex *tqG,double complex *tdqG,int l1,double *r,double eps,double veps,QPD2 *qd)
{
  void ew_qG1_v1_cFl(double complex *Fl,double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd);

  double complex Fl,dFl[3];
  int cc,l2,lc;

  lc=-(int)round( ((r[0]+(double)l1*qd->vd1[0])*qd->vd2[0]+(r[1]+(double)l1*qd->vd1[1])*qd->vd2[1])
                    /(qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1]) );

  cc=0;
  l2=lc;
  cc+=1;
  ew_qG1_v1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
  *tqG+=Fl;
  tdqG[0]+=dFl[0];
  tdqG[1]+=dFl[1];
  tdqG[2]+=dFl[2];

  // l2 > lc
  for(l2=lc+1;l2<lc+EW_LIMIT;l2++){
    cc+=1;
    ew_qG1_v1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;
    tdqG[0]+=dFl[0];
    tdqG[1]+=dFl[1];
    tdqG[2]+=dFl[2];

    if( d3hm_sf_diff_z(&Fl,tqG,eps) ) break;
  }
  if(l2==lc+EW_LIMIT){
    return -1;
  }

  // l2 < 0
  for(l2=lc-1;l2>lc-EW_LIMIT;l2--){
    cc+=1;
    ew_qG1_v1_cFl(&Fl,dFl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;
    tdqG[0]+=dFl[0];
    tdqG[1]+=dFl[1];
    tdqG[2]+=dFl[2];

    if( d3hm_sf_diff_z(&Fl,tqG,eps) ) break;
  }
  if(l2==lc-EW_LIMIT){
    return -1;
  }

  return cc;
}

void ew_qG1_v1_cFl(double complex *Fl,double complex *dFl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd)
{
  double complex w,ce,cf;
  double rl,i_rl,Rl[2],ke,i_ve,rl_2e,arg,ef;

  Rl[0]=(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0];
  Rl[1]=(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1];
  rl=sqrt((r[0]+Rl[0])*(r[0]+Rl[0])+(r[1]+Rl[1])*(r[1]+Rl[1])+r[2]*r[2]);
  i_rl=1.0/rl;

  ke=qd->k*veps;
  i_ve=1.0/veps;
  rl_2e=rl*0.5*i_ve;

  arg=qd->vk[0]*Rl[0]+qd->vk[1]*Rl[1];
  ce=cos(arg)+I*sin(arg);
  ef=exp(ke*ke-rl_2e*rl_2e);
  w=Faddeeva_w(ke+I*rl_2e,eps);
  cf=ce*i_rl*i_rl*ef*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));

  *Fl=ce*i_rl*ef*creal(w);
  dFl[0]=(r[0]+Rl[0])*cf;
  dFl[1]=(r[1]+Rl[1])*cf;
  dFl[2]=cf;
}

int ew_qG2_v1(double complex *qG2,double complex *dqG2,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG2_v1_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err);

  double complex odqG;
  double vg1[2],vg2[2],detd,idd;
  int n1,cc,err,nc;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0
  idd=1.0/detd;

  nc=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  cc=0;
  *qG2=0.0;
  dqG2[0]=0.0;
  dqG2[1]=0.0;
  dqG2[2]=0.0;

  n1=nc;
  cc+=ew_qG2_v1_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+EW_LIMIT;n1++){
    odqG=dqG2[2];
    cc+=ew_qG2_v1_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&dqG2[2],&odqG,eps) ) break;
    }
  }
  if(n1==nc+EW_LIMIT){
    return -2;
  }

  // n1 < nc
  for(n1=nc-1;n1>nc-EW_LIMIT;n1--){
    odqG=dqG2[2];
    cc+=ew_qG2_v1_sum_n2(qG2,dqG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&dqG2[2],&odqG,eps) ) break;
    }
  }
  if(n1==nc-EW_LIMIT){
    return -2;
  }

  *qG2*=0.25*I*fabs(idd);
  dqG2[0]*=-0.25*fabs(idd);
  dqG2[1]*=-0.25*fabs(idd);
  dqG2[2]*= 0.25*fabs(idd);

  return cc;
}

int ew_qG2_v1_sum_n2(double complex *tqG,double complex *tdqG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err)
{
  int ew_qG2_v1_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd);

  double complex Fn,dFn[3];
  int cc,n2,nc,flg;

  nc=(int)round(-((-qd->vk[0]+(double)n1*vg1[0])*vg2[0]+(-qd->vk[1]+(double)n1*vg1[1])*vg2[1])/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  flg=0;
  cc=0;
  n2=nc;
  cc+=1;
  *err=ew_qG2_v1_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd);  if(*err<0) return 0;
  *tqG+=Fn;
  tdqG[0]+=dFn[0];
  tdqG[1]+=dFn[1];
  tdqG[2]+=dFn[2];

  // n2 > nc
  for(n2=nc+1;n2<nc+EW_LIMIT;n2++){
    cc+=1;
    *err=ew_qG2_v1_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc+EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-EW_LIMIT;n2--){
    cc+=1;
    *err=ew_qG2_v1_cFn(&Fn,dFn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;
    tdqG[0]+=dFn[0];
    tdqG[1]+=dFn[1];
    tdqG[2]+=dFn[2];

    if(*err>0){
      if( d3hm_sf_diff_z(&dFn[2],&tdqG[2],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc-EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  if(flg==0) *err=1;
  else *err=0;

  return cc;
}

int ew_qG2_v1_cFn(double complex *Fn,double complex *dFn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd)
{
  double complex ce,w,ff;
  double an[2],bn2,bn,arg,azve,bnaz,ebz,ef,f1,f2;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];
  arg=an[0]*r[0]+an[1]*r[1];
  ce=cos(arg)+I*sin(arg);

  if(bn2>0.0){
    bn=sqrt(bn2);
    azve=r[2]/(2.0*veps);
    arg=bn*r[2];
    ef=exp(bn2*veps*veps-azve*azve);
    ff=cos(arg)+I*sin(arg);
    w=Faddeeva_w(bn*veps+I*azve,eps);

    *Fn=2.0*ce/bn*(I*ef*cimag(w)+ff);
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    dFn[2]=2.0*ce*(ef*creal(w)-ff);

    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2); // |bn|
    azve=r[2]/(2.0*veps);
    bnaz=bn*r[2];
    ebz=exp(bnaz);
    f1=ebz*erfc(bn*veps+azve);
    f2=1.0/ebz*erfc(bn*veps-azve);

    *Fn=ce/(I*bn)*(f1+f2);
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    dFn[2]=ce*(f1-f2);

    return 1;
  }
  else { // bn=0.0;
    return -3;
  }
}

int d3hm_qpgf_d2v2_ew_sd(double complex *qG,double complex *dqG,double complex *d2qG,double *r,double eps,QPD2 *qd)
{
  int ew_qG1_v2(double complex *qG1,double complex *dqG1,double complex *d2qG1,double *r,double eps,double veps,QPD2 *qd);
  int ew_qG2_v2(double complex *qG2,double complex *dqG2,double complex *d2qG2,double *r,double eps,double veps,QPD2 *qd);

  double complex qG1,dqG1[3],d2qG1[6],qG2,dqG2[3],d2qG2[6];
  double detd,veps,adp,adm,r2_min;
  int c1,c2,l1,l2,i;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0

  adp=sqrt(pow(qd->vd1[0]+qd->vd2[0],2)+pow(qd->vd1[1]+qd->vd2[1],2));
  adm=sqrt(pow(qd->vd1[0]-qd->vd2[0],2)+pow(qd->vd1[1]-qd->vd2[1],2));
  veps=sqrt(fabs(detd)/(4.0*M_PI)*adp/adm);

  l1=0;
  l2=0;
  r2_min=pow(r[0]+(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0],2)
      +  pow(r[1]+(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1],2)+r[2]*r[2];

  if( qd->k*qd->k*veps*veps-r2_min/(4.0*veps*veps)>EW_EXPM )
    veps=sqrt(EW_EXPM+sqrt(qd->k*qd->k*r2_min+EW_EXPM*EW_EXPM))/(sqrt(2.0)*qd->k);

  // qG1
  c1=ew_qG1_v2(&qG1,dqG1,d2qG1,r,eps,veps,qd);
  *qG=qG1;
  for(i=0;i<3;i++) dqG[i]=dqG1[i];
  for(i=0;i<6;i++) d2qG[i]=d2qG1[i];
  if(c1<0) return c1;

  // qG2
  c2=ew_qG2_v2(&qG2,dqG2,d2qG2,r,eps,veps,qd);
  *qG+=qG2;
  for(i=0;i<3;i++) dqG[i]+=dqG2[i];
  for(i=0;i<6;i++) d2qG[i]+=d2qG2[i];
  if(c2<0) return c2;

  return c1+c2;
}

int ew_qG1_v2(double complex *qG1,double complex *dqG1,double complex *d2qG1,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG1_v2_sum_l2(double complex *tqG,double complex *tdqG,double complex *td2qG,int l1,double *r,double eps,double veps,QPD2 *qd);

  double complex od2qG[6];
  int l1,cc,lc,i;

  lc=0;

  cc=0;
  *qG1=0.0;
  for(i=0;i<3;i++) dqG1[i]=0.0;
  for(i=0;i<6;i++) d2qG1[i]=0.0;

  l1=lc;
  cc+=ew_qG1_v2_sum_l2(qG1,dqG1,d2qG1,l1,r,eps,veps,qd);  if(cc<0) return cc;

  // l1 > lc
  for(l1=lc+1;l1<lc+EW_LIMIT;l1++){
    for(i=0;i<2;i++) od2qG[i]=d2qG1[i];
    cc+=ew_qG1_v2_sum_l2(qG1,dqG1,d2qG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(&d2qG1[0],&od2qG[0],eps) )
      if( d3hm_sf_comp_z(&d2qG1[1],&od2qG[1],eps) ) break;
  }
  if(l1==lc+EW_LIMIT){
    return -1;
  }

  // l1 < lc
  for(l1=lc-1;l1>lc-EW_LIMIT;l1--){
    for(i=0;i<2;i++) od2qG[i]=d2qG1[i];
    cc+=ew_qG1_v2_sum_l2(qG1,dqG1,d2qG1,l1,r,eps,veps,qd); if(cc<0) return cc;

    if( d3hm_sf_comp_z(&d2qG1[0],&od2qG[0],eps) )
      if( d3hm_sf_comp_z(&d2qG1[1],&od2qG[1],eps) ) break;
  }
  if(l1==lc-EW_LIMIT){
    return -1;
  }

  *qG1/=4.0*M_PI;
  dqG1[0]/=-4.0*M_PI;
  dqG1[1]/=-4.0*M_PI;
  dqG1[2]*=-r[2]/(4.0*M_PI);
  d2qG1[0]/=-8.0*M_PI;
  d2qG1[1]/=-8.0*M_PI;
  d2qG1[2]/=-8.0*M_PI;
  d2qG1[3]/=-8.0*M_PI;
  d2qG1[4]*=-r[2]/(8.0*M_PI);
  d2qG1[5]*=-r[2]/(8.0*M_PI);

  return cc;
}

int ew_qG1_v2_sum_l2(double complex *tqG,double complex *tdqG,double complex *td2qG,int l1,double *r,double eps,double veps,QPD2 *qd)
{
  void ew_qG1_v2_cFl(double complex *Fl,double complex *dFl,double complex *d2Fl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd);

  double complex Fl,dFl[3],d2Fl[6];
  int cc,l2,lc,i;

  lc=-(int)round( ((r[0]+(double)l1*qd->vd1[0])*qd->vd2[0]+(r[1]+(double)l1*qd->vd1[1])*qd->vd2[1])
                    /(qd->vd2[0]*qd->vd2[0]+qd->vd2[1]*qd->vd2[1]) );

  cc=0;
  l2=lc;
  cc+=1;
  ew_qG1_v2_cFl(&Fl,dFl,d2Fl,l1,l2,r,eps,veps,qd);
  *tqG+=Fl;
  for(i=0;i<3;i++) tdqG[i]+=dFl[i];
  for(i=0;i<6;i++) td2qG[i]+=d2Fl[i];

  // l2 > lc
  for(l2=lc+1;l2<lc+EW_LIMIT;l2++){
    cc+=1;
    ew_qG1_v2_cFl(&Fl,dFl,d2Fl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;
    for(i=0;i<3;i++) tdqG[i]+=dFl[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fl[i];

    if( d3hm_sf_diff_z(&d2Fl[0],&td2qG[0],eps) )
      if( d3hm_sf_diff_z(&d2Fl[1],&td2qG[1],eps) ) break;
  }
  if(l2==lc+EW_LIMIT){
    return -1;
  }

  // l2 < 0
  for(l2=lc-1;l2>lc-EW_LIMIT;l2--){
    cc+=1;
    ew_qG1_v2_cFl(&Fl,dFl,d2Fl,l1,l2,r,eps,veps,qd);
    *tqG+=Fl;
    for(i=0;i<3;i++) tdqG[i]+=dFl[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fl[i];

    if( d3hm_sf_diff_z(&d2Fl[0],&td2qG[0],eps) )
      if( d3hm_sf_diff_z(&d2Fl[1],&td2qG[1],eps) ) break;
  }
  if(l2==lc-EW_LIMIT){
    return -1;
  }

  return cc;
}

void ew_qG1_v2_cFl(double complex *Fl,double complex *dFl,double complex *d2Fl,int l1,int l2,double *r,double eps,double veps,QPD2 *qd)
{
  double complex w,ce,cf1,cf2,cf3,cef;
  double rl,rl2,i_rl,i_rl2,Rl[2],ke,i_ve,rl_2e,arg,ef,ra[3];

  Rl[0]=(double)l1*qd->vd1[0]+(double)l2*qd->vd2[0];
  Rl[1]=(double)l1*qd->vd1[1]+(double)l2*qd->vd2[1];
  rl2=(r[0]+Rl[0])*(r[0]+Rl[0])+(r[1]+Rl[1])*(r[1]+Rl[1])+r[2]*r[2];
  rl=sqrt(rl2);
  i_rl=1.0/rl;
  i_rl2=i_rl*i_rl;
  ra[0]=r[0]+Rl[0];
  ra[1]=r[1]+Rl[1];
  ra[2]=r[2];

  ke=qd->k*veps;
  i_ve=1.0/veps;
  rl_2e=rl*0.5*i_ve;

  arg=qd->vk[0]*Rl[0]+qd->vk[1]*Rl[1];
  ce=cos(arg)+I*sin(arg);
  ef=exp(ke*ke-rl_2e*rl_2e);
  cef=ce*ef;
  w=Faddeeva_w(ke+I*rl_2e,eps);

  cf1=cef*i_rl2*(i_rl*creal(w)-qd->k*cimag(w)+i_ve/sqrt(M_PI));
  cf2=cef*i_rl2*(qd->k*qd->k*i_rl*2.0*creal(w)-i_ve*i_ve*i_ve/sqrt(M_PI));
  cf3=-3.0*i_rl2*2.0*cf1+cf2;

  *Fl=cef*i_rl*creal(w);

  dFl[0]=ra[0]*cf1;
  dFl[1]=ra[1]*cf1;
  dFl[2]=cf1;

  d2Fl[0]=2.0*cf1+ra[0]*ra[0]*cf3;
  d2Fl[1]=2.0*cf1+ra[1]*ra[1]*cf3;
  d2Fl[2]=2.0*cf1+ra[2]*ra[2]*cf3;
  d2Fl[3]=ra[0]*ra[1]*cf3;
  d2Fl[4]=ra[1]*cf3;
  d2Fl[5]=ra[0]*cf3;
}

int ew_qG2_v2(double complex *qG2,double complex *dqG2,double complex *d2qG2,double *r,double eps,double veps,QPD2 *qd)
{
  int ew_qG2_v2_sum_n2(double complex *tqG,double complex *tdqG,double complex *td2qG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err);

  double complex od2qG[6],ci;
  double vg1[2],vg2[2],detd,idd,cd;
  int n1,cc,err,nc,i;

  detd=qd->vd1[0]*qd->vd2[1]-qd->vd2[0]*qd->vd1[1];
  if(detd==0.0) return -4; // determinant d == 0
  idd=1.0/detd;

  nc=(int)round( (qd->vk[0]*qd->vd1[0]+qd->vk[1]*qd->vd1[1])/(2.0*M_PI) );

  vg1[0]= 2.0*M_PI*idd*qd->vd2[1];
  vg1[1]=-2.0*M_PI*idd*qd->vd2[0];
  vg2[0]=-2.0*M_PI*idd*qd->vd1[1];
  vg2[1]= 2.0*M_PI*idd*qd->vd1[0];

  cc=0;
  *qG2=0.0;
  for(i=0;i<3;i++) dqG2[i]=0.0;
  for(i=0;i<6;i++) d2qG2[i]=0.0;

  n1=nc;
  cc+=ew_qG2_v2_sum_n2(qG2,dqG2,d2qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

  // n1 > nc
  for(n1=nc+1;n1<nc+EW_LIMIT;n1++){
    for(i=0;i<2;i++) od2qG[i]=d2qG2[i];
    cc+=ew_qG2_v2_sum_n2(qG2,dqG2,d2qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&d2qG2[0],&od2qG[0],eps) )
        if( d3hm_sf_comp_z(&d2qG2[1],&od2qG[1],eps) ) break;
    }
  }
  if(n1==nc+EW_LIMIT){
    return -2;
  }

  // n1 < nc
  for(n1=nc-1;n1>nc-EW_LIMIT;n1--){
    for(i=0;i<2;i++) od2qG[i]=d2qG2[i];
    cc+=ew_qG2_v2_sum_n2(qG2,dqG2,d2qG2,n1,r,vg1,vg2,eps,veps,qd,&err); if(err<0) return err;

    if(err>0){
      if( d3hm_sf_comp_z(&d2qG2[0],&od2qG[0],eps) )
        if( d3hm_sf_comp_z(&d2qG2[1],&od2qG[1],eps) ) break;
    }
  }
  if(n1==nc-EW_LIMIT){
    return -2;
  }

  cd=0.25*fabs(idd);
  ci=cd*I;

  *qG2*=ci;
  dqG2[0]*=-cd;
  dqG2[1]*=-cd;
  dqG2[2]*= cd;
  d2qG2[0]*=-ci;
  d2qG2[1]*=-ci;
  d2qG2[2]*=-cd;
  d2qG2[3]*=-ci;
  d2qG2[4]*=ci;
  d2qG2[5]*=ci;

  return cc;
}

int ew_qG2_v2_sum_n2(double complex *tqG,double complex *tdqG,double complex *td2qG,int n1,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd,int *err)
{
  int ew_qG2_v2_cFn(double complex *Fn,double complex *dFn,double complex *d2Fn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd);

  double complex Fn,dFn[3],d2Fn[6];
  int cc,n2,nc,i,flg;

  nc=(int)round(-((-qd->vk[0]+(double)n1*vg1[0])*vg2[0]+(-qd->vk[1]+(double)n1*vg1[1])*vg2[1])/(vg2[0]*vg2[0]+vg2[1]*vg2[1]));

  flg=0;
  cc=0;
  n2=nc;
  cc+=1;
  *err=ew_qG2_v2_cFn(&Fn,dFn,d2Fn,n1,n2,r,vg1,vg2,eps,veps,qd);  if(*err<0) return 0;
  *tqG+=Fn;
  for(i=0;i<3;i++) tdqG[i]+=dFn[i];
  for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];

  // n2 > nc
  for(n2=nc+1;n2<nc+EW_LIMIT;n2++){
    cc+=1;
    *err=ew_qG2_v2_cFn(&Fn,dFn,d2Fn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;
    for(i=0;i<3;i++) tdqG[i]+=dFn[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];

    if(*err>0){
      if( d3hm_sf_diff_z(&d2Fn[0],&td2qG[0],eps) )
        if( d3hm_sf_diff_z(&d2Fn[1],&td2qG[1],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc+EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  // n2 < nc
  for(n2=nc-1;n2>nc-EW_LIMIT;n2--){
    cc+=1;
    *err=ew_qG2_v2_cFn(&Fn,dFn,d2Fn,n1,n2,r,vg1,vg2,eps,veps,qd); if(*err<0) return 0;
    *tqG+=Fn;
    for(i=0;i<3;i++) tdqG[i]+=dFn[i];
    for(i=0;i<6;i++) td2qG[i]+=d2Fn[i];

    if(*err>0){
      if( d3hm_sf_diff_z(&d2Fn[0],&td2qG[0],eps) )
        if( d3hm_sf_diff_z(&d2Fn[1],&td2qG[1],eps) ) break;
    }
    else flg++;
  }
  if(n2==nc-EW_LIMIT){
    *err=-2;
    return 0; // summation limit
  }

  if(flg==0) *err=1;
  else *err=0;

  return cc;
}

int ew_qG2_v2_cFn(double complex *Fn,double complex *dFn,double complex *d2Fn,int n1,int n2,double *r,double *vg1,double *vg2,double eps,double veps,QPD2 *qd)
{
  double complex ce,w,ff,F,G;
  double an[2],bn2,bn,arg,azve,bnaz,ebz,ef,f1,f2;

  an[0]=-qd->vk[0]+(double)n1*vg1[0]+(double)n2*vg2[0];
  an[1]=-qd->vk[1]+(double)n1*vg1[1]+(double)n2*vg2[1];
  bn2=qd->k*qd->k-an[0]*an[0]-an[1]*an[1];
  arg=an[0]*r[0]+an[1]*r[1];
  ce=cos(arg)+I*sin(arg);
  azve=r[2]/(2.0*veps);
  ef=exp(bn2*veps*veps-azve*azve);

  if(bn2>0.0){
    bn=sqrt(bn2);
    arg=bn*r[2];
    ff=cos(arg)+I*sin(arg);
    w=Faddeeva_w(bn*veps+I*azve,eps);

    F=2.0*(I*ef*cimag(w)+ff);
    G=2.0*(ef*creal(w)-ff);

    *Fn=ce/bn*F;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    dFn[2]=ce*G;

    d2Fn[0]=dFn[0]*an[0];
    d2Fn[1]=dFn[1]*an[1];
    d2Fn[2]=ce*(I*bn*F+2.0/(sqrt(M_PI)*veps)*ef);
    d2Fn[3]=dFn[0]*an[1];
    d2Fn[4]=an[1]*dFn[2];
    d2Fn[5]=an[0]*dFn[2];

    return 0;
  }
  else if(bn2<0.0){
    bn=sqrt(-bn2); // |bn|
    bnaz=bn*r[2];
    ebz=exp(bnaz);
    f1=ebz*erfc(bn*veps+azve);
    f2=1.0/ebz*erfc(bn*veps-azve);

    F=f1+f2;
    G=f1-f2;

    *Fn=-I*ce/bn*F;
    dFn[0]=*Fn*an[0];
    dFn[1]=*Fn*an[1];
    dFn[2]=ce*G;

    d2Fn[0]=dFn[0]*an[0];
    d2Fn[1]=dFn[1]*an[1];
    d2Fn[2]=ce*(-bn*F+2.0/(sqrt(M_PI)*veps)*ef);
    d2Fn[3]=dFn[0]*an[1];
    d2Fn[4]=an[1]*dFn[2];
    d2Fn[5]=an[0]*dFn[2];

    return 1;
  }
  else { // bn=0.0;
    return -3;
  }
}

