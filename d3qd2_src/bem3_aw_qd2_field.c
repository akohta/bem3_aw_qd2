#include "bem3_aw_qd2.h"

int phi_s_aqd2(double complex *phi,double *rt,int type,AQD2 *ad)
{
  int domain_id(double *rt,AQD2 *qd);// bem3_aw_qd2.c, non-periodicity
  
  double complex CC[9],kc,ce;
  double F,r0[3],arg;
  int did,s,sd,n,l1,l2;

  if(ad->bd.ps==1){
    domain_id_l_aqd2(&l1,&l2,rt,ad);
    r0[0]=rt[0]-(double)l1*ad->bd.qd.vd1[0]-(double)l2*ad->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*ad->bd.qd.vd1[1]-(double)l2*ad->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }

  *phi=0.0;
  F=0.0;
  kc=(double complex)ad->k0[did];

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];
    if(ad->bd.ps==1 && did==0) q0_coef_rt(CC,r0,sd,kc,type,&(ad->bd));
    else coef_rt(CC,r0,sd,kc,type,&(ad->bd));

    for(n=0;n<4;n++) *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];

    F+=creal(CC[8]);
  }

  if(did==0) *phi/=1.0+F;
  else *phi/=F;
  
  if(l1!=0 || l2!=0){
    arg= ad->bd.qd.vk[0]*((double)l1*ad->bd.qd.vd1[0]+(double)l2*ad->bd.qd.vd2[0])
        +ad->bd.qd.vk[1]*((double)l1*ad->bd.qd.vd1[1]+(double)l2*ad->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    *phi*=ce;
  }
  
  return did;
}

int phi_t_aqd2(double complex *phi,double *rt,int type,AQD2 *ad)
{
  double complex p,v[3];
  int did;

  did=phi_s_aqd2(phi,rt,type,ad);
  if(did==0){
    calc_apw_pv(&p,v,rt,&(ad->pw)); // aw_pw.h
    *phi+=-p/ad->pw.wd.k2;
  }

  return did;
}

int phi_i_aqd2(double complex *phi,double *rt,int type,AQD2 *ad)
{
  double complex p,v[3];

  calc_apw_pv(&p,v,rt,&(ad->pw)); // aw_pw.h
  *phi=-p/ad->pw.wd.k2;
  return domain_id_aqd2(rt,ad);
}

int p_s_aqd2(double complex *p,double *rt,int type,AQD2 *ad)
{
  int did;

  did=phi_s_aqd2(p,rt,type,ad);
  *p*=-ad->k2[did];
  return did;
}

int p_t_aqd2(double complex *p,double *rt,int type,AQD2 *ad)
{
  int did;

  did=phi_t_aqd2(p,rt,type,ad);
  *p*=-ad->k2[did];
  return did;
}

int p_i_aqd2(double complex *p,double *rt,int type,AQD2 *ad)
{
  double complex v[3];

  calc_apw_pv(p,v,rt,&(ad->pw)); // aw_pw.h 
  return domain_id_aqd2(rt,ad);
}

int pv_s_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  int pv_s_bieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad);
  int pv_s_dbieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad);
  
  int did;

  did=pv_s_dbieq_aqd2(p,pv,rt,type,ad);
  if(did<0){ // near boundary
    did=pv_s_bieq_aqd2(p,pv,rt,type,ad);
  }

  return did; 
}

int pv_t_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  double complex pvi[3],pi;
  int did,l;

  did=pv_s_aqd2(p,pv,rt,type,ad);

  if(did==0){
    calc_apw_pv(&pi,pvi,rt,&(ad->pw));
    *p+=pi;
    for(l=0;l<3;l++) pv[l]+=pvi[l];
  }

  return did;
}

int pv_i_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
   calc_apw_pv(p,pv,rt,&(ad->pw));
  return domain_id_aqd2(rt,ad);
}

int pv_s_bieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  int domain_id(double *rt,AQD2 *qd);// bem3_aw_qd2.c, non-periodicity
  
  double complex CC[9],ce;
  double F,cf,r0[3],arg;
  int did,s,sd,l1,l2,n,i;

  if(ad->bd.ps==1){
    domain_id_l_aqd2(&l1,&l2,rt,ad);
    r0[0]=rt[0]-(double)l1*ad->bd.qd.vd1[0]-(double)l2*ad->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*ad->bd.qd.vd1[1]-(double)l2*ad->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }

  *p=0.0;
  for(i=0;i<3;i++) pv[i]=0.0;
  F=0.0;
  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];
    
    if(did==0 && ad->bd.ps==1) q0_coef_rt(CC,r0,sd,ad->k0[did],type,&(ad->bd));
    else coef_rt(CC,r0,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(i=0;i<3;i++) pv[i]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][i]-CC[n+4]*ad->bd.sb[did].pv[s][n][i];
    }

    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    *p*=cf;
    for(i=0;i<3;i++) pv[i]*=cf;
  }
  else{
    cf=1.0/F;
    *p*=cf;
    for(i=0;i<3;i++) pv[i]*=cf;
  }
  *p*=-ad->k2[did];
  
  if(l1!=0 || l2!=0){
    arg= ad->bd.qd.vk[0]*((double)l1*ad->bd.qd.vd1[0]+(double)l2*ad->bd.qd.vd2[0])
        +ad->bd.qd.vk[1]*((double)l1*ad->bd.qd.vd1[1]+(double)l2*ad->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    *p*=ce;
    for(i=0;i<3;i++) pv[i]*=ce;
  }
  return did;
}

int pv_t_bieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  double complex pi,pvi[3];
  int did,l;

  did=pv_s_bieq_aqd2(p,pv,rt,type,ad);
  if(did==0){
    calc_apw_pv(&pi,pvi,rt,&(ad->pw)); // aw_pw.h
    *p+=pi;
    for(l=0;l<3;l++) pv[l]+=pvi[l];
  }

  return did;
}

int pv_i_bieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  int did;

  did=domain_id_aqd2(rt,ad);
  calc_apw_pv(p,pv,rt,&(ad->pw)); // aw_pw.h
  return did;
}

int pv_s_dbieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  int domain_id(double *rt,AQD2 *qd);// bem3_aw_qd2.c, non-periodicity
  
  double complex CC[9],dC[3][9],kc,P,dP[3],ce;
  double F,dF[3],r0[3],arg,i_C;
  int did,l1,l2,sd,s,i,n;

  if(ad->bd.ps==1){
    domain_id_l_aqd2(&l1,&l2,rt,ad);
    r0[0]=rt[0]-(double)l1*ad->bd.qd.vd1[0]-(double)l2*ad->bd.qd.vd2[0];
    r0[1]=rt[1]-(double)l1*ad->bd.qd.vd1[1]-(double)l2*ad->bd.qd.vd2[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }
  else {
    l1=0;
    l2=0;
    r0[0]=rt[0];
    r0[1]=rt[1];
    r0[2]=rt[2];
    did=domain_id(r0,ad);
  }

  kc=ad->k0[did];

  F=0.0;
  P=0.0;
  for(i=0;i<3;i++){
    dP[i]=0.0;
    dF[i]=0.0;
  }

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    if(did==0 && ad->bd.ps==1) q0_dcoef_rt_grad(CC,dC,r0,sd,kc,type,&(ad->bd));
    else dcoef_rt_grad(CC,dC,r0,sd,kc,type,&(ad->bd));

    for(i=0;i<3;i++){
      for(n=0;n<4;n++){
        if(i==0) P+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
        dP[i]+=dC[i][n+0]*ad->bd.sb[did].dP[s][n]-dC[i][n+4]*ad->bd.sb[did].P[s][n];
      }
      dF[i]+=creal(dC[i][8]);
    }
    F+=creal(CC[8]);
  }

  if(fabs(dF[0])>CBD_CDF || fabs(dF[1])>CBD_CDF || fabs(dF[2])>CBD_CDF){
    *p=0.0;
    for(i=0;i<3;i++) pv[i]=0.0;
    if(did==0) return -OPENDID;
    else return -did;
  }

  if(did==0) i_C=1.0/(1.0+F);
  else i_C=1.0/F;
  *p=-ad->k2[did]*P*i_C; // p=-k2*phi
  for(i=0;i<3;i++) pv[i]=-(dP[i]-P*dF[i])*i_C; // v=-nabla phi
  
  if(l1!=0 || l2!=0){
    arg= ad->bd.qd.vk[0]*((double)l1*ad->bd.qd.vd1[0]+(double)l2*ad->bd.qd.vd2[0])
        +ad->bd.qd.vk[1]*((double)l1*ad->bd.qd.vd1[1]+(double)l2*ad->bd.qd.vd2[1]);
    ce=cos(arg)+I*sin(arg);
    *p*=ce;
    for(i=0;i<3;i++) pv[i]*=ce;
  }

  return did;
}

int pv_t_dbieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  double complex pi,pvi[3];
  int did,l;

  did=pv_s_dbieq_aqd2(p,pv,rt,type,ad);
  if(did<0) return did;

  if(did==0){
    calc_apw_pv(&pi,pvi,rt,&(ad->pw)); // aw_pw.h
    *p+=pi;
    for(l=0;l<3;l++) pv[l]+=pvi[l];
  }

  return did;
}

int pv_i_dbieq_aqd2(double complex *p,double complex *pv,double *rt,int type,AQD2 *ad)
{
  return pv_i_aqd2(p,pv,rt,type,ad);
}

void phi_s_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  double complex CC[9];
  double F,rt[3];
  int s,sd,n,td;

  *phi=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    if(did==0 && ad->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++) *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
    F+=creal(CC[8]);
  }

  if(did==0) *phi/=1.0+F;
  else *phi/=F;
}

void phi_t_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  double complex CC[9],pi,v[3];
  double F,rt[3];
  int s,sd,n,td;

  *phi=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];
    
    if(did==0 && ad->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++) *phi+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
    F+=creal(CC[8]);
  }

  if(did==0){
    *phi/=1.0+F;
    calc_apw_pv(&pi,v,rt,&(ad->pw)); // aw_pw.h
    *phi+=-pi/ad->pw.wd.k2;
  }
  else{
    *phi/=F;
  }
}

void phi_i_bd_aqd2(double complex *phi,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  double complex v[3];
  double rt[3];
  int td;

  if(did==0){
    td=ad->bd.sb[did].sid[t];
    r_bd(rt,td,zeta_t,eta_t,&(ad->bd));
    calc_apw_pv(phi,v,rt,&(ad->pw)); // aw_pw.h
    *phi/=-ad->pw.wd.k2;
  }
  else *phi=0.0;
}

void pv_s_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  double complex CC[9];
  double F,rt[3],cf;
  int s,sd,l,n,td;

  *p=0.0;
  for(l=0;l<3;l++) pv[l]=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    if(did==0 && ad->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(l=0;l<3;l++) pv[l]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][l]-CC[n+4]*ad->bd.sb[did].pv[s][n][l];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  }
  else{
    cf=1.0/F;
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  } 
}

void pv_t_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  double complex CC[9],pi,vi[3];
  double F,rt[3],cf;
  int s,sd,l,n,td;

  *p=0.0;
  for(l=0;l<3;l++) pv[l]=0.0;
  F=0.0;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  for(s=1;s<=ad->bd.sb[did].Ne;s++){
    sd=ad->bd.sb[did].sid[s];

    if(did==0 && ad->bd.ps==1) q0_coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));
    else coef_bd(CC,rt,td,zeta_t,eta_t,sd,ad->k0[did],type,&(ad->bd));

    for(n=0;n<4;n++){
      *p+=CC[n+0]*ad->bd.sb[did].dP[s][n]-CC[n+4]*ad->bd.sb[did].P[s][n];
      for(l=0;l<3;l++) pv[l]+=CC[n+0]*ad->bd.sb[did].dpv[s][n][l]-CC[n+4]*ad->bd.sb[did].pv[s][n][l];
    }
    F+=creal(CC[8]);
  }

  if(did==0){
    cf=1.0/(1.0+F);
    calc_apw_pv(&pi,vi,rt,&(ad->pw)); // aw_pw.h
    *p=-(*p)*cf*ad->k2[did]+pi;
    for(l=0;l<3;l++) pv[l]=cf*pv[l]+vi[l];
  }
  else{
    cf=1.0/F;
    *p*=-cf*ad->k2[did];
    for(l=0;l<3;l++) pv[l]*=cf;
  } 
}

void pv_i_bd_aqd2(double complex *p,double complex *pv,int did,int t,double zeta_t,double eta_t,int type,AQD2 *ad)
{
  int td;
  double rt[3];
  int i;

  td=ad->bd.sb[did].sid[t];
  r_bd(rt,td,zeta_t,eta_t,&(ad->bd));

  if(did==0){
    calc_apw_pv(p,pv,rt,&(ad->pw)); // aw_pw.h
  }
  else {
    *p=0.0;
    for(i=0;i<3;i++) pv[i]=0.0;
  }
}
