#include "bem3_aw_qd2.h"

void solve_bieq_aqd2(AQD2 *qd)
{
  void initalize_cmd(CMD *cm);
  void finalize_cmd(CMD *cm);
  void create_cmatrix(CMD *cm,AQD2 *qd);
  void create_tmatrix_csr(CMD *cm,AQD2 *qd);
  void solve_tmatrix_csr(CMD *cm,AQD2 *qd);
  void solve_pv_bv(CMD *cm,AQD2 *qd);
  void solve_dpv_bv(CMD *cm,AQD2 *qd);
  void q_solve_pv_bv(CMD *cm,AQD2 *qd);

  time_t start,end,ms,me;
  CMD cm;

  printf("\nsolve acoustic wave boundary value \n");
  time(&start);

  // all coefficient matrix
  printf("  coefficient matrix          "); fflush(stdout);
  time(&ms);
  cm.type=2; // integration setting 0:4p GL,1:9p or 7p(triangular) GL, 2:GLN p GL, 3: GHN p GL
  cm.MN=qd->MN;
  initalize_cmd(&cm);
  create_cmatrix(&cm,qd);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  // velocity potential boundary value
  printf("  solve VP boundary value     "); fflush(stdout);
  time(&ms);
  cm.nn=(size_t)qd->bd.NN;
  cm.na=cm.nn*2;
  create_tmatrix_csr(&cm,qd);
  solve_tmatrix_csr(&cm,qd);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  // particle velocity boundary value
  printf("  solve PV boundary value     "); fflush(stdout);
  time(&ms);
  if(qd->bd.ps==1) q_solve_pv_bv(&cm,qd); // quasi-periodic boundary condition
  else solve_pv_bv(&cm,qd); // normal
  solve_dpv_bv(&cm,qd);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  finalize_cmd(&cm);
  time(&end);
  printf("Total elapsed time : %g (sec)\n",difftime(end,start));
}

///////////////////////////////////////////////////////////////////////
void initalize_cmd(CMD *cm)
{
  int i,MN;

  MN=cm->MN;
  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->thfn");
  cm->tdgfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdgfn");
  cm->tdhfn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdhfn");
  cm->tdffn=(char **)m_alloc2(MN+1,sizeof(char *),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdffn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.dat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.dat",i);
    cm->tdgfn[i]=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdgfn[i]");
    cm->tdhfn[i]=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdhfn[i]");
    cm->tdffn[i]=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->tdffn[i]");
    sprintf(cm->tdgfn[i],"tmpdG_%05d.dat",i);
    sprintf(cm->tdhfn[i],"tmpdH_%05d.dat",i);
    sprintf(cm->tdffn[i],"tmpdF_%05d.dat",i);
  }

  cm->aval=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->aval");
  cm->aptr=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->aptr");
  cm->aidx=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->aidx");
  sprintf(cm->aval,"tmpAval.dat");
  sprintf(cm->aptr,"tmpAprt.dat");
  sprintf(cm->aidx,"tmpAidx.dat");
  cm->b=(char *)m_alloc2(16,sizeof(char ),"bem3_aw_qd2_solve_bieq.c, initialize_cmd(), cm->b");
  sprintf(cm->b,"tmpB.dat");
}

void finalize_cmd(CMD *cm)
{
  int i;

  // delete temporary file
  for(i=0;i<=cm->MN;i++){
    remove(cm->tgfn[i]);
    remove(cm->thfn[i]);
    remove(cm->tdgfn[i]);
    remove(cm->tdhfn[i]);
    remove(cm->tdffn[i]);
  }
  remove(cm->aval);
  remove(cm->aptr);
  remove(cm->aidx);
  remove(cm->b);

  // free memory
  for(i=0;i<=cm->MN;i++){
    free(cm->tgfn[i]);    free(cm->thfn[i]);
    free(cm->tdgfn[i]);    free(cm->tdhfn[i]);  free(cm->tdffn[i]);
  }
  free(cm->tgfn);  free(cm->thfn);
  free(cm->tdgfn);  free(cm->tdhfn); free(cm->tdffn);

  cm->MN=0;

  free(cm->aval);
  free(cm->aptr);
  free(cm->aidx);
  cm->nn=0;
  cm->na=0;
  cm->nnz=0;
  free(cm->b);
}

void create_cmatrix(CMD *cm,AQD2 *qd)
{
  void create_cmatrix_domain(int did,CMD *cm,AQD2 *qd);

  int i;

  for(i=0;i<=qd->MN;i++) create_cmatrix_domain(i,cm,qd);
}

void create_cmatrix_domain(int did,CMD *cm,AQD2 *qd)
{
  FILE *fg,*fh,*fdg,*fdh,*fdf;
  double complex *tG,*tH,*tdG,*tdH,CC[9],dCz[9],dCe[9],kc;
  double F,vtz[3],vte[3],dFz,dFe,*tdF;
  size_t Ne,N,t,s,tn,tl,i;
  int td,sd;

  Ne=(size_t)qd->bd.sb[did].Ne;
  N=4*Ne;

  tG=(double complex *)m_alloc2(N*N,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),tG"); // malloc
  tH=(double complex *)m_alloc2(N*N,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),tH"); // malloc
  if((fg=fopen(cm->tgfn[did],"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),*fg. Failed to create %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),*fh. Failed to create %s file.\n",cm->thfn[did]);    exit(1);  }

  tdG=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),tdG"); // malloc
  tdH=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),tdH"); // malloc
  tdF=(double *)m_alloc2(3*N,sizeof(double),"bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),tdF"); // malloc
  if((fdg=fopen(cm->tdgfn[did],"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),*fdg. Failed to create %s file.\n",cm->tdgfn[did]);    exit(1);  }
  if((fdh=fopen(cm->tdhfn[did],"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),*fdh. Failed to create %s file.\n",cm->tdhfn[did]);    exit(1);  }
  if((fdf=fopen(cm->tdffn[did],"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_cmatrix_domain(),*fdf. Failed to create %s file.\n",cm->tdffn[did]);    exit(1);  }

  kc=(double complex)qd->k0[did];
  #pragma omp parallel for schedule(dynamic) private(td,tl,tn,vtz,vte,F,dFz,dFe,s,sd,CC,dCz,dCe,i) // parallel
  for(t=1;t<=Ne;t++){
    td=qd->bd.sb[did].sid[t];
    if( ELT3==check_element_type(td,&(qd->bd)) ) tl=3;
    else tl=4;

    for(tn=0;tn<tl;tn++){
      tz_te_bd_node(vtz,vte,td,tn,&(qd->bd));

      F=0.0;
      dFz=0.0;
      dFe=0.0;
      for(s=1;s<=Ne;s++){
        sd=qd->bd.sb[did].sid[s];

        if(did==0 && qd->bd.ps==1) q0_dcoef_bd_node_t2_cc(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(qd->bd));
        else dcoef_bd_node_t2(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(qd->bd));

        for(i=0;i<4;i++){
          tG[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+0];
          tH[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+4];

          tdG[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+4];
          tdG[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+4];
        }
        F+=creal(CC[8]);
        dFz+=creal(dCz[8]);
        dFe+=creal(dCe[8]);
      }

      if(did==0){
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=1.0+F;
        tdF[(t-1)*4*3+tn*3+0]=1.0+F;
      }
      else {
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=F;
        tdF[(t-1)*4*3+tn*3+0]=F;
      }
      tdF[(t-1)*4*3+tn*3+1]=dFz;
      tdF[(t-1)*4*3+tn*3+2]=dFe;
    }
  }

  fwrite(tG,sizeof(double complex),N*N,fg);
  fwrite(tH,sizeof(double complex),N*N,fh);
  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);

  fwrite(tdG,sizeof(double complex),2*N*N,fdg);
  fwrite(tdH,sizeof(double complex),2*N*N,fdh);
  fwrite(tdF,sizeof(double),3*N,fdf);
  fclose(fdg);
  fclose(fdh);
  fclose(fdf);
  free(tdG);
  free(tdH);
  free(tdF);
}

void create_tmatrix_csr(CMD *cm,AQD2 *qd)
{
  void create_matrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,AQD2 *qd);

  FILE *av,*ai,*ap,*b;

  size_t did;

  if((av=fopen(cm->aval,"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_tmatrix_csr(),*av. Failed to create %s file.\n",cm->aval);    exit(1);  }
  if((ai=fopen(cm->aidx,"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_tmatrix_csr(),*ai. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((ap=fopen(cm->aptr,"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_tmatrix_csr(),*ap. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((b=fopen(cm->b,"wb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_tmatrix_csr(),*b. Failed to create %s file.\n",cm->b);    exit(1);  }

  cm->nnz=0; // initialize nnz
  for(did=0;did<=cm->MN;did++) create_matrix_csr_dac(did,av,ap,ai,b,cm,qd);

  fwrite(&(cm->nnz),sizeof(size_t),1,ap);

  fclose(av);
  fclose(ai);
  fclose(ap);
  fclose(b);
}

void create_matrix_csr_dac(int did,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,AQD2 *qd)
{
  FILE *fg,*fh;

  double complex *tG,*tH,*tA,tB,k2m,k2s,hk2;
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l,cc,*ti;
  int td,sd;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),*fg. Failed to open %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),*fh. Failed to open %s file.\n",cm->thfn[did]);    exit(1);  }

  Ne=(size_t)qd->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),tH"); // malloc

  tA=(double complex *)m_alloc2(cm->na,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),tA"); // malloc
  ti=(size_t *)m_alloc2(cm->na,sizeof(size_t),"bem3_aw_qd2_solve_bieq.c, create_matrix_csr_dac(),ti"); // malloc

  for(t=1;t<=Ne;t++){
    td=qd->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      fread(tG,sizeof(double complex),Ne*4,fg);
      fread(tH,sizeof(double complex),Ne*4,fh);
      if( tn==3 && ELT3==check_element_type(td,&(qd->bd)) )   continue;

      fwrite(&(cm->nnz),sizeof(size_t),1,ap); // write A pointer
      for(l=0;l<cm->na;l++) tA[l]=0.0;
      tB=0.0;

      for(s=1;s<=Ne;s++){
        sd=qd->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)qd->bd.md[asd]; // main domain id
        sdid=(size_t)qd->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(qd->bd));                         

        if(did==mdid){ // main domain
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++) {
              tA[ cm->nn*0 + qd->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + qd->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++) {
              tA[ cm->nn*0 + qd->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*1 + qd->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
        } // end main domain
        else { // sub domain
          k2m=qd->k2[mdid];
          k2s=qd->k2[sdid];
          hk2=k2m/k2s;
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++){
              tA[ cm->nn*0 + qd->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + qd->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
              if(mdid==0) tB+=-hk2*tH[(s-1)*4+l]*qd->bd.Pi[asd][l]-tG[(s-1)*4+l]*qd->bd.dPi[asd][l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++){
              tA[ cm->nn*0 + qd->bd.eni[asd][l] ]=hk2*tH[(s-1)*4+l];
              tA[ cm->nn*1 + qd->bd.eni[asd][l] ]=    tG[(s-1)*4+l];
              if(mdid==0) tB+=-hk2*tH[(s-1)*4+l]*qd->bd.Pi[asd][l]-tG[(s-1)*4+l]*qd->bd.dPi[asd][l];
            }
          }
        } // end sub domain
      } // end for s

      // compress and store data
      cc=0;
      for(l=0;l<cm->na;l++){
        if( creal(tA[l])==0.0 && cimag(tA[l])==0.0) continue;
        tA[cc]=tA[l];
        ti[cc]=l;
        cc+=1;
      }
      fwrite(tA,sizeof(double complex),cc,av);
      fwrite(ti,sizeof(size_t),cc,ai);
      cm->nnz+=cc;

      fwrite(&tB,sizeof(double complex),1,b);

    } // end for tn
  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
  free(tA);
  free(ti);
}

void solve_tmatrix_csr(CMD *cm,AQD2 *qd)
{
  int b_mkl_solver_pardiso_d(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // d3qd1_mkl_solver.c
  int b_mkl_solver_pardiso_s(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // d3qd1_mkl_solver.c
  void tmatrix_bd_store(double complex *X,AQD2 *ad);

  int err;
  double complex *x; // results
  x=(double complex *)m_alloc2(cm->na,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, solve_tmatrix_csr(),x");

  if(PARDISO_PREC==0) err=b_mkl_solver_pardiso_d(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  else err=b_mkl_solver_pardiso_s(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  if(err!=0){
    printf("bem3_aw_qd2_solve_bieq.c, solve_tmatrix_csr(), b_mkl_solver_pardiso() returned error. err=%d. Exit...\n",err);
    exit(1);
  }

  // store data
  tmatrix_bd_store(x,qd);

  free(x);
}

void tmatrix_bd_store(double complex *X,AQD2 *qd)
{
  double complex k2m,k2s,hk2;
  size_t d,s,l,nn;
  int sd,asd,etype,mdid,sdid;

  nn=qd->bd.NN;

  for(d=0;d<=qd->MN;d++){
    for(s=1;s<=qd->bd.sb[d].Ne;s++){
      sd=qd->bd.sb[d].sid[s];
      asd=abs(sd);
      mdid=qd->bd.md[asd];
      sdid=qd->bd.sd[asd];
      etype=check_element_type(sd,&(qd->bd));

      if(d==mdid){ // main domain
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            qd->bd.sb[d]. P[s][l]=X[nn*0+qd->bd.eni[asd][l]];
            qd->bd.sb[d].dP[s][l]=X[nn*1+qd->bd.eni[asd][l]];
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            qd->bd.sb[d]. P[s][l]=X[nn*0+qd->bd.eni[asd][l]];
            qd->bd.sb[d].dP[s][l]=X[nn*1+qd->bd.eni[asd][l]];
          }
        }
      }
      else { // subdomain
        k2m=qd->k2[mdid];
        k2s=qd->k2[sdid];
        hk2=k2m/k2s;
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            qd->bd.sb[d]. P[s][l]= hk2*X[nn*0+qd->bd.eni[asd][l]];
            qd->bd.sb[d].dP[s][l]=-1.0*X[nn*1+qd->bd.eni[asd][l]];
            if(mdid==0){
              qd->bd.sb[d]. P[s][l]+= hk2*qd->bd. Pi[asd][l];
              qd->bd.sb[d].dP[s][l]+=-1.0*qd->bd.dPi[asd][l];
            }
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            qd->bd.sb[d]. P[s][l]= hk2*X[nn*0+qd->bd.eni[asd][l]];
            qd->bd.sb[d].dP[s][l]=-1.0*X[nn*1+qd->bd.eni[asd][l]];
            if(mdid==0){
              qd->bd.sb[d]. P[s][l]+= hk2*qd->bd. Pi[asd][l];
              qd->bd.sb[d].dP[s][l]+=-1.0*qd->bd.dPi[asd][l];
            }
          }
        }
      } // end subdomain
    }       // end s
  } // end d
}

void solve_pv_bv(CMD *cm,AQD2 *qd)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dpdt_fcm(double complex *dPz,double complex *dPe,int t,int tn,int Ne,
      double complex **P,double complex **dP,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dPz,dPe;
  double sig,tdF[3],vtz[3],vte[3],vn[3],a[9];
  //size_t Ne,d,t,atd,tn,s,sn,i;
  int d,Ne,t,tn,at,atd,i;

  for(d=0;d<=qd->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, solve_pv_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, solve_pv_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, solve_pv_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=qd->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, solve_pv_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, solve_pv_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=qd->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(qd->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(qd->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*qd->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dPdtz,dPdte
        dpdt_fcm(&dPz,&dPe,t,tn,Ne,qd->bd.sb[d].P,qd->bd.sb[d].dP,tdG,tdH,tdF);
        // particle velocity
        for(i=0;i<3;i++) qd->bd.sb[d].pv[t][tn][i]=-(qd->bd.sb[d].dP[t][tn]*a[i*3+0]+dPz*a[i*3+1]+dPe*a[i*3+2]);

      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d
}

void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte)
{
  double i_det;
  int i;

  i_det=1.0/(vn[0]*vtz[1]*vte[2]+vn[1]*vtz[2]*vte[0]+vn[2]*vtz[0]*vte[1]
             -( vn[2]*vtz[1]*vte[0]+vn[1]*vtz[0]*vte[2]+vn[0]*vtz[2]*vte[1]));

  ma[0*3+0]= vtz[1]*vte[2]-vtz[2]*vte[1];
  ma[0*3+1]=- vn[1]*vte[2]+ vn[2]*vte[1];
  ma[0*3+2]=  vn[1]*vtz[2]- vn[2]*vtz[1];

  ma[1*3+0]=-vtz[0]*vte[2]+vtz[2]*vte[0];
  ma[1*3+1]=  vn[0]*vte[2]- vn[2]*vte[0];
  ma[1*3+2]=- vn[0]*vtz[2]+ vn[2]*vtz[0];

  ma[2*3+0]= vtz[0]*vte[1]-vtz[1]*vte[0];
  ma[2*3+1]=- vn[0]*vte[1]+ vn[1]*vte[0];
  ma[2*3+2]=  vn[0]*vtz[1]- vn[1]*vtz[0];

  for(i=0;i<9;i++) ma[i]*=i_det;
}

void dpdt_fcm(double complex *dPz,double complex *dPe,int t,int tn,int Ne,double complex **P,double complex **dP,double complex *tdG,double complex *tdH,double *tdF)
{
  int s,sn;

  *dPz=0.0;
  *dPe=0.0;
  for(s=1;s<=Ne;s++){
    for(sn=0;sn<4;sn++){
      *dPz+=tdG[Ne*4*0+4*(s-1)+sn]*dP[s][sn]-tdH[Ne*4*0+4*(s-1)+sn]*P[s][sn];
      *dPe+=tdG[Ne*4*1+4*(s-1)+sn]*dP[s][sn]-tdH[Ne*4*1+4*(s-1)+sn]*P[s][sn];
    }
  }
  *dPz=(*dPz-P[t][tn]*tdF[1])/tdF[0];
  *dPe=(*dPe-P[t][tn]*tdF[2])/tdF[0];
}


void solve_dpv_bv(CMD *cm,AQD2 *qd)
{
  FILE *fg,*fh;
  MKL_Complex16 *A,*B;
  double complex *tG,*tH;
  int Ne,d,t,atd,tn,nn,nc,nr,s,ns,asd,i;
  MKL_INT *ipiv,nrhs,lda,ldb,info;

  nrhs=3;

  for(d=0;d<=qd->MN;d++){
    Ne=qd->bd.sb[d].Ne;
    if(Ne==0) continue;

    if((fg=fopen(cm->tgfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),*fg. Failed to open %s file.\n",cm->tgfn[d]);    exit(1);  }
    if((fh=fopen(cm->thfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),*fh. Failed to open %s file.\n",cm->thfn[d]);    exit(1);  }

    // matrix size
    nn=0;
    for(s=1;s<=Ne;s++){
      asd=abs(qd->bd.sb[d].sid[s]);
      if( ELT3==check_element_type(asd,&(qd->bd)) ) nn+=3;
      else nn+=4;
    }
    lda=nn;
    ldb=nrhs;

    A =(MKL_Complex16 *)m_alloc2(nn*lda,sizeof(MKL_Complex16),"bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),A"); // malloc
    B =(MKL_Complex16 *)m_alloc2(nn*ldb,sizeof(MKL_Complex16),"bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),B"); // malloc
    ipiv=(MKL_INT *)m_alloc2(nn,sizeof(MKL_INT),"bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),ipiv"); // malloc
    tG=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),tG"); // malloc
    tH=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(),tH"); // malloc

    nc=0;
    for(t=1;t<=Ne;t++){
      atd=abs(qd->bd.sb[d].sid[t]);

      for(tn=0;tn<4;tn++){
        fread(tG,sizeof(double complex),Ne*4,fg);
        fread(tH,sizeof(double complex),Ne*4,fh);
        if( tn==3 && ELT3==check_element_type(atd,&(qd->bd)) )  continue;

        for(i=0;i<ldb;i++){
          B[ldb*nc+i].real=0.0;
          B[ldb*nc+i].imag=0.0;
        }
        nr=0;
        for(s=1;s<=Ne;s++){
          asd=abs(qd->bd.sb[d].sid[s]);
          if( ELT4==check_element_type(asd,&(qd->bd)) ){
            for(ns=0;ns<4;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[ldb*nc+i].real+=creal(tH[4*(s-1)+ns]*qd->bd.sb[d].pv[s][ns][i]);
                B[ldb*nc+i].imag+=cimag(tH[4*(s-1)+ns]*qd->bd.sb[d].pv[s][ns][i]);
              }
            }
          }
          else {
            for(ns=0;ns<3;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[ldb*nc+i].real+=creal(tH[4*(s-1)+ns]*qd->bd.sb[d].pv[s][ns][i]);
                B[ldb*nc+i].imag+=cimag(tH[4*(s-1)+ns]*qd->bd.sb[d].pv[s][ns][i]);
              }
            }
          }
        } // end for s
        nc+=1;
      } // end for tn

    } // end for t

    // solve mkl lapack
    info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, nn, nrhs, A, lda, ipiv, B, ldb );
    if(info!=0){
      printf("bem3_aw_qd2_solve_bieq.c, solve_dpv_bv(), LAPACKE_zgesv(). info=%lld error. Exit...\n",info);
      exit(1);
    }

    // store data
    nc=0;
    for(s=1;s<=Ne;s++){
      asd=abs(qd->bd.sb[d].sid[s]);
      if( ELT4==check_element_type(asd,&(qd->bd)) ){
        for(ns=0;ns<4;ns++){
          for(i=0;i<3;i++) {
            qd->bd.sb[d].dpv[s][ns][i]=B[ldb*nc+i].real+B[ldb*nc+i].imag*I;
          }
          nc+=1;
        }
      }
      else {
        for(ns=0;ns<3;ns++){
          for(i=0;i<3;i++) {
            qd->bd.sb[d].dpv[s][ns][i]=B[ldb*nc+i].real+B[ldb*nc+i].imag*I;
          }
          nc+=1;
        }
      }
    }

    fclose(fg);
    fclose(fh);
    free(A);
    free(B);
    free(ipiv);
    free(tG);
    free(tH);
  } // end for d
}

void q_solve_pv_bv(CMD *cm,AQD2 *qd)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dpdt_fcm(double complex *dPz,double complex *dPe,int t,int tn,int Ne,
      double complex **P,double complex **dP,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dPz,dPe,pi,pvi[3],*tpv,rhs[3];
  double sig,tdF[3],vtz[3],vte[3],vn[3],a[9],rtz[3],rte[3];
  int d,Ne,t,tn,at,atd,i;

  tpv=(double complex *)m_alloc2(qd->bd.NN*3,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),tpv"); // malloc

  for(d=1;d<=qd->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=qd->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"bem3_aw_qd2_solve_bieq.c, q_solve_pv_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=qd->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(qd->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(qd->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*qd->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dPdtz,dPdte
        dpdt_fcm(&dPz,&dPe,t,tn,Ne,qd->bd.sb[d].P,qd->bd.sb[d].dP,tdG,tdH,tdF);
        // particle velocity
        for(i=0;i<3;i++) qd->bd.sb[d].pv[t][tn][i]=-(qd->bd.sb[d].dP[t][tn]*a[i*3+0]+dPz*a[i*3+1]+dPe*a[i*3+2]);

        // for scattered field
        if(qd->bd.md[atd]==0){
          for(i=0;i<3;i++) rhs[i]=0.0;
          calc_apw_pv(&pi,pvi,qd->bd.ren[atd][tn],&(qd->pw));
          for(i=0;i<3;i++){
            rtz[i]=qd->rho0[0]*vtz[i];
            rte[i]=qd->rho0[0]*vte[i];
            rhs[0]=vn[0]*(qd->bd.sb[d].pv[t][tn][0]-pvi[0])
                  +vn[1]*(qd->bd.sb[d].pv[t][tn][1]-pvi[1])
                  +vn[2]*(qd->bd.sb[d].pv[t][tn][2]-pvi[2]);
            rhs[1]=vtz[0]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][0]-qd->rho0[0]*pvi[0])
                  +vtz[1]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][1]-qd->rho0[0]*pvi[1])
                  +vtz[2]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][2]-qd->rho0[0]*pvi[2]);
            rhs[2]=vte[0]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][0]-qd->rho0[0]*pvi[0])
                  +vte[1]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][1]-qd->rho0[0]*pvi[1])
                  +vte[2]*(qd->rho0[d]*qd->bd.sb[d].pv[t][tn][2]-qd->rho0[0]*pvi[2]);
          }
          inv_matrix_33(a,vn,rtz,rte);
          for(i=0;i<3;i++) tpv[qd->bd.eni[atd][tn]*3+i]=a[i*3+0]*rhs[0]+a[i*3+1]*rhs[1]+a[i*3+2]*rhs[2];
        }
      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d

  // store scattered field data
  for(t=1;t<=qd->bd.sb[0].Ne;t++){
    at=qd->bd.sb[0].sid[t];
    atd=abs(at);
    for(tn=0;tn<4;tn++){
      if( tn==3 && ELT3==check_element_type(atd,&(qd->bd)) )  continue;
      for(i=0;i<3;i++) qd->bd.sb[0].pv[t][tn][i]=tpv[qd->bd.eni[atd][tn]*3+i];
     }
  }

  free(tpv);
}
