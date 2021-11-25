#include "bem3_aw_qd2.h"

void read_aqd2(int argc,char **argv,AQD2 *qd)
{
  void filename_chk(int argc,char **argv);
  void read_periodicity_data(char *ped_fn,AQD2 *qd);
  void read_medium_data(char *med_fn,AQD2 *qd);
  void read_mesh_data(char *msh_fn,AQD2 *qd);
  
  filename_chk(argc,argv);
  read_data_apw(argv[1],&(qd->pw)); // aw_pw.h
  read_periodicity_data(argv[2],qd);
  read_medium_data(argv[3],qd);
  read_mesh_data(argv[4],qd);
  
  if(argc==13){
    qd->rv[0]=atof(argv[ 6]);
    qd->rv[1]=atof(argv[ 7]);
    qd->rv[2]=atof(argv[ 8]);
    qd->th   =atof(argv[ 9]);
    qd->tv[0]=atof(argv[10]);
    qd->tv[1]=atof(argv[11]);
    qd->tv[2]=atof(argv[12]);
  }
  else {
    qd->rv[0]=1.0;
    qd->rv[1]=0.0;
    qd->rv[2]=0.0;
    qd->th=0.0;
    qd->tv[0]=0.0;
    qd->tv[1]=0.0;
    qd->tv[2]=0.0;
  }
}

void print_aqd2(AQD2 *qd)
{
  void print_periodicity_data(AQD2 *qd);
  void print_medium_data(AQD2 *qd);
  void print_mesh_data(AQD2 *qd);

  print_data_apw(&(qd->pw));
  print_periodicity_data(qd);
  print_medium_data(qd);
  print_mesh_data(qd);
  
  if(qd->th!=0.0 || vabs_d(qd->tv)!=0.0){
    printf("-- rotation and translation settings --\n");
    if(qd->th!=0.0){
      printf("vector defining rotation axis :(% 8.7g,% 8.7g,% 8.7g)\n",qd->rv[0],qd->rv[1],qd->rv[2]);
      printf("rotation angle           [rad]: %8.7g\n",qd->th);
    }
    if(vabs_d(qd->tv)!=0.0){
      printf("translation vector            :(%8.7g,%8.7g,%8.7g)\n",qd->tv[0],qd->tv[1],qd->tv[2]);
    }
  }
  printf("\n");
}

void initialize_aqd2(AQD2 *qd)
{
  void rotation_translation_obj(double *rv,double th,double *tv,AQD2 *ad);
  void init_elem_const(BOUD *bd);
  void malloc_sub_domain(AQD2 *ad);
  void init_sub_domain(AQD2 *ad);
  void init_boundary_data(AQD2 *ad);

  int i;

  // incident field
  setup_apw(&(qd->pw));
  // QPDT
  qd->bd.qd.k=qd->pw.wd.k0;
  qd->bd.qd.vk[0]=qd->pw.wd.kx;
  qd->bd.qd.vk[1]=qd->pw.wd.ky;
  qd->bd.qd.vk[2]=qd->pw.wd.kz;
  if(qd->bd.ps!=0){
    if(qd->bd.qd.vk[0]==0.0 || qd->bd.qd.vk[1]==0.0 ){ // kx==0.0 or ky==0.0
      printf("kx and ky must not be equal to 0 (must satisfy quasi-peridic boundary condition). Exit...\n");
      exit(1);
    }
    if(qd->bd.qd.vk[2]==0.0){ // kz=0.0
      printf("kz must not be equal to 0 (The wave source and the objcet must not overlap). Exit...\n");
    }
  }
  // medium
  qd->rho0[0]=qd->pw.rho0;
  qd->c0[0]=qd->pw.c0;
  qd->k0[0]=qd->pw.wd.k0;
  qd->K0[0]=qd->pw.wd.K0;
  qd->k1[0]=qd->pw.wd.k1;
  qd->k2[0]=qd->pw.wd.k2;
  for(i=1;i<=qd->MN;i++){
    qd->k0[i]=qd->pw.wd.omega/qd->c0[i];
    qd->K0[i]=qd->rho0[i]*qd->c0[i]*qd->c0[i];
    qd->k1[i]=I*qd->pw.wd.omega/qd->K0[i];
    qd->k2[i]=I*qd->pw.wd.omega*qd->rho0[i];
  }
  // rotation and translation object
  if(qd->th!=0.0 || vabs_d(qd->tv)!=0.0) rotation_translation_obj(qd->rv,qd->th,qd->tv,qd);  
  // element constant
  init_elem_const(&(qd->bd));
  // sub domain
  malloc_sub_domain(qd);
  init_sub_domain(qd);
  // boundary data
  init_boundary_data(qd);

}

void finalize_aqd2(AQD2 *qd)
{
  void mfree_sub_domain(AQD2 *qd);
  void mfree_elem(BOUD *bd);
  void mfree_node(BOUD *bd);

  mfree_sub_domain(qd);

  mfree_elem(&(qd->bd));
  mfree_node(&(qd->bd));

  free(qd->rho0);
  free(qd->c0);
  free(qd->k0);
  free(qd->K0);
  free(qd->k1);
  free(qd->k2);
}

int domain_id_aqd2(double *ro,AQD2 *qd)
{
  int domain_id(double *rt,AQD2 *qd);
  double fid_calc_solid_angle(int type,double r[4][3],int *flg);

  double rv[4][3],omega,rt[3];
  double *Og=(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd2.c, domain_id_aqd2()");
  int i,j,k,d,flg,l1,l2;

  if(qd->bd.ps==0) return domain_id(ro,qd);

  domain_id_l_aqd2(&l1,&l2,ro,qd);
  rt[0]=ro[0]-(double)l1*qd->bd.qd.vd1[0]-(double)l2*qd->bd.qd.vd2[0];
  rt[1]=ro[1]-(double)l1*qd->bd.qd.vd1[1]-(double)l2*qd->bd.qd.vd2[1];
  rt[2]=ro[2];

  for(d=0;d<qd->MN+1;d++){
    for(i=1;i<=qd->bd.sb[d].Ne;i++){
      if(qd->bd.sb[d].sid[i]>0){
        // read node data
        for(j=0;j<4;j++)
          for(k=0;k<3;k++) rv[j][k]=qd->bd.rn[qd->bd.ed[qd->bd.sb[d].sid[i]][j]][k]-rt[k];
        omega=fid_calc_solid_angle(check_element_type(qd->bd.sb[d].sid[i],&(qd->bd)),rv,&flg);
        if(flg<0){ // on boundary
          free(Og);
          return d;
        }
        Og[d]+=omega;
        Og[qd->bd.sd[qd->bd.sb[d].sid[i]]]-=omega;

      } // end if
    }
    if(d==0 && fabs(Og[d])<2.0*M_PI){ // opened region
      free(Og);
      return d;
    }
    else if(Og[d]>2.0*M_PI){ // closed region
      free(Og);
      return d;
    }
  }

  free(Og);
  return -1; // error 
}

void domain_id_l_aqd2(int *l1,int *l2,double *rt,AQD2 *qd)
{
  double i_detd;
  
  i_detd=1.0/(qd->bd.qd.vd1[0]*qd->bd.qd.vd2[1]-qd->bd.qd.vd2[0]*qd->bd.qd.vd1[1]);
  *l1=(int)floor( (rt[0]*qd->bd.qd.vd2[1]-rt[1]*qd->bd.qd.vd2[0])*i_detd+0.5);
  *l2=(int)floor(-(rt[0]*qd->bd.qd.vd1[1]-rt[1]*qd->bd.qd.vd1[0])*i_detd+0.5);
}

void dat_write_aqd2(char *fname,AQD2 *qd)
{
  FILE *fp;
  int i,j,d;

  if((fp=fopen(fname,"wb"))==NULL){    printf("bem3_aw_qd1.c, dat_write_aqd1(), Failed to create the %s file.\n",fname);    exit(1);  }

  fwrite(qd,sizeof(AQD2),1,fp);
  // material def
  fwrite(qd->rho0,sizeof(double),qd->MN+1,fp);
  fwrite(qd->c0,sizeof(double),qd->MN+1,fp);
  fwrite(qd->k0,sizeof(double),qd->MN+1,fp);
  fwrite(qd->K0,sizeof(double),qd->MN+1,fp);
  fwrite(qd->k1,sizeof(double complex),qd->MN+1,fp);
  fwrite(qd->k2,sizeof(double complex),qd->MN+1,fp);
  // incident field
  fwrite(&(qd->pw),sizeof(Apw),1,fp);
  // BOUD
  for(i=0;i<=qd->bd.Nn;i++) fwrite(qd->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) fwrite(qd->bd.ed[i],sizeof(int),4,fp);
  for(i=0;i<=qd->bd.Ne;i++) fwrite(qd->bd.eni[i],sizeof(int),4,fp);
  fwrite(qd->bd.md,sizeof(int),qd->bd.Ne+1,fp);
  fwrite(qd->bd.sd,sizeof(int),qd->bd.Ne+1,fp);
  fwrite(qd->bd.gd,sizeof(int),qd->bd.Ne+1,fp);
  for(i=0;i<=qd->bd.Ne;i++) for(j=0;j<4;j++) fwrite(qd->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) for(j=0;j<4;j++) fwrite(qd->bd.wen[i][j],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) fwrite(qd->bd. Pi[i],sizeof(double complex),4,fp);
  for(i=0;i<=qd->bd.Ne;i++) fwrite(qd->bd.dPi[i],sizeof(double complex),4,fp);
  // sub domain data
  for(d=0;d<=qd->MN;d++){
    fwrite(&(qd->bd.sb[d].Ne),sizeof(int),1,fp);
    fwrite(qd->bd.sb[d].sid,sizeof(int),qd->bd.sb[d].Ne+1,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) fwrite(qd->bd.sb[d]. P[i],sizeof(double complex),4,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) fwrite(qd->bd.sb[d].dP[i],sizeof(double complex),4,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(qd->bd.sb[d]. pv[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fwrite(qd->bd.sb[d].dpv[i][j],sizeof(double complex),3,fp);
  }

  fclose(fp);
}

void dat_read_aqd2(char *fname,AQD2 *qd)
{
  void malloc_node(BOUD *bd);  
  void malloc_elem(BOUD *bd); 
  void init_elem_const(BOUD *bd);
  void malloc_sub_domain(AQD2 *md);

  FILE *fp;
  int i,j,d,tmp;

  if((fp=fopen(fname,"rb"))==NULL){    printf("bem3_aw_qd1.c, dat_read_aqd1(), Failed to open the %s file.\n",fname);    exit(1);  }

  fread(qd,sizeof(AQD2),1,fp);
  // material def
  qd->rho0=(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd1.c, dat_read_aqd1(),ad->rho0");
  qd->c0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd1.c, dat_read_aqd1(),ad->c0");
  qd->k0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd1.c, dat_read_aqd1(),ad->k0");
  qd->K0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd1.c, dat_read_aqd1(),ad->K0");
  qd->k1=(double complex *)m_alloc2(qd->MN+1,sizeof(double complex),"bem3_aw_qd1.c, dat_read_aqd1(),ad->k1");
  qd->k2=(double complex *)m_alloc2(qd->MN+1,sizeof(double complex),"bem3_aw_qd1.c, dat_read_aqd1(),ad->k2");
  fread(qd->rho0,sizeof(double),qd->MN+1,fp);
  fread(qd->c0,sizeof(double),qd->MN+1,fp);
  fread(qd->k0,sizeof(double),qd->MN+1,fp);
  fread(qd->K0,sizeof(double),qd->MN+1,fp);
  fread(qd->k1,sizeof(double complex),qd->MN+1,fp);
  fread(qd->k2,sizeof(double complex),qd->MN+1,fp);
  // incident field
  fread(&(qd->pw),sizeof(Apw),1,fp);
  // BOUD
  malloc_node(&(qd->bd)); // malloc
  malloc_elem(&(qd->bd)); // malloc
  for(i=0;i<=qd->bd.Nn;i++) fread(qd->bd.rn[i],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) fread(qd->bd.ed[i],sizeof(int),4,fp);
  for(i=0;i<=qd->bd.Ne;i++) fread(qd->bd.eni[i],sizeof(int),4,fp);
  fread(qd->bd.md,sizeof(int),qd->bd.Ne+1,fp);
  fread(qd->bd.sd,sizeof(int),qd->bd.Ne+1,fp);
  fread(qd->bd.gd,sizeof(int),qd->bd.Ne+1,fp);
  for(i=0;i<=qd->bd.Ne;i++) for(j=0;j<4;j++) fread(qd->bd.ren[i][j],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) for(j=0;j<4;j++) fread(qd->bd.wen[i][j],sizeof(double),3,fp);
  for(i=0;i<=qd->bd.Ne;i++) fread(qd->bd. Pi[i],sizeof(double complex),4,fp);
  for(i=0;i<=qd->bd.Ne;i++) fread(qd->bd.dPi[i],sizeof(double complex),4,fp);
  init_elem_const(&(qd->bd)); // setup
  // sub domain data
  malloc_sub_domain(qd); // malloc
  for(d=0;d<=qd->MN;d++){
    fread(&tmp,sizeof(int),1,fp);
    fread(qd->bd.sb[d].sid,sizeof(int),qd->bd.sb[d].Ne+1,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) fread(qd->bd.sb[d]. P[i],sizeof(double complex),4,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) fread(qd->bd.sb[d].dP[i],sizeof(double complex),4,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(qd->bd.sb[d]. pv[i][j],sizeof(double complex),3,fp);
    for(i=0;i<=qd->bd.sb[d].Ne;i++) for(j=0;j<4;j++) fread(qd->bd.sb[d].dpv[i][j],sizeof(double complex),3,fp);
  }
  fclose(fp);
}

///////////////////////////////////////////////////////////////////////
void filename_chk(int argc,char **argv)
{
  if(argc!=6 && argc!=13){
    printf("This program needs command line arguments as follows.\n");
    printf("%s plane_wave_datafile_name periodicity_datafile_name medium_datafile_name mesh_datafile_name output_datafile_name\n",argv[0]);
    printf("[rv_x rv_y rv_z theta tr_x tr_y tr_z](optional)\n");
    printf("rv : vector defining rotation axis, theta : rotation angle ( using Rodrigues' rotation formula ), tr : translation vector\n");
    exit(0);
  }
}

void rotation_translation_obj(double *rv,double th,double *tv,AQD2 *ad)
{
  double ct,st,r[3],M[9],nv[3];
  size_t s,i;

  nv[0]=rv[0];
  nv[1]=rv[1];
  nv[2]=rv[2];
  vuni_d(nv);

  // rotation matrix
  st=sin(th);
  ct=cos(th);
  M[0]=ct+nv[0]*nv[0]*(1.0-ct);
  M[1]=nv[0]*nv[1]*(1.0-ct)-nv[2]*st;
  M[2]=nv[2]*nv[0]*(1.0-ct)+nv[1]*st;
  M[3]=nv[0]*nv[1]*(1.0-ct)+nv[2]*st;
  M[4]=ct+nv[1]*nv[1]*(1.0-ct);
  M[5]=nv[1]*nv[2]*(1.0-ct)-nv[0]*st;
  M[6]=nv[2]*nv[0]*(1.0-ct)-nv[1]*st;
  M[7]=nv[1]*nv[2]*(1.0-ct)+nv[0]*st;
  M[8]=ct+nv[2]*nv[2]*(1.0-ct);

  for(s=1;s<=ad->bd.Nn;s++){
    for(i=0;i<3;i++) r[i]=M[3*i+0]*ad->bd.rn[s][0]+M[3*i+1]*ad->bd.rn[s][1]+M[3*i+2]*ad->bd.rn[s][2]+tv[i];
    for(i=0;i<3;i++) ad->bd.rn[s][i]=r[i];
  }
}

void read_periodicity_data(char *ped_fn,AQD2 *qd)
{
  FILE *fp;
  char buf[255]="";
  double td0,td1;
  int ti;

  if((fp=fopen(ped_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",ped_fn);    exit(1);  }
  strcpy(qd->ped_fn,ped_fn);
  fgets(buf,256,fp);
  fgets(buf,256,fp);
  fscanf(fp,"%d\n",&ti); qd->bd.ps=ti;
  fscanf(fp,"%lf %lf",&td0,&td1); qd->bd.qd.vd1[0]=td0; qd->bd.qd.vd1[1]=td1;
  fscanf(fp,"%lf %lf",&td0,&td1); qd->bd.qd.vd2[0]=td0; qd->bd.qd.vd2[1]=td1;
}

void print_periodicity_data(AQD2 *qd)
{
  printf("--- periodicity data ---\n");
  printf("periodicity data file name = %s\n",qd->ped_fn);
  printf("periodicity setting        = %10d\n",qd->bd.ps);
  printf("lattice constant d1        = (% 7g,% 7g)\n",qd->bd.qd.vd1[0],qd->bd.qd.vd1[1]);
  printf("lattice constant d2        = (% 7g,% 7g)\n",qd->bd.qd.vd2[0],qd->bd.qd.vd2[1]);
  printf("\n");
}

void read_medium_data(char *med_fn,AQD2 *qd)
{
  FILE *fp;
  double td;
  char buf[256]="";
  int i,ti;

  if((fp=fopen(med_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",med_fn);    exit(1);  }
  strcpy(qd->med_fn,med_fn);
  fgets(buf,256,fp);
  fgets(buf,256,fp);
  fscanf(fp,"%d\n",&ti);  qd->MN=ti;

  // malloc
  qd->rho0=(double *)m_alloc2(qd->MN+1,sizeof(double),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->rho0");
  qd->c0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->c0");
  qd->k0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->k0");
  qd->K0  =(double *)m_alloc2(qd->MN+1,sizeof(double),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->K0");
  qd->k1=(double complex *)m_alloc2(qd->MN+1,sizeof(double complex),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->k1");
  qd->k2=(double complex *)m_alloc2(qd->MN+1,sizeof(double complex),"memory allocation error! bem3_aw_qd2.c,read_medium_data(),qd->k2");

  fgets(buf,256,fp);
  for(i=1;i<=qd->MN;i++){
    fscanf(fp,"%lf",&td); qd->rho0[i]=td;
    fscanf(fp,"%lf",&td); qd->c0[i]=td;
  }

  fclose(fp);
}

void print_medium_data(AQD2 *qd)
{
  int i;
  printf("--- medium data ---\n");
  printf("medium data file name      = %s\n",qd->med_fn);
  for(i=1;i<=qd->MN;i++){
    printf("domain id %2d rho0 [kg/m^3] =%8.6g\n",i,qd->rho0[i]);
    printf("             c0      [m/s] =%8.6g\n",qd->c0[i]);
  }
  printf("\n");
}

void read_mesh_data(char *msh_fn,AQD2 *qd)
{
  void malloc_node(BOUD *bd);
  void malloc_elem(BOUD *bd);

  FILE *fp;
  char buf[256]="";
  double td;
  int ti,i,j,ti2,etype,tmpi,nc;

  if((fp=fopen(msh_fn,"rt"))==NULL){    printf("Can not open the '%s' file. Exit...\n",msh_fn);    exit(1);  }
  strcpy(qd->msh_fn,msh_fn);
  fgets(buf,256,fp);
  fscanf(fp,"%lf",&td);
  // check file version
  if(td<MSHVER){
    printf("This program supports mesh file version %g later. Reading file version is %g. Exit...\n",MSHVER,td);
    fclose(fp);
    exit(1);
  }
  fscanf(fp,"%d",&ti);
  //check data format
  if(ti!=MSHASCI){
    printf("This program supports 'ASCII' data format mesh file. Exit...\n");
    fclose(fp);
    exit(1);
  }
  fscanf(fp,"%d\n",&ti);
  //check data precision
  if(ti!=MSHPREC){
    printf("This program supports double precision mesh data. Exit...\n");
    fclose(fp);
    exit(1);
  }
  fgets(buf,256,fp);
  fgets(buf,256,fp);

  fscanf(fp,"%d\n",&ti); //printf("Nn=%d\n",ti);
  qd->bd.Nn=ti;
  malloc_node(&(qd->bd));
  for(i=1;i<=qd->bd.Nn;i++){
    fscanf(fp,"%d",&ti);    if(ti!=i)       printf("bad id %d\n",ti);
    fscanf(fp,"%lf",&td);   qd->bd.rn[i][0]=td;
    fscanf(fp,"%lf",&td);   qd->bd.rn[i][1]=td;
    fscanf(fp,"%lf\n",&td); qd->bd.rn[i][2]=td;
  }
  fgets(buf,256,fp);
  fgets(buf,256,fp);
  //for(i=1;i<=md->bd.Nn;i++)    printf("r[%d]=(%15.14e,%15.14e,%15.14e)\n",i,md->bd.rn[i][0],md->bd.rn[i][1],md->bd.rn[i][2]); // test

  fscanf(fp,"%d",&ti);
  qd->bd.Ne=ti/2; //printf("Ne=%d\n",md->bd.Ne);
  malloc_elem(&(qd->bd));

  nc=0;
  for(i=1;i<=qd->bd.Ne;i++){
    // element id
    fscanf(fp,"%d",&ti);
    if(ti!=i*2-1){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    fscanf(fp,"%d",&ti);
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    fscanf(fp,"%d",&ti);
    for(j=0;j<ti;j++){
      fscanf(fp,"%d",&ti2);
      if(j==0){ // domain id ( gmsh physical entity)
        if(ti2==OPENDID) ti2=0;
        if(qd->MN>=ti2) qd->bd.md[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data. exit..\n",ti2);
          exit(1);
        }
      }
      else if(j==1){ // group id ( elementary geometrical entity )
        qd->bd.gd[i]=ti2;
      }
    }
    // node id
    fscanf(fp,"%d",&ti);  qd->bd.ed[i][0]=ti;
    fscanf(fp,"%d",&ti);  qd->bd.ed[i][1]=ti;
    fscanf(fp,"%d",&ti);  qd->bd.ed[i][2]=ti;
    if(etype==ELT3) qd->bd.ed[i][3]=0;
    else {
      fscanf(fp,"%d",&ti);  qd->bd.ed[i][3]=ti;
    }
    // element node id
    if(etype==ELT3){
      qd->bd.eni[i][0]=nc++;
      qd->bd.eni[i][1]=nc++;
      qd->bd.eni[i][2]=nc++;
      qd->bd.eni[i][3]=-1;
    }
    else {
      qd->bd.eni[i][0]=nc++;
      qd->bd.eni[i][1]=nc++;
      qd->bd.eni[i][2]=nc++;
      qd->bd.eni[i][3]=nc++;
    }

    // element id
    fscanf(fp,"%d",&ti);
    if(ti!=i*2){
      printf("bad id :%d. Exit...\n",ti);
      exit(1);
    }
    // element type
    fscanf(fp,"%d",&ti);
    etype=ti;
    if(ti!=ELT3 && ti!=ELT4){
      printf("bad element type. element type must be %d or %d. Exit...\n",ELT3,ELT4);
      exit(1);
    }
    // number of tags
    fscanf(fp,"%d",&ti);
    for(j=0;j<ti;j++){
      fscanf(fp,"%d",&ti2);
      if(j==0){ // domain id
        if(ti2==OPENDID) ti2=0;
        if(qd->MN>=ti2) qd->bd.sd[i]=ti2;
        else {
          printf("domain id %d is not defined medium data. check domain and medium data! exit..\n",ti2);
          exit(1);
        }
      }
    }
    // check node id
    if(etype==ELT3){
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    else {
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][0]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][3]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][2]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
      fscanf(fp,"%d",&ti);
      if(qd->bd.ed[i][1]!=ti){
        printf("node id miss matched. check element id %d. Exit...\n",i*2);
        exit(1);
      }
    }
    // exchange open domain to main domain
    if(qd->bd.sd[i]==0){
      qd->bd.sd[i]=qd->bd.md[i];
      qd->bd.md[i]=0;
      if(etype==ELT3){
        tmpi=qd->bd.ed[i][1];
        qd->bd.ed[i][1]=qd->bd.ed[i][2];
        qd->bd.ed[i][2]=tmpi;
      }
      else {
        tmpi=qd->bd.ed[i][3];
        qd->bd.ed[i][3]=qd->bd.ed[i][1];
        qd->bd.ed[i][1]=tmpi;
      }
    }
  }

  fclose(fp);
  qd->bd.NN=nc;
}

void print_mesh_data(AQD2 *qd)
{
  printf("--- mesh data ---\n");
  printf("mesh data file name        = %s\n",qd->msh_fn);
  printf("node number                = %8d\n",qd->bd.Nn);
  printf("defined element number     = %8d\n",qd->bd.Ne*2);
  printf("\n");
}

void malloc_node(BOUD *bd)
{
  int i,N=bd->Nn;

  bd->rn=(double **)m_alloc2(N+1,sizeof(double*),"bem3_aw_qd2.c, malloc_node(), bd->rn");
  for(i=0;i<=N;i++){
    bd->rn[i]=(double *)m_alloc2(3,sizeof(double),"bem3_aw_qd2.c, malloc_node(), bd->rn[i]");
  }
}

void mfree_node(BOUD *bd)
{
  int i,N=bd->Nn;

  for(i=0;i<=N;i++) free(bd->rn[i]);
  free(bd->rn);
  bd->Nn=0;
}

void malloc_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  bd->ed =(int **)m_alloc2(Ne+1,sizeof(int *),"bem3_aw_qd2.c, malloc_elem(), bd->ed");
  bd->eni=(int **)m_alloc2(Ne+1,sizeof(int *),"bem3_aw_qd2.c, malloc_elem(), bd->eni");
  for(i=0;i<=Ne;i++){
    bd->ed [i]=(int *)m_alloc2(4,sizeof(int ),"bem3_aw_qd2.c, malloc_elem(), bd->ed[i]");
    bd->eni[i]=(int *)m_alloc2(4,sizeof(int ),"bem3_aw_qd2.c, malloc_elem(), bd->eni[i]");
  }

  bd->md=(int *)m_alloc2(Ne+1,sizeof(int),"bem3_aw_qd2.c, malloc_elem(), bd->md");
  bd->sd=(int *)m_alloc2(Ne+1,sizeof(int),"bem3_aw_qd2.c, malloc_elem(), bd->sd");
  bd->gd=(int *)m_alloc2(Ne+1,sizeof(int),"bem3_aw_qd2.c, malloc_elem(), bd->gd");

  // element constant
  bd->cr =(double ***)m_alloc2(Ne+1,sizeof(double **),"bem3_aw_qd2.c, malloc_elem(), bd->cr");
  bd->cw =(double ***)m_alloc2(Ne+1,sizeof(double **),"bem3_aw_qd2.c, malloc_elem(), bd->cw");
  bd->ren=(double ***)m_alloc2(Ne+1,sizeof(double **),"bem3_aw_qd2.c, malloc_elem(), bd->ren");
  bd->wen=(double ***)m_alloc2(Ne+1,sizeof(double **),"bem3_aw_qd2.c, malloc_elem(), bd->wen");
  bd-> Pi=(double complex **)m_alloc2(Ne+1,sizeof(double complex *),"bem3_aw_qd2.c, malloc_elem(), bd->Pi");
  bd->dPi=(double complex **)m_alloc2(Ne+1,sizeof(double complex *),"bem3_aw_qd2.c, malloc_elem(), bd->dPi");
  for(i=0;i<=Ne;i++){
    bd->cr[i]=(double **)m_alloc2(3,sizeof(double *),"bem3_aw_qd2.c, malloc_elem(), bd->cr[i]");
    bd->cw[i]=(double **)m_alloc2(3,sizeof(double *),"bem3_aw_qd2.c, malloc_elem(), bd->cw[i]");
    for(j=0;j<3;j++){
      bd->cr[i][j]=(double *)m_alloc2(4,sizeof(double),"bem3_aw_qd2.c, malloc_elem(), bd->cr[i][j]");
      bd->cw[i][j]=(double *)m_alloc2(3,sizeof(double),"bem3_aw_qd2.c, malloc_elem(), bd->cw[i][j]");
    }

    bd->ren[i]=(double **)m_alloc2(4,sizeof(double *),"bem3_aw_qd2.c, malloc_elem(), bd->ren[i]");
    bd->wen[i]=(double **)m_alloc2(4,sizeof(double *),"bem3_aw_qd2.c, malloc_elem(), bd->wen[i]");
    bd-> Pi[i]=(double complex *)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_elem(), bd->Pi[i]");
    bd->dPi[i]=(double complex *)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_elem(), bd->dPi[i]");
    for(j=0;j<4;j++){
      bd->ren[i][j]=(double *)m_alloc2(3,sizeof(double),"bem3_aw_qd2.c, malloc_elem(), bd->ren[i][j]");
      bd->wen[i][j]=(double *)m_alloc2(3,sizeof(double),"bem3_aw_qd2.c, malloc_elem(), bd->wen[i][j]");
    }
  }
}

void mfree_elem(BOUD *bd)
{
  int i,j,Ne=bd->Ne;

  for(i=0;i<=Ne;i++){
    free(bd->ed[i]);
    free(bd->eni[i]);
  }
  free(bd->ed);
  free(bd->eni);

  free(bd->md);
  free(bd->sd);
  free(bd->gd);

  for(i=0;i<=Ne;i++){
    for(j=0;j<3;j++){
      free(bd->cr[i][j]);
      free(bd->cw[i][j]);
    }
    free(bd->cr[i]);
    free(bd->cw[i]);

    for(j=0;j<4;j++){
      free(bd->ren[i][j]);
      free(bd->wen[i][j]);
    }
    free(bd->ren[i]);
    free(bd->wen[i]);
    free(bd->Pi[i]);
    free(bd->dPi[i]);
  }
  free(bd->cr);
  free(bd->cw);
  free(bd->ren);
  free(bd->wen);
  free(bd->Pi);
  free(bd->dPi);

  bd->Ne=0;
}

void init_elem_const(BOUD *bd)
{
  int i,j,d,Ne,a,b;
  double rc[3][4];

  Ne=bd->Ne;

  // geometric constant
  for(i=1;i<=Ne;i++){
    for(d=0;d<3;d++)
      for(j=0;j<4;j++) rc[d][j]=bd->rn[bd->ed[i][j]][d];

    if(bd->ed[i][3]!=0){ // bi-linear element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=0.25*( rc[d][0]+rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][1]=0.25*(-rc[d][0]+rc[d][1]+rc[d][2]-rc[d][3]);
        bd->cr[i][d][2]=0.25*(-rc[d][0]-rc[d][1]+rc[d][2]+rc[d][3]);
        bd->cr[i][d][3]=0.25*( rc[d][0]-rc[d][1]+rc[d][2]-rc[d][3]);
      }

      a=1;      b=2;
      bd->cw[i][0][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][0][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][0][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=2;      b=0;
      bd->cw[i][1][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][1][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][1][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
      a=0;      b=1;
      bd->cw[i][2][0]=0.125*( (rc[a][0]-rc[a][2])*(rc[b][1]-rc[b][3]) - (rc[a][1]-rc[a][3])*(rc[b][0]-rc[b][2]) );
      bd->cw[i][2][1]=0.125*( (rc[a][2]-rc[a][3])*(rc[b][0]-rc[b][1]) - (rc[a][0]-rc[a][1])*(rc[b][2]-rc[b][3]) );
      bd->cw[i][2][2]=0.125*( (rc[a][1]-rc[a][2])*(rc[b][0]-rc[b][3]) - (rc[a][0]-rc[a][3])*(rc[b][1]-rc[b][2]) );
    }
    else { // linear triangular element
      for(d=0;d<3;d++){
        bd->cr[i][d][0]=1.0/3.0*( rc[d][0]+rc[d][1]+rc[d][2]);
        bd->cr[i][d][1]=1.0/3.0*(-rc[d][0]+2.0*rc[d][1]-rc[d][2]);
        bd->cr[i][d][2]=1.0/sqrt(3.0)*( rc[d][2]-rc[d][0]);
        bd->cr[i][d][3]=0.0;
      }
      // cw : V
      bd->cw[i][0][0]=( (rc[1][1]-rc[1][0])*(rc[2][2]-rc[2][0]) - (rc[2][1]-rc[2][0])*(rc[1][2]-rc[1][0]) );
      bd->cw[i][0][1]=0.0;
      bd->cw[i][0][2]=0.0;
      bd->cw[i][1][0]=( (rc[2][1]-rc[2][0])*(rc[0][2]-rc[0][0]) - (rc[0][1]-rc[0][0])*(rc[2][2]-rc[2][0]) );
      bd->cw[i][1][1]=0.0;
      bd->cw[i][1][2]=0.0;
      bd->cw[i][2][0]=( (rc[0][1]-rc[0][0])*(rc[1][2]-rc[1][0]) - (rc[1][1]-rc[1][0])*(rc[0][2]-rc[0][0]) );
      bd->cw[i][2][1]=0.0;
      bd->cw[i][2][2]=0.0;
    }
  }
  // element constant
  // gaussian quadrature node and weight
  bd->zt_44[0]=-P44_N;  bd->zt_44[1]= P44_N;  bd->zt_44[2]= P44_N;  bd->zt_44[3]=-P44_N;
  bd->et_44[0]=-P44_N;  bd->et_44[1]=-P44_N;  bd->et_44[2]= P44_N;  bd->et_44[3]= P44_N;
  bd->wt_44[0]= P44_W;  bd->wt_44[1]= P44_W;  bd->wt_44[2]= P44_W;  bd->wt_44[3]= P44_W;

  bd->zt_49[0]=-P49_N;   bd->zt_49[1]= P49_N;   bd->zt_49[2]= P49_N;
  bd->zt_49[3]=-P49_N;   bd->zt_49[4]= 0.0;     bd->zt_49[5]= P49_N;
  bd->zt_49[6]= 0.0;     bd->zt_49[7]=-P49_N;   bd->zt_49[8]= 0.0;
  bd->et_49[0]=-P49_N;   bd->et_49[1]=-P49_N;   bd->et_49[2]= P49_N;
  bd->et_49[3]= P49_N;   bd->et_49[4]=-P49_N;   bd->et_49[5]= 0.0;
  bd->et_49[6]= P49_N;   bd->et_49[7]= 0.0;     bd->et_49[8]= 0.0;
  bd->wt_49[0]= P49_W0;  bd->wt_49[1]= P49_W0;  bd->wt_49[2]= P49_W0;
  bd->wt_49[3]= P49_W0;  bd->wt_49[4]= P49_W1;  bd->wt_49[5]= P49_W1;
  bd->wt_49[6]= P49_W1;  bd->wt_49[7]= P49_W1;  bd->wt_49[8]= P49_W2;

  bd->zt_34[0]=-P34_N0;  bd->zt_34[1]= 2.0*P34_N0;   bd->zt_34[2]=-P34_N0;  bd->zt_34[3]= 0.0;
  bd->et_34[0]=-P34_N1;  bd->et_34[1]= 0.0;          bd->et_34[2]= P34_N1;  bd->et_34[3]= 0.0;
  bd->wt_34[0]= P34_W0;  bd->wt_34[1]= P34_W0;                   bd->wt_34[2]= P34_W0;  bd->wt_34[3]=-P34_W1;

  bd->zt_37[0]=-P37_N0;  bd->zt_37[1]= 2.0*P37_N0;  bd->zt_37[2]=-P37_N0;
  bd->zt_37[3]= P37_N1;  bd->zt_37[4]=-2.0*P37_N1;  bd->zt_37[5]= P37_N1;  bd->zt_37[6]= 0.0;
  bd->et_37[0]=-P37_N2;  bd->et_37[1]= 0.0;                     bd->et_37[2]= P37_N2;
  bd->et_37[3]= P37_N3;  bd->et_37[4]= 0.0;                     bd->et_37[5]=-P37_N3;  bd->et_37[6]= 0.0;
  bd->wt_37[0]= P37_W0;     bd->wt_37[1]= P37_W0;  bd->wt_37[2]= P37_W0;
  bd->wt_37[3]= P37_W1;     bd->wt_37[4]= P37_W1;  bd->wt_37[5]= P37_W1;     bd->wt_37[6]= P37_W2;

  // gauss-legendre GLN point rule
  gauleg(-1.0,1.0,bd->xli,bd->wli,GLN);
  gauleg(-1.0,1.0,bd->xhi,bd->whi,GHN);
}

void malloc_sub_domain(AQD2 *qd)
{
  int *Nc,i,j,k;
  Nc=(int *)m_alloc2(qd->MN+1,sizeof(int),"bem3_aw_qd2.c, malloc_sub_domain(), Nc");
  for(i=0;i<=qd->MN;i++) Nc[i]=0;

  for(i=1;i<=qd->bd.Ne;i++){
    Nc[qd->bd.md[i]]++;
    Nc[qd->bd.sd[i]]++;
  }
  //for(i=0;i<md->MN+1;i++)    printf("id %d Nc:%d\n",i,Nc[i]); // test

  qd->bd.sb=(SUBD *)m_alloc2(qd->MN+1,sizeof(SUBD),"bem3_aw_qd2.c, malloc_sub_domain(), ad->bd.sb");
  for(i=0;i<=qd->MN;i++){
    qd->bd.sb[i].Ne=Nc[i];
    qd->bd.sb[i].sid=(int *)m_alloc2(Nc[i]+1,sizeof(int),"bem3_aw_qd2.c, malloc_sub_domain(), ad->bd.sb[i],sid");
    qd->bd.sb[i].P =(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].P");
    qd->bd.sb[i].dP=(double complex **)m_alloc2(Nc[i]+1,sizeof(double complex *),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].dP");
    qd->bd.sb[i].pv =(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].pv");
    qd->bd.sb[i].dpv=(double complex ***)m_alloc2(Nc[i]+1,sizeof(double complex **),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].dpv");
    for(j=0;j<=Nc[i];j++){
      qd->bd.sb[i].P [j]=(double complex *)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].P[j]");
      qd->bd.sb[i].dP[j]=(double complex *)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].dP[j]");
      qd->bd.sb[i].pv [j]=(double complex **)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].pv[j]");
      qd->bd.sb[i].dpv[j]=(double complex **)m_alloc2(4,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].dpv[j]");
      for(k=0;k<4;k++){
        qd->bd.sb[i].pv [j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].pv[j][k]");
        qd->bd.sb[i].dpv[j][k]=(double complex *)m_alloc2(3,sizeof(double complex),"bem3_aw_qd2.c, malloc_sub_domain(), qd->bd.sb[i].dpv[j][k]");
      }
    }
  }

  free(Nc);
}

void mfree_sub_domain(AQD2 *qd)
{
  int i,j,k;

  for(i=0;i<=qd->MN;i++){
    for(j=0;j<=qd->bd.sb[i].Ne;j++){
      free(qd->bd.sb[i].P [j]);
      free(qd->bd.sb[i].dP[j]);
      for(k=0;k<4;k++){
        free(qd->bd.sb[i].pv [j][k]);
        free(qd->bd.sb[i].dpv[j][k]);
      }
      free(qd->bd.sb[i].pv [j]);
      free(qd->bd.sb[i].dpv[j]);
    }
    free(qd->bd.sb[i].sid);
    free(qd->bd.sb[i].P);
    free(qd->bd.sb[i].dP);
    free(qd->bd.sb[i].pv);
    free(qd->bd.sb[i].dpv);
  }

  free(qd->bd.sb);
}

void init_sub_domain(AQD2 *qd)
{
  int d,i,c;

  for(d=0;d<=qd->MN;d++){
    c=1;
    for(i=1;i<=qd->bd.Ne;i++){
      if(qd->bd.md[i]==d){
        qd->bd.sb[d].sid[c]=i;
        c++;
      }
      else if(qd->bd.sd[i]==d){
        qd->bd.sb[d].sid[c]=-i;
        c++;
      }
    }
  }
}

void init_boundary_data(AQD2 *qd)
{
  double complex p,dp,cp;
  double cr[3][4],cw[3][3],r[3],w[3];
  int i,j,l,m;

  cp=-1.0/qd->pw.wd.k2;

  // element node and incident field data
  for(i=1;i<=qd->bd.Ne;i++){
    for(l=0;l<3;l++){
      for(m=0;m<4;m++) cr[l][m]=qd->bd.cr[i][l][m];
      for(m=0;m<3;m++) cw[l][m]=qd->bd.cw[i][l][m];
    }

    if(qd->bd.ed[i][3]==0){ // linear triangular element
      for(j=0;j<4;j++){
        lit_rw_zeta_eta(r,w,qd->bd.zt_34[j],qd->bd.et_34[j],cr,cw);
        for(l=0;l<3;l++){
          qd->bd.ren[i][j][l]=r[l];
          qd->bd.wen[i][j][l]=w[l];
        }
        vuni_d(w);
        calc_apw_dpdn(&p,&dp,r,w,&(qd->pw)); // aw_pw.h
        qd->bd. Pi[i][j]=cp* p;
        qd->bd.dPi[i][j]=cp*dp;
      }
    }
    else { // bi-linear element
      for(j=0;j<4;j++){
        bil_rw_zeta_eta(r,w,qd->bd.zt_44[j],qd->bd.et_44[j],cr,cw);
        for(l=0;l<3;l++){
          qd->bd.ren[i][j][l]=r[l];
          qd->bd.wen[i][j][l]=w[l];
        }
        vuni_d(w);
        calc_apw_dpdn(&p,&dp,r,w,&(qd->pw)); // aw_pw.h
        qd->bd. Pi[i][j]=cp* p;
        qd->bd.dPi[i][j]=cp*dp;
      }
    }
  }
}

int domain_id(double *rt,AQD2 *qd)
{
  double fid_calc_solid_angle(int type,double r[4][3],int *flg);

  double rv[4][3],omega;
  double *Og=(double *)m_alloc2(qd->MN+1,sizeof(double),"bem3_aw_qd1.c, domain_id()");
  int i,j,k,d,flg;

  for(d=0;d<qd->MN+1;d++){
    for(i=1;i<=qd->bd.sb[d].Ne;i++){
      if(qd->bd.sb[d].sid[i]>0){
        // read node data
        for(j=0;j<4;j++)
          for(k=0;k<3;k++) rv[j][k]=qd->bd.rn[qd->bd.ed[qd->bd.sb[d].sid[i]][j]][k]-rt[k];
        omega=fid_calc_solid_angle(check_element_type(qd->bd.sb[d].sid[i],&(qd->bd)),rv,&flg);
        if(flg<0){ // on boundary
          free(Og);
          return d;
        }
        Og[d]+=omega;
        Og[qd->bd.sd[qd->bd.sb[d].sid[i]]]-=omega;

      } // end if
    }
    if(d==0 && fabs(Og[d])<2.0*M_PI){ // opened region
      free(Og);
      return d;
    }
    else if(Og[d]>2.0*M_PI){ // closed region
      free(Og);
      return d;
    }
  }

  free(Og);
  return -1; // error
}

double fid_calc_solid_angle(int type,double r[4][3],int *flg)
{
  double n01[3],n12[3],n20[3],nt[3],sa,ca;
  double Omega,a0,a1,a2,D;
  int st,i;

  Omega=0.0;
  *flg=0;

  st=0;
  if(type==ELT3) for(i=0;i<3;i++)  st+=vuni_d(r[i]);
  else for(i=0;i<4;i++) st+=vuni_d(r[i]);
  if(st<0){
    *flg=-1;    return 0; // on boundary
  }

  vcrs_d(n01,r[0],r[1]);
  vcrs_d(n12,r[1],r[2]);
  vcrs_d(n20,r[2],r[0]);

  vcrs_d(nt,n01,n20);
  sa=vabs_d(nt);
  ca=-vdot_d(n01,n20);
  a0=atan2(sa,ca);

  vcrs_d(nt,n12,n01);
  sa=vabs_d(nt);
  ca=-vdot_d(n12,n01);
  a1=atan2(sa,ca);

  vcrs_d(nt,n20,n12);
  sa=vabs_d(nt);
  ca=-vdot_d(n20,n12);
  a2=atan2(sa,ca);

  D=vdot_d(r[0],n12);
  if(D>0.0) Omega+=a0+a1+a2-M_PI;
  else if(D<0.0) Omega-=a0+a1+a2-M_PI;
  else { // on boundary
    *flg=-1;    return 0;
  }

  if(ELT4==type){
    n01[0]=-n20[0];    n01[1]=-n20[1];    n01[2]=-n20[2];

    vcrs_d(n12,r[2],r[3]);
    vcrs_d(n20,r[3],r[0]);

    vcrs_d(nt,n01,n20);
    sa=vabs_d(nt);
    ca=-vdot_d(n01,n20);
    a0=atan2(sa,ca);

    vcrs_d(nt,n12,n01);
    sa=vabs_d(nt);
    ca=-vdot_d(n12,n01);
    a1=atan2(sa,ca);

    vcrs_d(nt,n20,n12);
    sa=vabs_d(nt);
    ca=-vdot_d(n20,n12);
    a2=atan2(sa,ca);

    D=vdot_d(r[0],n12);
    if(D>0.0) Omega+=a0+a1+a2-M_PI;
    else if(D<0.0) Omega-=a0+a1+a2-M_PI;
    else { // on boundary
      *flg=-1;    return 0;
    }
  }
  return Omega;
}

