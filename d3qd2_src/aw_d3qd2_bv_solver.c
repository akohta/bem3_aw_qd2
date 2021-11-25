#include "bem3_aw_qd2.h"

int main(int argc,char **argv)
{
  AQD2 qd;
  
  read_aqd2(argc,argv,&qd);
  print_aqd2(&qd);
  initialize_aqd2(&qd);
  solve_bieq_aqd2(&qd);
  dat_write_aqd2(argv[5],&qd); 
  finalize_aqd2(&qd);
  return 0;
}
