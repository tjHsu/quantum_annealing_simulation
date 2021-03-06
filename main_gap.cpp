#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <complex>
#include <cmath>
#include <time.h>
#include "spin_system.h"

using namespace std;
// #define MKL_Complex16 complex<double>
// #include "omp.h"
// #include "mkl.h"


int main(int argc, char* argv[]){

  spin_system test;
  int N_sys=8;
  int N_env=0;
  //double Time=atof(argv[1]);
  double J_start=atoi(argv[1]);
  double tau=atof(argv[2]);
  double Temperature=atof(argv[3]);
  //double G=atof(argv[4]);
  // test.initialize(N_sys,N_env,Time,tau,Temperature);

  int N=N_sys+N_env;
  int nofstates=(int) pow(2,N);

  ofstream success_probability_out("output_ProbabilityVsT.dat");

  // double T[10]={1e1,1e2,2e2,5e2,1e3,2e3,5e3,1e4,2e4,1e5};
  double T[40]={1e1,1.25e1,1.5e1,1.75e1,1.875e1,2e1,2.25e1,2.5e1,2.75e1,2.875e1,3e1,3.25e1,3.5e1,4e1,5e1,6e1,7e1,8e1,1e2,1.25e2,1.5e2,1.75e2,2e2,2.25e2,2.5e2,2.75e2,3e2,3.25e2,3.5e2,4e2,5e2,6e2,7e2,8e2,1e3,1.5e3,2e3,3e3,5e3,1e4};
  // double J[6]={0,0.05,0.1,0.2,0.5,1};
  success_probability_out<<"Total_steps(tau=0.1) ";
  // for (int i = J_start; i < J_start+1; i++) {
  //   success_probability_out<<"J="<<J[i]<<" ";
  // }
  // success_probability_out<<endl;

  time_t start_all=time(0);
  clock_t t_all;
  t_all = clock();
  int env_on=1;
  double min_tmp=100;
  double max_tmp=0;
  double gap_tmp=0;
  double T_tmp=0;
  for (int i = 0; i < 40; i++) {
    success_probability_out<<T[i]<<" ";
    for (int j = 1; j < 101; j++) {
      cout<<"Run with Time_steps= "<<T[i]<<", J_index= "<<j<<"."<<endl;

      test.initialize(N_sys,N_env,(T[i]/10),tau,Temperature,0, env_on, j);

      if (1==env_on) {
        env_on=0;
      }
      test.skip_zeroterm();

      time_t start=time(0);
      clock_t t;
      t = clock();
      test.run();
      t =clock()-t;
      time_t end=time(0);
      double time=difftime(end,start);
      double norm=0.;
      double sum_GS=0.;
      double sum_ES=0.;
      for(int i=0; i<nofstates;i++){
        norm+=test.coefficient_return[i];
        if (0==(i-test.gs_sol)%256) {
          sum_GS+=test.coefficient_return[i];
        } else {
          sum_ES+=test.coefficient_return[i];
        }
      }


      cout<<"This Run took me "<<t<<" clicks ( "<<((float) t)/CLOCKS_PER_SEC<<" processing seconds.)"<<endl;
      cout<<"costed me "<<time<<" second in real time"<<endl;
      cout<<"norm after run = "<<norm<<endl;
      cout<<"pobability of GS: "<<sum_GS<<endl;
      cout<<"pobability of non-GS: "<<sum_ES<<endl;
      cout<<endl;


      success_probability_out<<test.success_probability_return<<" ";
      if (test.success_probability_return>max_tmp)
        max_tmp=test.success_probability_return;
      if (test.success_probability_return<min_tmp)
        min_tmp=test.success_probability_return;

    }
    if (max_tmp-min_tmp>gap_tmp){
      gap_tmp=max_tmp-min_tmp;
      T_tmp=T[i];
    }
    max_tmp=0;
    min_tmp=100;
    success_probability_out<<endl;
  }

  t_all =clock()-t_all;
  time_t end_all=time(0);
  double time_all=difftime(end_all,start_all);
  cout<<"My program overall took me "<<t_all<<" clicks ( "<<((float) t_all)/CLOCKS_PER_SEC<<" processing seconds.)"<<endl;
  cout<<"costed me "<<time_all<<" second in real time"<<endl;
  cout<<"T_tmp= "<<T_tmp<<", Gap_tmp= "<<gap_tmp<<endl;





  cout<<"reutrn "<<0<<endl;
  return 0;
}
