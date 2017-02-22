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

using namespace std;
#define MKL_Complex16 complex<double>
#include "omp.h"
#include "mkl.h"

class spin_system {
private:
  int N;
  int N_sys;
  int N_env;
  double T;
  double tau;
  int nofstates;
  int nofstates_sys;
  int nofstates_env;

  double* ss_operator_real;
  double* ss_operator_imaginary;
  double* ds_operator_real;
  double* ds_operator_imaginary;
  double* J_x;
  double* J_y;
  double* J_z;
  double G;
  double* Jx_env;
  double* Jy_env;
  double* Jz_env;
  double* Jx_se;
  double* Jy_se;
  double* Jz_se;
  double* h_x;
  double* h_y;
  double* h_z;
  double  h_x_start;
  double Gamma;
  double Delta;

  /*Basis for Heat Bath Environment*/
  complex<double>* z;
  double* w;

  double* psi_real;
  double* psi_imaginary;
  double* psi_tmp_real;
  double* psi_tmp_imaginary;
  double* psi_sys_real;
  double* psi_sys_imaginary;

  // /*Not yet finished: optimizing by skipping the zero term.*/
  // double* Jz_k_marked;
  // double* Jz_l_marked;
  // int count_z;
  // double* Jy_k_marked;
  // double* Jy_l_marked;
  // int count_y;
  // double* Jx_k_marked;
  // double* Jx_l_marked;
  // int count_x;

  // /*uncommend if want to use spin_allinone();*/
  // double* psi_tmp_x_real;
  // double* psi_tmp_x_imaginary;
  // double* psi_tmp_y_real;
  // double* psi_tmp_y_imaginary;
  // double* psi_tmp_z_real;
  // double* psi_tmp_z_imaginary;

  void single_spin_op(double);
  void double_spin_op_x(double);
  void double_spin_op_y(double);
  void double_spin_op_z(double);
  void set_initial_sys_state(char const *);
  double energy(double);
  double energy_all(double);
  double energy_env(double);
  double energy_se(double);
  double spin(char,int);
  // void spin_allinone();
  void environment(int, double);
  void Jenv_generate(int,double);
  void Jse_generate(int, int, double);
  void generate(int, double*, double*, char const *, char const *, char const *, char const *);
  void direct_product(int, double*, double*, complex<double>*, double*, double*);

  void read(int,double*,char const *);



public:
  void initialize(int, int ,double ,double, double );
  void run();

  double* spin_return;
  double* energy_sys_return;
  double* energy_env_return;
  double* energy_se_return;
  double* energy_all_return;
  double* coefficient_return;
};

int main(int argc, char* argv[]){

  spin_system test;
  int N_sys=8;
  int N_env=8;
  double Time=atof(argv[1]);
  double tau=atof(argv[2]);
  double Temperature=atof(argv[3]);
  test.initialize(N_sys,N_env,Time,tau,Temperature);

  int N=N_sys+N_env;
  int nofstates=(int) pow(2,N);

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
    if (0==(i-119)%256) {
      sum_GS+=test.coefficient_return[i];
    } else {
      sum_ES+=test.coefficient_return[i];
    }
  }

  cout<<endl;
  cout<<"It took me "<<t<<" clicks ( "<<((float) t)/CLOCKS_PER_SEC<<" processing seconds.)"<<endl;
  cout<<"costed me "<<time<<" second in real time"<<endl;
  cout<<"norm after run = "<<norm<<endl;
  cout<<"pobability of GS: "<<sum_GS<<endl;
  cout<<"pobability of non-GS: "<<sum_ES<<endl;
  cout<<endl;
  cout<<"reutrn "<<0<<endl;
  return 0;
}


/* Initialize the annealing system with environment.
  Input:
    N_sys_user_defined: number of spins of system
    N_env_user_defined: number of spins of environment
    T_user_defined: The Annealing time
    tau_user_defined: the timestep of simulation
*/
void spin_system::initialize(int N_sys_user_defined, int N_env_user_defined, double T_user_defined, double tau_user_defined, double Temperature_user_defined){
  double Temperature = Temperature_user_defined;
  N_sys = N_sys_user_defined;
  N_env = N_env_user_defined;
  N = N_sys+N_env;
  T = T_user_defined;
  tau = tau_user_defined;
  nofstates = (int) pow(2,N);// number of states construct by the number of total spins

  /* initialize the coupling factor J for double spin Hamiltonian,
    and the factor h for single spin Hamiltonian.
    J is for system itself, J_env is for environment itself, and J_se is the interaction
  */
  J_x = new double [N_sys*N_sys]();
  J_y = new double [N_sys*N_sys]();
  J_z = new double [N_sys*N_sys]();
  Jx_env = new double [N_env*N_env]();
  Jy_env = new double [N_env*N_env]();
  Jz_env = new double [N_env*N_env]();
  Jx_se = new double [N_env*N_sys]();
  Jy_se = new double [N_env*N_sys]();
  Jz_se = new double [N_env*N_sys]();
  h_x = new double [N]();
  h_y = new double [N]();
  h_z = new double [N]();
  h_x_start= 1;
  // for (int i = 0; i < N_sys*N_sys; i++){
  //   J_x[i] = 0.;
  //   J_y[i] = 0.;
  //   J_z[i] = 0.;
  // }
  // for (int i = 0; i < N_env*N_env; i++){
  //   Jx_env[i] = 0.;
  //   Jy_env[i] = 0.;
  //   Jz_env[i] = 0.;
  // }
  // for (int i = 0; i < N_env*N_sys; i++){
  //   Jx_se[i] = 0.;
  //   Jy_se[i] = 0.;
  //   Jz_se[i] = 0.;
  // }
  // for (int i = 0; i < N; i++){
  //   h_x[i]=0.;
  //   h_y[i]=0.;
  //   h_z[i]=0.;
  // }


  read(N_sys*N_sys,J_z,"J4.txt");
  read(N_sys*N_sys,J_x,"J4x.txt");
  read(N_sys*N_sys,J_y,"J4y.txt");
  read(N_env*N_env,Jx_env,"J4x_env.txt");
  read(N_env*N_env,Jy_env,"J4y_env.txt");
  read(N_env*N_env,Jz_env,"J4z_env.txt");
  read(N_sys*N_env,Jx_se,"J4x_se.txt");
  read(N_sys*N_env,Jy_se,"J4y_se.txt");
  read(N_sys*N_env,Jz_se,"J4z_se.txt");
  read(N,h_z,"h4.txt");
  read(N,h_x,"h4x.txt");
  read(N,h_y,"h4y.txt");


  /* initialize the array for the enivironment's partition factor w[], and its eignevector in computational basis z[] */
  z= new complex<double> [(int)pow(2,N_env)*(int)pow(2,N_env)]();
  // for (int i = 0; i < (int)pow(2,N_env)*(int)pow(2,N_env); i++) {
  //   z[i]=0;
  // }
  w= new double [(int)pow(2,N_env)]();
  // for (int i = 0; i < (int)pow(2,N_env); i++) {
  //   w[i]=0;
  // }
  /*G is the global factor usually set between -1 and 1.*/
  // G=1.0;
  // Jenv_generate(N_env,G);//randomly generate J_env
  // G=0.05;
  // Jse_generate(N_sys,N_env,G);//randomly generate J_se
  /*the second parameter is Temperature*/
  environment(N_env,Temperature);//get w[] and z[] with lapack diagonalization.




  /* initialize the wave function in the ground state
  */

  psi_real = new double [nofstates]();
  psi_imaginary = new double [nofstates]();
  // for (int i = 0; i < nofstates; i++) {
  //     psi_real[i]      = 0;
  //     psi_imaginary[i] = 0;
  // }
  psi_tmp_real = new double [nofstates]();
  psi_tmp_imaginary = new double [nofstates]();
  // for (int i = 0; i < nofstates; i++) {
  //   psi_tmp_real[i]      = 0;
  //   psi_tmp_imaginary[i] = 0;
  // }

  psi_sys_real = new double [(int)pow(2,N_sys)]();
  psi_sys_imaginary = new double [(int)pow(2,N_sys)]();
  // for (int i = 0; i < (int)pow(2,N_sys); i++) {
  //   psi_sys_real[i]      = 0;
  //   psi_sys_imaginary[i] = 0;
  // }
  set_initial_sys_state("allx");

  /* initialize the  matrix. We have two kind of Hamiltonian operator here.
    First is the single spin operator, ss_operator.
    Second is the double spin operator, ds_operator.
    Ref. Article: Computational Methods for Simulating Quantum Computers euqation (68) & equation (70).
  */
  ss_operator_real      = new double [4]();
  ss_operator_imaginary = new double [4]();
  ds_operator_real      = new double [8]();
  ds_operator_imaginary = new double [8]();
  // for (int i = 0; i < 4; i++) {
  //   ss_operator_real[i]      = 0;
  //   ss_operator_imaginary[i] = 0;
  // }
  // for (int i = 0; i < 8; i++) {
  //   ds_operator_real[i]      = 0;
  //   ds_operator_imaginary[i] = 0;
  // }

  // /* uncommend if want to use spin_allinone(); */
  // psi_tmp_x_real=new double [nofstates];
  // psi_tmp_x_imaginary=new double [nofstates];
  // psi_tmp_y_real=new double [nofstates];
  // psi_tmp_y_imaginary=new double [nofstates];
  // psi_tmp_z_real=new double [nofstates];
  // psi_tmp_z_imaginary=new double [nofstates];


  /*set the return measuremt*/
  int total_steps=(int) (T/tau);

  spin_return = new double [N*3*(total_steps+1)];//3 state for x,y,z
  for (int i = 0; i < N*3*(total_steps+1); i++){
    spin_return[i]=0.;
  }
  energy_sys_return = new double [total_steps+1];
  energy_env_return = new double [total_steps+1];
  energy_se_return  = new double [total_steps+1];
  energy_all_return = new double [total_steps+1];

  for (int i = 0; i < (total_steps+1); i++) {
    energy_sys_return[i]=0.;
    energy_env_return[i]=0.;
    energy_se_return[i]=0.;
    energy_all_return[i]=0.;
  }
  coefficient_return = new double [nofstates];
  for (int i = 0; i < nofstates; i++) {
    coefficient_return[i]=0.;
  }


  Gamma=1.; //time evolution of the initial Hamiltonian. Goes from 1 to 0
  Delta=1.; //time evolution of the desired Hamiltonian. Goes from 0 to 1

}


/* Set up the initial basis state
  Input:
    d: define which mode to set up
      "read": read "psi_real.dat" and "psi_imagine_dat"
      "1upRan": set spin_0 as up in x and others are random.
      "allx": as the ground state of sigma_x.
      "allRand": all random.
*/
void spin_system::set_initial_sys_state(char const * d){


  if ("read"==d){
    read((int)pow(2,N_sys),psi_sys_real,"psi_real.dat");
    read((int)pow(2,N_sys),psi_sys_imaginary,"psi_imagine.dat");
    cout<<"Read for set_initial_sys_state() function"<<endl;
  }
  else if ("1upRand"==d) {
    srand (time(NULL));
    double normalize_factor=0.;
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
      if(i%2==0){
        psi_sys_real[i]      = (double) rand()/RAND_MAX;//pow(nofstates, -0.5);
        psi_sys_imaginary[i] = (double) rand()/RAND_MAX;
        normalize_factor += psi_sys_real[i]*psi_sys_real[i] + psi_sys_imaginary[i]*psi_sys_imaginary[i];
      } else {
        psi_sys_real[i]      = 0;
        psi_sys_imaginary[i] = 0;
      }
    }
    ofstream Psi_r_out("psi_real.dat");
    ofstream Psi_i_out("psi_imagine.dat");
    normalize_factor = sqrt(normalize_factor);
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
      if(i%2==0){
        psi_sys_real[i]      =psi_sys_real[i]/normalize_factor;
        psi_sys_imaginary[i] =psi_sys_imaginary[i]/normalize_factor;
      }
      Psi_r_out<<psi_sys_real[i]<<endl;
      Psi_i_out<<psi_sys_imaginary[i]<<endl;
    }

  }
  else if ("allx"==d){
    ofstream Psi_r_out("psi_real.dat");
    ofstream Psi_i_out("psi_imagine.dat");
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
        psi_sys_real[i]      =pow((int) pow(2,N_sys), -0.5);
        psi_sys_imaginary[i] =0;
        Psi_r_out<<psi_sys_real[i]<<endl;
        Psi_i_out<<psi_sys_imaginary[i]<<endl;
    }
  }
  else if ("allRand"==d){
    srand (time(NULL));
    double normalize_factor=0.;
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
        psi_sys_real[i]      = (double) rand()/RAND_MAX;//pow(nofstates, -0.5);
        psi_sys_imaginary[i] = (double) rand()/RAND_MAX;
        normalize_factor += psi_sys_real[i]*psi_sys_real[i] + psi_sys_imaginary[i]*psi_sys_imaginary[i];
    }
    ofstream Psi_r_out("psi_real.dat");
    ofstream Psi_i_out("psi_imagine.dat");
    normalize_factor = sqrt(normalize_factor);
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
      psi_sys_real[i]      =psi_sys_real[i]/normalize_factor;
      psi_sys_imaginary[i] =psi_sys_imaginary[i]/normalize_factor;
      Psi_r_out<<psi_sys_real[i]<<endl;
      Psi_i_out<<psi_sys_imaginary[i]<<endl;
    }

  }

  else {
    cout<<endl<<"Wrong Parameter for function: generate_initial_state(char)"<<endl;
    cout<<"!!!BE AWARE of the correctness of the result!!!"<<endl<<endl;

  }

}

/* Operating sigma_x,_y,_z*h to the spins.
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[], ss_operator_real[], ss_operator_imaginary[]
  #:now ONLY deal with h in SYSTEM, not including Environment.

*/
void spin_system::single_spin_op(double t){

  for (int k = 0; k < N; k++) {
    double hx=0.;
    double hy=0.;
    double hz=0.;
    if(k>=N_sys){
      hx= h_x[k];
      hy= h_y[k];
      hz= h_z[k];
    } else {
      hx=h_x_start*Gamma+h_x[k]*Delta;
      hy=h_y[k]*Delta;
      hz=h_z[k]*Delta;
    }
    int i1=(int) pow(2,k);

    /* update the single spin Hamiltonian matrix with t.
      In principal, we can save some computing time here,
      because of reading convenience I didn't do so.
    */
    double norm=0;
    norm=sqrt(hx*hx+hy*hy+hz*hz);
    //set the initial transverse field
    // if (k<N_sys) {
    //   norm=sqrt((Gamma*h_x_start+Delta*h_x[k])*(Gamma*h_x_start+Delta*h_x[k])+Delta*h_y[k]*Delta*h_y[k]+Delta*h_z[k]*Delta*h_z[k]);
    // } else {
    //   norm=0.;
    // }



    if (norm-0<1e-14) {
      ss_operator_real[0]      = 1;//cos(tau*norm*0.5);
      ss_operator_real[1]      = 0;
      ss_operator_real[2]      = 0;
      ss_operator_real[3]      = 1;//cos(tau*norm*0.5);
      ss_operator_imaginary[0] = 0;
      ss_operator_imaginary[1] = 0;
      ss_operator_imaginary[2] = 0;
      ss_operator_imaginary[3] = 0;
    }

    else {
      double a=cos(tau*0.5*norm);
      double b=sin(tau*0.5*norm)/norm;
      ss_operator_real[0]      = a;//cos(tau*0.5*norm);
      ss_operator_real[1]      = hy*b;//sin(tau*0.5*norm)/norm;
      ss_operator_real[2]      = -1*hy*b;//sin(tau*0.5*norm)/norm;
      ss_operator_real[3]      = a;//cos(tau*0.5*norm);
      ss_operator_imaginary[0] = hz*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[1] = hx*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[2] = hx*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[3] = -1*hz*b;//sin(tau*0.5*norm)/norm;
    }
    #pragma omp parallel for default(none) shared(i1)
    for (int l = 0; l < nofstates; l+=2) {
      /* get index for the second place we need to operate ss_operator on.
        this is a copy from the previous code.
        need to be moderated more properly without divde operator.
      */
      int i2= l & i1;
      int i = l -i2+i2/i1;
      int j = i+i1;

      /* operate the matrix to the spin.
      */
      double psi_real_temp_i      = 0;
      double psi_imaginary_temp_i = 0;
      double psi_real_temp_j      = 0;
      double psi_imaginary_temp_j = 0;

      psi_real_temp_i = ss_operator_real[0]*psi_real[i] + ss_operator_real[1]*psi_real[j]
                      - ss_operator_imaginary[0]*psi_imaginary[i] - ss_operator_imaginary[1]*psi_imaginary[j];
      psi_imaginary_temp_i = ss_operator_real[0]*psi_imaginary[i] + ss_operator_real[1]*psi_imaginary[j]
                         + ss_operator_imaginary[0]*psi_real[i] + ss_operator_imaginary[1]*psi_real[j];

      psi_real_temp_j = ss_operator_real[2]*psi_real[i] + ss_operator_real[3]*psi_real[j]
                      - ss_operator_imaginary[2]*psi_imaginary[i] - ss_operator_imaginary[3]*psi_imaginary[j];
      psi_imaginary_temp_j = ss_operator_real[2]*psi_imaginary[i] + ss_operator_real[3]*psi_imaginary[j]
                         + ss_operator_imaginary[2]*psi_real[i] + ss_operator_imaginary[3]*psi_real[j];
      psi_real[i]       = psi_real_temp_i;
      psi_imaginary[i]  = psi_imaginary_temp_i;
      psi_real[j]       = psi_real_temp_j;
      psi_imaginary[j]  = psi_imaginary_temp_j;

    }
  }
}

/* Operating sigma_x*sigma_x
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[],ds_operator_real[],ds_operator_imaginary[]
*/
void spin_system::double_spin_op_x(double t){
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jx_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jx_se[k+(l-N_sys)*N_sys];
      } else {
        J=Delta*J_x[k+l*N_sys];
      }
  // for (int i = 0; i < count_x; i++) { //optimize use
    // int k=Jx_k_marked[i];
    // int l=Jx_l_marked[i];
      if(abs(J)>1e-15){

        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=0;//J_z[k+l*N]/1.;//4.;
        double b=J;//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=J;//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double b_block=b*tau*0.5;
        double c_block=c*tau*0.5;

        double cos_b_block=cos(b_block);
        double cos_c_block=cos(c_block);
        double sin_b_block=sin(b_block);
        double sin_c_block=sin(c_block);

        ds_operator_real[0]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[2]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[5]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[7]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[1] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[3] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[4] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[6] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj)
        for (int m = 0; m < nofstates; m+=4) {
          /* get index for the second place we need to operate ss_operator on.
            this is a copy from the previous code.
            need to be moderated more properly without divde operator.
          */
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;


          /* Following the similar manner in single_spin_op,
            I create several temperory variables to store values of psi
          */
          double psi_real_temp_n0      = 0;
          double psi_imaginary_temp_n0 = 0;
          double psi_real_temp_n1      = 0;
          double psi_imaginary_temp_n1 = 0;
          double psi_real_temp_n2      = 0;
          double psi_imaginary_temp_n2 = 0;
          double psi_real_temp_n3      = 0;
          double psi_imaginary_temp_n3 = 0;

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  /*-ds_operator_imaginary[0]*psi_imaginary[n0]*/ - ds_operator_imaginary[1]*psi_imaginary[n3];
          psi_imaginary_temp_n0 = /*ds_operator_imaginary[0]*psi_real[n0] +*/ ds_operator_imaginary[1]*psi_real[n3]
                                  +ds_operator_real[0]*psi_imaginary[n0] /*+ ds_operator_real[1]*psi_imaginary[n3]*/;

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] /*+ ds_operator_real[3]*psi_real[n2]*/
                                  /*-ds_operator_imaginary[2]*psi_imaginary[n1]*/ - ds_operator_imaginary[3]*psi_imaginary[n2];
          psi_imaginary_temp_n1 = /*ds_operator_imaginary[2]*psi_real[n1] +*/ ds_operator_imaginary[3]*psi_real[n2]
                                  +ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;

          psi_real_temp_n2      = /*ds_operator_real[4]*psi_real[n1] + */ds_operator_real[5]*psi_real[n2]
                                  -ds_operator_imaginary[4]*psi_imaginary[n1] /*- ds_operator_imaginary[5]*psi_imaginary[n2]*/;
          psi_imaginary_temp_n2 = ds_operator_imaginary[4]*psi_real[n1] /*+ ds_operator_imaginary[5]*psi_real[n2]*/
                                  /*+ds_operator_real[4]*psi_imaginary[n1] */+ ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  -ds_operator_imaginary[6]*psi_imaginary[n0] /*- ds_operator_imaginary[7]*psi_imaginary[n3]*/;
          psi_imaginary_temp_n3 = ds_operator_imaginary[6]*psi_real[n0] /*+ ds_operator_imaginary[7]*psi_real[n3]*/
                                  /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];

          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      }
    }
  }
}

/* Operating sigma_y*sigma_y
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[],ds_operator_real[],ds_operator_imaginary[]
*/
void spin_system::double_spin_op_y(double t){
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jy_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jy_se[k+(l-N_sys)*N_sys];
      } else {
        J=Delta*J_y[k+l*N_sys];
      }

  // for (int i = 0; i < count_y; i++) {
    // int k=Jy_k_marked[i];
    // int l=Jy_l_marked[i];
      if(abs(J)>1e-15){

        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=0;//J_z[k+l*N]/1.;//4.;
        double b=-J;//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=J;//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double b_block=b*tau*0.5;
        double c_block=c*tau*0.5;

        double cos_b_block=cos(b_block);
        double cos_c_block=cos(c_block);
        double sin_b_block=sin(b_block);
        double sin_c_block=sin(c_block);



        ds_operator_real[0]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[2]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[5]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[7]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[1] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[3] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[4] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[6] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj)
        for (int m = 0; m < nofstates; m+=4) {
          /* get index for the second place we need to operate ss_operator on.
            this is a copy from the previous code.
            need to be moderated more properly without divde operator.
          */
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;

          /* Following the similar manner in single_spin_op,
            I create several temperory variables to store values of psi
          */
          double psi_real_temp_n0      = 0;
          double psi_imaginary_temp_n0 = 0;
          double psi_real_temp_n1      = 0;
          double psi_imaginary_temp_n1 = 0;
          double psi_real_temp_n2      = 0;
          double psi_imaginary_temp_n2 = 0;
          double psi_real_temp_n3      = 0;
          double psi_imaginary_temp_n3 = 0;

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  /*-ds_operator_imaginary[0]*psi_imaginary[n0]*/ - ds_operator_imaginary[1]*psi_imaginary[n3];
          psi_imaginary_temp_n0 = /*ds_operator_imaginary[0]*psi_real[n0] +*/ ds_operator_imaginary[1]*psi_real[n3]
                                  +ds_operator_real[0]*psi_imaginary[n0] /*+ ds_operator_real[1]*psi_imaginary[n3]*/;

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] /*+ ds_operator_real[3]*psi_real[n2]*/
                                  /*-ds_operator_imaginary[2]*psi_imaginary[n1]*/ - ds_operator_imaginary[3]*psi_imaginary[n2];
          psi_imaginary_temp_n1 = /*ds_operator_imaginary[2]*psi_real[n1] +*/ ds_operator_imaginary[3]*psi_real[n2]
                                  +ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;

          psi_real_temp_n2      = /*ds_operator_real[4]*psi_real[n1] + */ds_operator_real[5]*psi_real[n2]
                                  -ds_operator_imaginary[4]*psi_imaginary[n1] /*- ds_operator_imaginary[5]*psi_imaginary[n2]*/;
          psi_imaginary_temp_n2 = ds_operator_imaginary[4]*psi_real[n1] /*+ ds_operator_imaginary[5]*psi_real[n2]*/
                                  /*+ds_operator_real[4]*psi_imaginary[n1] */+ ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  -ds_operator_imaginary[6]*psi_imaginary[n0] /*- ds_operator_imaginary[7]*psi_imaginary[n3]*/;
          psi_imaginary_temp_n3 = ds_operator_imaginary[6]*psi_real[n0] /*+ ds_operator_imaginary[7]*psi_real[n3]*/
                                  /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];

          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      }
    }
  }
}

/* Operating sigma_z*sigma_z
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[],ds_operator_real[],ds_operator_imaginary[]
*/
void spin_system::double_spin_op_z(double t){

  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jz_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jz_se[k+(l-N_sys)*N_sys];
      } else {
        J=Delta*J_z[k+l*N_sys];
      }
  // for (int i = 0; i < count_z; i++) {
    // int k=Jz_k_marked[i];
    // int l=Jz_l_marked[i];
      if(abs(J)>1e-15){


        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=J;//4.;
        double b=0;//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=0;//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double a_block=a*tau;
        double cos_a_block=cos(a_block);
        double sin_a_block=sin(a_block);
        ds_operator_real[0]      = cos_a_block;//cos( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_real[1]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_real[2]      = cos_a_block;//cos(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_real[3]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_real[4]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_real[5]      = cos_a_block;//cos(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_real[6]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_real[7]      = cos_a_block;//cos( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_imaginary[0] = sin_a_block;//sin( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_imaginary[1] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_imaginary[2] = -sin_a_block;//sin(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_imaginary[3] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_imaginary[4] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_imaginary[5] = -sin_a_block;//sin(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_imaginary[6] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_imaginary[7] = sin_a_block;//sin( a*Delta*tau)*cos(b*Delta*tau);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj)
        for (int m = 0; m < nofstates; m+=4) {


          /* get index for the second place we need to operate ss_operator on.
            this is a copy from the previous code.
            need to be moderated more properly without divde operator.
          */
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;

          /////////checking!!!////////////
          // int tid=omp_get_thread_num();
          // if(3==ch)
          // cout<<"tid = "<<tid<<", n0, n1, n2, n3 = "<<n0<<", "<<n1<<", "<<n2<<", "<<n3<<endl;
          /////////checking!!!////////////
          /* Following the similar manner in single_spin_op,

            I create several temperory variables to store values of psi
          */
          double psi_real_temp_n0      = 0;
          double psi_imaginary_temp_n0 = 0;
          double psi_real_temp_n1      = 0;
          double psi_imaginary_temp_n1 = 0;
          double psi_real_temp_n2      = 0;
          double psi_imaginary_temp_n2 = 0;
          double psi_real_temp_n3      = 0;
          double psi_imaginary_temp_n3 = 0;

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  -ds_operator_imaginary[0]*psi_imaginary[n0] /*- ds_operator_imaginary[1]*psi_imaginary[n3]*/;
          psi_imaginary_temp_n0 = ds_operator_imaginary[0]*psi_real[n0] /*+ ds_operator_imaginary[1]*psi_real[n3]*/
                                  +ds_operator_real[0]*psi_imaginary[n0]/* + ds_operator_real[1]*psi_imaginary[n3]*/;

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] /*+ ds_operator_real[3]*psi_real[n2]*/
                                  -ds_operator_imaginary[2]*psi_imaginary[n1] /*- ds_operator_imaginary[3]*psi_imaginary[n2]*/;
          psi_imaginary_temp_n1 = ds_operator_imaginary[2]*psi_real[n1] /*+ ds_operator_imaginary[3]*psi_real[n2]*/
                                  +ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;

          psi_real_temp_n2      = /*ds_operator_real[4]*psi_real[n1] +*/ ds_operator_real[5]*psi_real[n2]
                                  /*-ds_operator_imaginary[4]*psi_imaginary[n1]*/ - ds_operator_imaginary[5]*psi_imaginary[n2];
          psi_imaginary_temp_n2 = /*ds_operator_imaginary[4]*psi_real[n1] */+ ds_operator_imaginary[5]*psi_real[n2]
                                  /*+ds_operator_real[4]*psi_imaginary[n1]*/ + ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  /*-ds_operator_imaginary[6]*psi_imaginary[n0]*/ - ds_operator_imaginary[7]*psi_imaginary[n3];
          psi_imaginary_temp_n3 = /*ds_operator_imaginary[6]*psi_real[n0] */+ ds_operator_imaginary[7]*psi_real[n3]
                                  /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];

          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      }
    }
  }
}

/* Calculating energy expectation value w/o construct matrix
  Input:
    t: the current time point
  Output:
    average_energy: the energy for the system
  Side effect:
    psi_tmp_real[],psi_tmp_imaginary[]
  #: Only Calculate the energy for the system's spins
  #: not fully implement wit hx,hy
*/
double spin_system::energy(double t){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  double average_energy=0.;
  double check_img=0.;

  double hx=-1*h_x_start;
  for (int k = 0; k < N_sys; k++) {
    int i1=(int) pow(2,k);
    double hz=-1*h_z[k];
    #pragma omp parallel for default(none) shared(i1,hx,hz)
    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l -i2+i2/i1;
      int j = i+i1;
      /*sigma_x*/
      psi_tmp_real[i]     += Gamma*hx*psi_real[j];
      psi_tmp_imaginary[i]+= Gamma*hx*psi_imaginary[j];
      psi_tmp_real[j]     += Gamma*hx*psi_real[i];
      psi_tmp_imaginary[j]+= Gamma*hx*psi_imaginary[i];

      /*sigma_z*/
      psi_tmp_real[i]     += Delta*hz*psi_real[i];
      psi_tmp_imaginary[i]+= Delta*hz*psi_imaginary[i];
      psi_tmp_real[j]     += -Delta*hz*psi_real[j];
      psi_tmp_imaginary[j]+= -Delta*hz*psi_imaginary[j];

    }
  }

  for (int k = 0; k <N_sys ; k++) {
    for (int l = k+1; l < N_sys; l++) {
      double Jx=-1*J_x[k+l*N_sys];
      double Jy=-1*J_y[k+l*N_sys];
      double Jz=-1*J_z[k+l*N_sys];
      if(abs(Jx)>1e-15||abs(Jy)>1e-15||abs(Jz)>1e-15){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz)
        for (int m = 0; m < nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;
          /*sigma_x*sigma_x*/
          psi_tmp_real[n0]      += Delta*Jx*psi_real[n3];
          psi_tmp_imaginary[n0] += Delta*Jx*psi_imaginary[n3];
          psi_tmp_real[n1]      += Delta*Jx*psi_real[n2];
          psi_tmp_imaginary[n1] += Delta*Jx*psi_imaginary[n2];
          psi_tmp_real[n2]      += Delta*Jx*psi_real[n1];
          psi_tmp_imaginary[n2] += Delta*Jx*psi_imaginary[n1];
          psi_tmp_real[n3]      += Delta*Jx*psi_real[n0];
          psi_tmp_imaginary[n3] += Delta*Jx*psi_imaginary[n0];
          /*sigma_y*sigma_y*/
          psi_tmp_real[n0]      += -Delta*Jy*psi_real[n3];
          psi_tmp_imaginary[n0] += -Delta*Jy*psi_imaginary[n3];
          psi_tmp_real[n1]      += Delta*Jy*psi_real[n2];
          psi_tmp_imaginary[n1] += Delta*Jy*psi_imaginary[n2];
          psi_tmp_real[n2]      += Delta*Jy*psi_real[n1];
          psi_tmp_imaginary[n2] += Delta*Jy*psi_imaginary[n1];
          psi_tmp_real[n3]      += -Delta*Jy*psi_real[n0];
          psi_tmp_imaginary[n3] += -Delta*Jy*psi_imaginary[n0];
          /*sigma_z*sigma_z*/
          psi_tmp_real[n0]      += Delta*Jz*psi_real[n0];
          psi_tmp_imaginary[n0] += Delta*Jz*psi_imaginary[n0];
          psi_tmp_real[n1]      += -Delta*Jz*psi_real[n1];
          psi_tmp_imaginary[n1] += -Delta*Jz*psi_imaginary[n1];
          psi_tmp_real[n2]      += -Delta*Jz*psi_real[n2];
          psi_tmp_imaginary[n2] += -Delta*Jz*psi_imaginary[n2];
          psi_tmp_real[n3]      += Delta*Jz*psi_real[n3];
          psi_tmp_imaginary[n3] += Delta*Jz*psi_imaginary[n3];


        }
      }
    }
  }
  for (int i = 0; i < nofstates; ++i) {
    average_energy += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in functoin energy()   "<<check_img<<endl;

  return average_energy;
}

/////////////////////////working on new functoin
/* Calculating energy expectation value w/o construct matrix
  Input:
    t: the current time point
  Output:
    average_energy: the energy for the system
  Side effect:
    psi_tmp_real[],psi_tmp_imaginary[]
  #: Only Calculate the energy for the Environment's spins
*/
double spin_system::energy_env(double t){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  double average_energy=0.;
  double check_img=0.;


  for (int k = N-N_env; k <N ; k++) {
    int i1=(int) pow(2,k);
    double hx=-1*h_x[k];
    double hy=-1*h_y[k];
    double hz=-1*h_z[k];
    // if (abs(hx)>1e-8||abs(hz)>1e-8) {
    //   cout<<"H_env has non zero hx or hz inside energy_env() function";
    // }
    #pragma omp parallel for default(none) shared(i1,hx,hz,hy)
    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l -i2+i2/i1;
      int j = i+i1;
      /*sigma_x*/
      psi_tmp_real[i]     += hx*psi_real[j];
      psi_tmp_imaginary[i]+= hx*psi_imaginary[j];
      psi_tmp_real[j]     += hx*psi_real[i];
      psi_tmp_imaginary[j]+= hx*psi_imaginary[i];
      /*sigma_y*/
      psi_tmp_real[i]     += hy*psi_imaginary[j];
      psi_tmp_imaginary[i]+= -hy*psi_real[j];
      psi_tmp_real[j]     += -hy*psi_imaginary[i];
      psi_tmp_imaginary[j]+= hy*psi_real[i];
      /*sigma_z*/
      psi_tmp_real[i]     += hz*psi_real[i];
      psi_tmp_imaginary[i]+= hz*psi_imaginary[i];
      psi_tmp_real[j]     += -hz*psi_real[j];
      psi_tmp_imaginary[j]+= -hz*psi_imaginary[j];

    }
  }

  for (int k = N-N_env; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double Jx=-1*Jx_env[(k-N_sys)+(l-N_sys)*N_env];
      double Jy=-1*Jy_env[(k-N_sys)+(l-N_sys)*N_env];
      double Jz=-1*Jz_env[(k-N_sys)+(l-N_sys)*N_env];
      if(abs(Jx)>1e-15||abs(Jy)>1e-15||abs(Jz)>1e-15){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz)
        for (int m = 0; m < nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;
          /*sigma_x*sigma_x*/
          psi_tmp_real[n0]      += Jx*psi_real[n3];
          psi_tmp_imaginary[n0] += Jx*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jx*psi_real[n2];
          psi_tmp_imaginary[n1] += Jx*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jx*psi_real[n1];
          psi_tmp_imaginary[n2] += Jx*psi_imaginary[n1];
          psi_tmp_real[n3]      += Jx*psi_real[n0];
          psi_tmp_imaginary[n3] += Jx*psi_imaginary[n0];
          /*sigma_y*sigma_y*/
          psi_tmp_real[n0]      += -Jy*psi_real[n3];
          psi_tmp_imaginary[n0] += -Jy*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jy*psi_real[n2];
          psi_tmp_imaginary[n1] += Jy*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jy*psi_real[n1];
          psi_tmp_imaginary[n2] += Jy*psi_imaginary[n1];
          psi_tmp_real[n3]      += -Jy*psi_real[n0];
          psi_tmp_imaginary[n3] += -Jy*psi_imaginary[n0];
          /*sigma_z*sigma_z*/
          psi_tmp_real[n0]      += Jz*psi_real[n0];
          psi_tmp_imaginary[n0] += Jz*psi_imaginary[n0];
          psi_tmp_real[n1]      += -Jz*psi_real[n1];
          psi_tmp_imaginary[n1] += -Jz*psi_imaginary[n1];
          psi_tmp_real[n2]      += -Jz*psi_real[n2];
          psi_tmp_imaginary[n2] += -Jz*psi_imaginary[n2];
          psi_tmp_real[n3]      += Jz*psi_real[n3];
          psi_tmp_imaginary[n3] += Jz*psi_imaginary[n3];


        }
      }
    }
  }
  for (int i = 0; i < nofstates; ++i) {
    average_energy += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in functoin energy()   "<<check_img<<endl;

  return average_energy;
}
////////////////////////////////////////////////

/////////////////////////working on new functoin
/* Calculating energy expectation value w/o construct matrix
  Input:
    t: the current time point
  Output:
    average_energy: the energy for the system
  Side effect:
    psi_tmp_real[],psi_tmp_imaginary[]
  #: Only Calculate the energy for H_se
*/

double spin_system::energy_se(double t){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  double average_energy=0.;
  double check_img=0.;

  // for H_se, there are not single interaction. so no h_z,h_y,h_z
  for (int k = 0; k <N_sys ; k++) {
    for (int l = N-N_env; l < N; l++) {
      double Jx=-1*Jx_se[k+(l-N_sys)*N_sys];
      double Jy=-1*Jy_se[k+(l-N_sys)*N_sys];
      double Jz=-1*Jz_se[k+(l-N_sys)*N_sys];
      if(abs(Jx)>1e-15||abs(Jy)>1e-15||abs(Jz)>1e-15){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz)
        for (int m = 0; m < nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;
          /*sigma_x*sigma_x*/
          psi_tmp_real[n0]      += Jx*psi_real[n3];
          psi_tmp_imaginary[n0] += Jx*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jx*psi_real[n2];
          psi_tmp_imaginary[n1] += Jx*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jx*psi_real[n1];
          psi_tmp_imaginary[n2] += Jx*psi_imaginary[n1];
          psi_tmp_real[n3]      += Jx*psi_real[n0];
          psi_tmp_imaginary[n3] += Jx*psi_imaginary[n0];
          /*sigma_y*sigma_y*/
          psi_tmp_real[n0]      += -Jy*psi_real[n3];
          psi_tmp_imaginary[n0] += -Jy*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jy*psi_real[n2];
          psi_tmp_imaginary[n1] += Jy*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jy*psi_real[n1];
          psi_tmp_imaginary[n2] += Jy*psi_imaginary[n1];
          psi_tmp_real[n3]      += -Jy*psi_real[n0];
          psi_tmp_imaginary[n3] += -Jy*psi_imaginary[n0];
          /*sigma_z*sigma_z*/
          psi_tmp_real[n0]      += Jz*psi_real[n0];
          psi_tmp_imaginary[n0] += Jz*psi_imaginary[n0];
          psi_tmp_real[n1]      += -Jz*psi_real[n1];
          psi_tmp_imaginary[n1] += -Jz*psi_imaginary[n1];
          psi_tmp_real[n2]      += -Jz*psi_real[n2];
          psi_tmp_imaginary[n2] += -Jz*psi_imaginary[n2];
          psi_tmp_real[n3]      += Jz*psi_real[n3];
          psi_tmp_imaginary[n3] += Jz*psi_imaginary[n3];


        }
      }
    }
  }
  for (int i = 0; i < nofstates; ++i) {
    average_energy += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in functoin energy()   "<<check_img<<endl;

  return average_energy;
}
////////////////////////////////////////////////

/////////////////////////working on new functoin

/* Calculating energy expectation value w/o construct matrix for all spins
  Input:
    t: the current time point
  Output:
    average_energy: the energy for the system
  Side effect:
    psi_tmp_real[],psi_tmp_imaginary[]
  #: Only Calculate the energy for the environment's spins
*/
double spin_system::energy_all(double t){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  double average_energy=0.;
  double check_img=0.;

  for (int k = 0; k < N; k++) {
    double hx=0.;
    double hy=0.;
    double hz=0.;

    if (k>=N_sys) {
      hx=-1*h_x[k];
      hy=-1*h_y[k];
      hz=-1*h_z[k];
    } else {
      hx=-1*h_x_start*Gamma;
      hy=0.;
      hz=-1*h_z[k]*Delta;
    }

    if(abs(hx)>1e-15||abs(hz)>1e-15||abs(hy)>1e-15){
      int i1=(int) pow(2,k);
      #pragma omp parallel for default(none) shared(i1,hx,hz,hy)
      for (int l = 0; l < nofstates; l+=2) {
        int i2= l & i1;
        int i = l -i2+i2/i1;
        int j = i+i1;
        /*sigma_x*/
        psi_tmp_real[i]     += hx*psi_real[j];
        psi_tmp_imaginary[i]+= hx*psi_imaginary[j];
        psi_tmp_real[j]     += hx*psi_real[i];
        psi_tmp_imaginary[j]+= hx*psi_imaginary[i];
        /*sigma_y*/
        psi_tmp_real[i]     += hy*psi_imaginary[j];
        psi_tmp_imaginary[i]+= -hy*psi_real[j];
        psi_tmp_real[j]     += -hy*psi_imaginary[i];
        psi_tmp_imaginary[j]+= hy*psi_real[i];
        /*sigma_z*/
        psi_tmp_real[i]     += hz*psi_real[i];
        psi_tmp_imaginary[i]+= hz*psi_imaginary[i];
        psi_tmp_real[j]     += -hz*psi_real[j];
        psi_tmp_imaginary[j]+= -hz*psi_imaginary[j];

      }
    }
  }

  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double Jx=0.;
      double Jy=0.;
      double Jz=0.;
      if(k>=N_sys){
        Jx=-1*Jx_env[(k-N_sys)+(l-N_sys)*N_env];
        Jy=-1*Jy_env[(k-N_sys)+(l-N_sys)*N_env];
        Jz=-1*Jz_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys){
        Jx=-1*Jx_se[k+(l-N_sys)*N_sys];
        Jy=-1*Jy_se[k+(l-N_sys)*N_sys];
        Jz=-1*Jz_se[k+(l-N_sys)*N_sys];
      } else {
        Jx=-1*J_x[k+l*N_sys]*Delta;
        Jy=-1*J_y[k+l*N_sys]*Delta;
        Jz=-1*J_z[k+l*N_sys]*Delta;
      }

      if(abs(Jx)>1e-15||abs(Jy)>1e-15||abs(Jz)>1e-15){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz)
        for (int m = 0; m < nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;
          /*sigma_x*sigma_x*/
          psi_tmp_real[n0]      += Jx*psi_real[n3];
          psi_tmp_imaginary[n0] += Jx*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jx*psi_real[n2];
          psi_tmp_imaginary[n1] += Jx*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jx*psi_real[n1];
          psi_tmp_imaginary[n2] += Jx*psi_imaginary[n1];
          psi_tmp_real[n3]      += Jx*psi_real[n0];
          psi_tmp_imaginary[n3] += Jx*psi_imaginary[n0];
          /*sigma_y*sigma_y*/
          psi_tmp_real[n0]      += -Jy*psi_real[n3];
          psi_tmp_imaginary[n0] += -Jy*psi_imaginary[n3];
          psi_tmp_real[n1]      += Jy*psi_real[n2];
          psi_tmp_imaginary[n1] += Jy*psi_imaginary[n2];
          psi_tmp_real[n2]      += Jy*psi_real[n1];
          psi_tmp_imaginary[n2] += Jy*psi_imaginary[n1];
          psi_tmp_real[n3]      += -Jy*psi_real[n0];
          psi_tmp_imaginary[n3] += -Jy*psi_imaginary[n0];
          /*sigma_z*sigma_z*/
          psi_tmp_real[n0]      += Jz*psi_real[n0];
          psi_tmp_imaginary[n0] += Jz*psi_imaginary[n0];
          psi_tmp_real[n1]      += -Jz*psi_real[n1];
          psi_tmp_imaginary[n1] += -Jz*psi_imaginary[n1];
          psi_tmp_real[n2]      += -Jz*psi_real[n2];
          psi_tmp_imaginary[n2] += -Jz*psi_imaginary[n2];
          psi_tmp_real[n3]      += Jz*psi_real[n3];
          psi_tmp_imaginary[n3] += Jz*psi_imaginary[n3];


        }
      }
    }
  }
  for (int i = 0; i < nofstates; ++i) {
    average_energy += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in functoin energy()   "<<check_img<<endl;

  return average_energy;
}
////////////////////////////////////////////////


/* Calculating the spin expectation value
  Input:
    d: 'x' for sigma_x, 'y' for sigma_y, 'z' for sigma_z
    which_spin: which spin we want to calculate (count from 0)
  Output:
    average_spin: the spin expectation value.
  Side Effect:
    psi_tmp_real[],psi_tmp_imaginary[]
*/
double spin_system::spin(char d,int which_spin){

  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  double average_spin=0.;
  double check_img=0.;

  for (int k = 0; k < N; k++) {
    if (k==which_spin) {
      int i1=(int) pow(2,k);
      #pragma omp parallel for default(none) shared(i1,d,std::cout)
      for (int l = 0; l < nofstates; l+=2) {
        int i2= l & i1;
        int i = l -i2+i2/i1;
        int j = i+i1;

        if(d=='x'){
          /*sigma_x*/
          psi_tmp_real[i]     += 0.5*psi_real[j];
          psi_tmp_imaginary[i]+= 0.5*psi_imaginary[j];
          psi_tmp_real[j]     += 0.5*psi_real[i];
          psi_tmp_imaginary[j]+= 0.5*psi_imaginary[i];
        }
        else if(d=='y'){
          /*sigma_y*/
          psi_tmp_real[i]     += 0.5*psi_imaginary[j];
          psi_tmp_imaginary[i]+= -0.5*psi_real[j];
          psi_tmp_real[j]     += -0.5*psi_imaginary[i];
          psi_tmp_imaginary[j]+= 0.5*psi_real[i];
        }
        else if(d=='z'){
         /*sigma_z*/
         psi_tmp_real[i]     += 0.5*psi_real[i];
         psi_tmp_imaginary[i]+= 0.5*psi_imaginary[i];
         psi_tmp_real[j]     += -0.5*psi_real[j];
         psi_tmp_imaginary[j]+= -0.5*psi_imaginary[j];
       }

        else {
          cout<<"WRONG input for calculating average spin"<<endl;
        }
      }
    }
  }

  for (int i = 0; i < nofstates; ++i) {
    // cout<<2++1*3<<endl;
    average_spin += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    // if (which_spin>=N_sys && (d=='x'||d=='y'||d=='z')) {
    //  cout<<which_spin<<": "<<d<<": "<<i<<": "<<average_spin<<endl;
    // }

    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }

  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in function spin()"<<check_img<<endl;
  return average_spin;


}

// /* calculate the average spin in one shot. can reduce the computing time.
//   Input:
//   Output:
//     change the array spin_return[]
// */
//   void spin_system::spin_allinone(){
//   for (int i = 0; i < N*3; i++){
//     spin_return[i]=0;
//   }
//
//   double check_img=0.;
//   for (int k = 0; k < N; k++) {
//     for (int i = 0; i < nofstates; i++) {
//       psi_tmp_real[i] = 0;
//       psi_tmp_imaginary[i] = 0;
//       psi_tmp_y_real[i] = 0;
//       psi_tmp_y_imaginary[i] = 0;
//       psi_tmp_z_real[i] = 0;
//       psi_tmp_z_imaginary[i] = 0;
//     }
//     int i1=(int) pow(2,k);
//     for (int d = 0; d < 3; d++) { //d = 0,1,2 means x,y,z
//         /* code */
//       #pragma omp parallel for default(none) shared(i1,std::cout,d)
//       for (int l = 0; l < nofstates; l+=2) {
//         int i2= l & i1;
//         int i = l -i2+i2/i1;
//         int j = i+i1;
//
//         if(0==d){
//           /*sigma_x*/
//           psi_tmp_real[i]     += 0.5*psi_real[j];
//           psi_tmp_imaginary[i]+= 0.5*psi_imaginary[j];
//           psi_tmp_real[j]     += 0.5*psi_real[i];
//           psi_tmp_imaginary[j]+= 0.5*psi_imaginary[i];
//         } else if (1==d){
//           /*sigma_y*/
//           psi_tmp_y_real[i]     += 0.5*psi_imaginary[j];
//           psi_tmp_y_imaginary[i]+= -0.5*psi_real[j];
//           psi_tmp_y_real[j]     += -0.5*psi_imaginary[i];
//           psi_tmp_y_imaginary[j]+= 0.5*psi_real[i];
//         } else {
//           /*sigma_z*/
//           psi_tmp_z_real[i]     += 0.5*psi_real[i];
//           psi_tmp_z_imaginary[i]+= 0.5*psi_imaginary[i];
//           psi_tmp_z_real[j]     += -0.5*psi_real[j];
//           psi_tmp_z_imaginary[j]+= -0.5*psi_imaginary[j];
//         }
//       }
//
//       // for (int i = 0; i < nofstates; ++i) {
//       //
//       //   spin_return[(k*3)]   += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
//       //   spin_return[(k*3)+1] += psi_real[i]*psi_tmp_y_real[i] - -1*psi_imaginary[i]*psi_tmp_y_imaginary[i];
//       //   spin_return[(k*3)+2] += psi_real[i]*psi_tmp_z_real[i] - -1*psi_imaginary[i]*psi_tmp_z_imaginary[i];
//       //   check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
//       // }
//       if(0==d){
//         for (int i = 0; i < nofstates; ++i)
//           spin_return[(k*3)]   += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
//       } else if(1==d) {
//         for (int i = 0; i < nofstates; ++i)
//           spin_return[(k*3)+1] += psi_real[i]*psi_tmp_y_real[i] - -1*psi_imaginary[i]*psi_tmp_y_imaginary[i];
//       } else {
//         for (int i = 0; i < nofstates; ++i)
//           spin_return[(k*3)+2] += psi_real[i]*psi_tmp_z_real[i] - -1*psi_imaginary[i]*psi_tmp_z_imaginary[i];
//       }
//     }
//     if (abs(check_img)>1e-13)
//       cout<<"Something went wrong in function spin()"<<check_img<<endl;
//   }
// }


/* Read in the mmber from a txt into an array
  Input:
    N: how many elements we want to read.
    Array: the array to store this elements
    filename: the name of the file
*/
void spin_system::read(int N, double* Array, char const * filename ){
  /* Set the input class myfile and open the file "filename".
    Check whether it is open with .is_open boolean.
  */
  // cout<<"inside read() N= "<<filename<<" "<<N<<endl;
  ifstream myfile;
  myfile.open(filename);
  if (!myfile.is_open()) {
    cout<<"Unable to open the file "<<filename<<endl;
  }

  /* Input the value in the file to tmp,
    then put them to the Array.
    Close the file at the end.
  */
  double tmp;
  int i=0;
  while (myfile>>tmp && i<N) {
    Array[i]=tmp;
    i++;
  }
  myfile.close();

  /* Check the read values.
  */
  for (i = 0; i < N; i++) {
    // cout<<filename<<" "<<i<<"th element ="<<Array[i]<< endl;
  }
}


/* Set up an environment by diagonalize the H matrix.
  Input:
    N_env: the # of spins for Environment.
    Temperature: temperature
  Side Effect:
    w[],z[]
*/
void spin_system::environment(int N_env, double Temperature){


  //Construct the H_env matrix for the later use for solving energy of the heat bath
  double Boltzmann_const= 1;//we use 1 here insted of 0.000086173303 ev/K;
  int nofstates=(int) pow(2,N_env);
  complex<double>* H_env;
  H_env=new complex<double> [nofstates*(nofstates+1)/2]();
  // for (int i = 0; i < nofstates*(nofstates+1)/2; i++) {
    // H_env[i]+=1.;
    // cout<<H_env[i]<<endl;
  // }
  /* 22.02.2017 Add single spin*/
  for (int k = 0; k < N_env; k++) {
    int i1=(int) pow(2,k);
    double hx=-1*h_x[k+N_sys];
    double hy=-1*h_y[k+N_sys];
    double hz=-1*h_z[k+N_sys];
    cout<<"hx hy hz ="<<hx<<" "<<hy<<" "<<hz<<" "<<endl;

    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l - i2 +i2/i1;
      int j = i+i1;
        H_env[i+i*(i+1)/2] += hz;
        H_env[j+j*(j+1)/2] += -hz;
        H_env[i+j*(j+1)/2].imag() += -hy;
        H_env[i+j*(j+1)/2] += hx;
    }
  }

  for (int k = 0; k < N_env; k++) {
    for (int l = k+1; l < N_env; l++) {
      double Jx=-1*Jx_env[k+l*N_env];
      double Jy=-1*Jy_env[k+l*N_env];
      double Jz=-1*Jz_env[k+l*N_env];
      if(Jx!=0||Jy!=0||Jz!=0){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        for (int m=0; m<nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2 - n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;

          H_env[n0+n0*(n0+1)/2] += 1.*Jz;
          H_env[n1+n1*(n1+1)/2] +=-1.*Jz;
          H_env[n2+n2*(n2+1)/2] +=-1.*Jz;
          H_env[n3+n3*(n3+1)/2] += 1.*Jz;

          H_env[n0+n3*(n3+1)/2] += 1.*Jx;
          H_env[n1+n2*(n2+1)/2] += 1.*Jx;
          // H[n2+n1*(n1+1)/2] += 1.*Jx;
          // H[n3+n0*(n0+1)/2] += 1.*Jx;

          H_env[n0+n3*(n3+1)/2] +=-1.*Jy;
          H_env[n1+n2*(n2+1)/2] += 1.*Jy;
          // H[n2+n1*(n1+1)/2] += 1.*Jy;
          // H[n3+n0*(n0+1)/2] +=-1.*Jy;

        }
      }
    }
  }

  //start preparing variable for lapack use
  //solving the energy to w[] array. then use it as the coefficient of |E_B>.
  int n=nofstates;
  double vl,vu;
  int il,iu;
  double abstol = 2*DLAMCH("S");
  int m;
  complex<double> ap[nofstates*(nofstates+1)/2];
  for (int i = 0; i < nofstates*(nofstates+1)/2; i++){
    ap[i]=H_env[i];
  }
  int lda = nofstates;
  complex<double> work[2*nofstates];
  double rwork[7*nofstates];
  int iwork[5*nofstates];
  int ifail[nofstates];
  int info;
  zhpevx("Vector","All","Upper", &n, ap, &vl, &vu, &il, &iu, &abstol, &m, w, z, &lda, work, rwork, iwork, ifail, &info );
  if(info!=0){
    cout<<"info = "<<info<<endl;
  }
  ofstream eng_out("eng_spec.dat");
  for (int i = 0; i < nofstates; i++) {
    eng_out<<w[i]<<endl;
  }

  double sum=0.;
  for (int i = 1; i < nofstates; i++) {
    w[i]=exp(-1*(w[i]-w[0])/(Temperature*Boltzmann_const));
    sum+=w[i];
  }
  w[0]=1;//exp(-1*(w[0]-w[0])/(Temperature*Boltzmann_const));
  sum+=w[0];
  cout<<"sum of z: "<<sum<<endl;
  for (int i = 0; i < nofstates; i++) {
    w[i]=w[i]/sum;
  }

  ofstream envr_out("env_partition_factor.dat");
  ofstream env_basis_out("env_basis.dat");
  for (int i = 0; i < nofstates; i++) {
    envr_out<<w[i]<<endl;
  }

  for (int i = 0; i < nofstates; i++) {
    for (int j = 0; j < nofstates; j++) {
      env_basis_out<<z[i*nofstates+j]<<" ";
    }
    env_basis_out<<endl;
  }

}

/* It make a driect product of basis states of E_i and Sys.This is a replace function for generate().
  Input:
    n: the #th of the E_i.(Which Ei's eigenvector you want to calculate with)
    array_real: the output real part
    array_imagine: the output imag part
    complex<double>* env: eigenvector of Env get from Lapack
    sys_real: the real part of psi of system
    sys_imag: the imag part of psi of system
  Side Effect:
    array_real[],array_imagine[]
  #: usually the side effect is on psi_real[] and psi_imaginary[]
*/
void spin_system::direct_product(int n, double* array_real, double* array_imagine, complex<double>* env, double* sys_real, double* sys_imag){

  int nos_sys=(int) pow(2,N_sys);
  int nos_env=(int) pow(2,N_env);
  for (int i = 0; i < nos_env; i++) {
    for (int j = 0; j < nos_sys; j++) {
      array_real[i*nos_sys+j]=env[n*nos_env+i].real()*sys_real[j]-env[n*nos_env+i].imag()*sys_imag[j];
      array_imagine[i*nos_sys+j]=env[n*nos_env+i].imag()*sys_real[j]+env[n*nos_env+i].real()*sys_imag[j];
      // cout<<array_real[i*nos_sys+j]<<" "<<array_imagine[i*nos_sys+j]<<endl;
    }
  }
}



/* !!!!!!!!WILL BE REMOVE in the future, since it might be a wrong implementation!!!!!!!
  Try to read the initial basis state from the system and the Enivironment
  And then combine then together into a new state.
  Input:
    int: the # of total spin.
    double*: the real part
    double*: the imaginary part
    char const*: the basis state of the system
    char const*: the basis state of the environment
*/
void spin_system::generate(int N, double* array_real, double* array_imagine, char const* filename_sysr, char const* filename_sysi, char const* filename_envr, char const* filename_envi ){

  double *sys_real, *sys_imag, *env_real, *env_imag;
  // int N_half=(int) pow(2,N/2);
  int nos_sys=(int) pow(2,N_sys);
  int nos_env=(int) pow(2,N_env);
  sys_real=new double[nos_sys];
  sys_imag=new double[nos_sys];
  env_real=new double[nos_env];
  env_imag=new double[nos_env];
  read(nos_sys, sys_real, filename_sysr);
  read(nos_sys, sys_imag, filename_sysi);
  read(nos_env, env_real, filename_envr);
  read(nos_env, env_imag, filename_envi);

  ofstream state_out("state_complete.dat");
  for (int i = 0; i < nos_env; i++) {
    for (int j = 0; j < nos_sys; j++) {
      array_real[i*nos_sys+j]=env_real[i]*sys_real[j]-env_imag[i]*sys_imag[j];
      array_imagine[i*nos_sys+j]=env_imag[i]*sys_real[j]+env_real[i]*sys_imag[j];
      state_out<<array_real[i*nos_sys+j]<<" "<<array_imagine[i*nos_sys+j]<<endl;
    }
  }
  double sum=0;
  int nos=(int) pow(2,N_sys+N_env);
  for (int i = 0; i < nos; i++) {
    sum+=array_real[i]*array_real[i]+array_imagine[i]*array_imagine[i];
  }
  // cout<<"sum inside generate functoin: "<<sum<<endl;
}

/* Make J_env txt file randomly between [-1,1]
  It is the coupling factor for environment spins
  Input:
    N: # of spins for Environment;
    G: strength
  Side Effect:
    Jx_env[],Jy_env[],Jz_env[]
*/
void spin_system::Jenv_generate(int N, double G){
  ofstream Jx_env_out("Jx_env.txt");
  ofstream Jy_env_out("Jy_env.txt");
  ofstream Jz_env_out("Jz_env.txt");
  srand(time(NULL));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (j>=i) {
        Jx_env_out<<0<<" ";
        Jy_env_out<<0<<" ";
        Jz_env_out<<0<<" ";
      } else {
        Jx_env_out<<G*(((double) rand()/RAND_MAX)*2-1)<<" ";
        Jy_env_out<<G*(((double) rand()/RAND_MAX)*2-1)<<" ";
        Jz_env_out<<G*(((double) rand()/RAND_MAX)*2-1)<<" ";
      }
    }
    Jx_env_out<<endl;
    Jy_env_out<<endl;
    Jz_env_out<<endl;
  }
  read(N*N,Jx_env,"Jx_env.txt");
  read(N*N,Jy_env,"Jy_env.txt");
  read(N*N,Jz_env,"Jz_env.txt");

}



/* Make J_se txt file randomly between [-1,1]
  and then multiply with a global interaction factor G,
  control how strong the system interact with the environment.
  Input:
    N_sys: # of spins for system
    N_env: # of spins for environment
    G: stenght factor
  Side Effect:
    Jx_se[],Jy_se[],Jz_se[]
*/
  void spin_system::Jse_generate(int N_sys, int N_env, double G){
    ofstream Jx_se_out("Jx_se.txt");
    ofstream Jy_se_out("Jy_se.txt");
    ofstream Jz_se_out("Jz_se.txt");
    srand(time(NULL));
    for (int i = 0; i < N_sys; i++) {
      for (int j = 0; j < N_env; j++) {
        Jx_se_out<<(((double) rand()/RAND_MAX)*2-1)*G <<" ";
        Jy_se_out<<(((double) rand()/RAND_MAX)*2-1)*G <<" ";
        Jz_se_out<<(((double) rand()/RAND_MAX)*2-1)*G <<" ";
      }
      Jx_se_out<<endl;
      Jy_se_out<<endl;
      Jz_se_out<<endl;
    }
    read(N_sys*N_env,Jx_se,"Jx_se.txt");
    read(N_sys*N_env,Jy_se,"Jy_se.txt");
    read(N_sys*N_env,Jz_se,"Jz_se.txt");
  }


/* The main process to run the simulation
  16.12.2016: I haven't add the time evolution part. It should be added after.
  26.12.2016: time evolution part added.
  03.02.2017: can handle with sigma_(x,y,z) and sigma_(x,y,z)*sigma_(x,y,z) respectly.
  20.02.2017: Run with environment heat bath.
*/
void spin_system::run(){

  int total_steps=0;
  total_steps=(int) T/tau;
  double* frequency;
  frequency=new double [total_steps+1]();
  /////test not go over whole J again and again/////
  //////////////////////////////////////////////////
  // count_z=0;
  // count_y=0;
  // count_x=0;
  // Jz_k_marked=new double [N_sys*(N_sys+1)/2];
  // Jz_l_marked=new double [N_sys*(N_sys+1)/2];
  // Jy_k_marked=new double [N_sys*(N_sys+1)/2];
  // Jy_l_marked=new double [N_sys*(N_sys+1)/2];
  // Jx_k_marked=new double [N_sys*(N_sys+1)/2];
  // Jx_l_marked=new double [N_sys*(N_sys+1)/2];
  // for (int k = 0; k <N_sys ; k++) {
  //   for (int l = k+1; l < N_sys; l++) {
  //     if (abs(J_z[k+l*N_sys])>1e-15){
  //       Jz_k_marked[count_z]=k;
  //       Jz_l_marked[count_z]=l;
  //       count_z+=1;
  //     }
  //     if (abs(J_y[k+l*N_sys])>1e-15){
  //       Jy_k_marked[count_y]=k;
  //       Jy_l_marked[count_y]=l;
  //       count_y+=1;
  //     }
  //     if (abs(J_x[k+l*N_sys])>1e-15){
  //       Jx_k_marked[count_x]=k;
  //       Jx_l_marked[count_x]=l;
  //       count_x+=1;
  //     }
  //   }
  // }
  //
  // cout<<"Elements in Jx, Jy, and Jz: "<<count_x<<" "<<count_y<<" "<<count_z<<" "<<endl;
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  double norm=0.;
  for (int E_i = 0; E_i < (int) pow(2,N_env); E_i++) {
    direct_product(E_i,psi_real,psi_imaginary,z,psi_sys_real,psi_sys_imaginary);
    for (int i = 0; i < nofstates; i++) {
      coefficient_return[i]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);

    }
  }
  for (int i = 0; i < nofstates; i++) {
    norm+=coefficient_return[i];
  }
  cout<<"norm before run: "<<norm<<endl;
  for (int i = 0; i < nofstates; i++) {
    coefficient_return[i]=0.;
  }

  for (int E_i = 0; E_i < (int) pow(2,N_env); E_i++) {
    if (abs(w[E_i]-0)<1e-8)
      continue;
    direct_product(E_i,psi_real,psi_imaginary,z,psi_sys_real,psi_sys_imaginary);

    for (int step = 0; step < total_steps+1; step++){ //+1 because i count the 0 point and the last poing as well.

      Delta=step*tau/T;
      Gamma=1-Delta;

      if (step%500==0)
        cout<<"E_i= "<<E_i<<", w[]= "<<w[E_i]<<", step: "<<step<<endl;

      energy_sys_return[step]+=w[E_i]*energy(step*tau);
      energy_env_return[step]+=w[E_i]*energy_env(step*tau);
      energy_se_return[step]+=w[E_i]*energy_se(step*tau);
      energy_all_return[step]+=w[E_i]*energy_all(step*tau);
      for (int s = 0; s < N; s++) {
        int index=step*N*3+s*3;
        spin_return[index]  +=w[E_i]*spin('x',s);
        spin_return[index+1]+=w[E_i]*spin('y',s);
        spin_return[index+2]+=w[E_i]*spin('z',s);
      }

      for (int i = 119; i < nofstates; i+=256) {
        frequency[step]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);
      }


      single_spin_op(step*tau);
      double_spin_op_x(step*tau);
      double_spin_op_y(step*tau);
      double_spin_op_z(step*tau);
      double_spin_op_y(step*tau);
      double_spin_op_x(step*tau);
      single_spin_op(step*tau);

    }
    for (int i = 0; i < nofstates; i++) {
      coefficient_return[i]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);
    }
  }
  // output the return value: coefficient, energy expectation value, and spin expectation value.
  ofstream Coefficient_out("coefficient.dat");
  for (int i = 0; i < nofstates; i++) {
    Coefficient_out<<coefficient_return[i]<<endl;
  }
  ofstream output("output_ind.dat");
  output<<"Time Energy_sys Energy_env Energy_se Energy_all Frequency ";
  for (int i = 0; i < N; i++) {
    output<<"Sx_"<<i<<" "<<"Sy_"<<i<<" "<<"Sz_"<<i<<" ";
  }
  output<<endl;
  for (int step = 0; step < total_steps+1; step++){
    output<<step*tau<<" ";
    output<<energy_sys_return[step]<<" ";
    output<<energy_env_return[step]<<" ";
    output<<energy_se_return[step]<<" ";
    output<<energy_all_return[step]<<" ";
    output<<frequency[step]<<" ";
    for (int i = 0; i < 3*N; i++) {
      output<<spin_return[step*3*N+i]<<" ";
    }
    output<<endl;
  }

}
