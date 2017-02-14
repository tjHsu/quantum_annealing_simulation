//#include "spinsys.h"


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

// #define Boltzmann_const 0.000086173303 //m2 kg s-2 K-1

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
  double* Jx_env;
  double* Jy_env;
  double* Jz_env;
  double* Jx_se;
  double* Jy_se;
  double* Jz_se;
  double* Jx_com;
  double* Jy_com;
  double* Jz_com;

  double G;
  double* h_x;
  double* h_y;
  double* h_z;
  double  h_x_start;
  double Gamma;
  double Delta;

//optimize test
  double* Jz_k_marked;
  double* Jz_l_marked;
  int count_z;
  double* Jy_k_marked;
  double* Jy_l_marked;
  int count_y;
  double* Jx_k_marked;
  double* Jx_l_marked;
  int count_x;

  void single_spin_op(double);
  void double_spin_op(double);
  void double_spin_op_x(double);
  void double_spin_op_y(double);
  void double_spin_op_z(double);
  void generate_initial_sys_state(char const *);
  void energy(double);
  double spin(char,int);

  void read(int,double*,char const *);
  void generate(int, double*, double*, char const *, char const *, char const *, char const *);



public:
  void initialize(int, int ,double ,double);
  void run();
  void environment(int, double);
  double* psi_real;
  double* psi_imaginary;
  double* psi_tmp_real;
  double* psi_tmp_imaginary;
  double* H_real;
  double* H_imaginary;
  double average_energy;
  double average_spin;
  void Jenv_generate(int,double);
  void Jse_generate(int, int, double);


};

int main(int argc, char* argv[]){

  spin_system test;
  test.initialize(8,8,100,0.1);


  int N=16;
  int nofstates=(int) pow(2,N);

  // test.initialize(8,0,100.,0.01);




  double norm=0.;
  for (int i = 0; i < nofstates; i++) {
    norm+=test.psi_real[i]*test.psi_real[i]+test.psi_imaginary[i]*test.psi_imaginary[i];
  }
  cout<<"norm before run = "<<norm<<endl;
  time_t start=time(0);
  clock_t t;
  t = clock();
  test.run();
  t =clock()-t;
  time_t end=time(0);
  double time=difftime(end,start);

  cout<<"It took me "<<t<<" clicks ( "<<((float) t)/CLOCKS_PER_SEC<<" processing seconds.)"<<endl;
  cout<<"costed me "<<time<<" second in real time"<<endl;

  norm =0.;

  double max=0.;
  int location;
  for (int i = 0; i < nofstates; i++) {
    double tmp=test.psi_real[i]*test.psi_real[i]+test.psi_imaginary[i]*test.psi_imaginary[i];
    norm+=tmp;
    if (tmp > max) {
      max = tmp;
      location = i;
    }
  }
  cout<<"norm after run = "<<norm<<endl;
  cout<<"max after run = "<<max<<endl;
  cout<<"location after run = "<<location<<endl;

  ofstream Coefficient_out("coefficient.dat");
  double sum_GS=0.;
  double sum_ES=0.;
  for (int i = 0; i < nofstates; i++) {
    double tmp=test.psi_real[i]*test.psi_real[i]+test.psi_imaginary[i]*test.psi_imaginary[i];
    if(tmp<1e-8){
      Coefficient_out<<0<<endl;
      sum_ES+=tmp;
    } else {
      if (0==(i-176)%256) {
        Coefficient_out<<tmp<<endl;
        // cout<<(i-176)/256<<"th GS= "<<tmp<<endl;
        sum_GS+=tmp;
      } else {
        Coefficient_out<<tmp<<endl;
        if (tmp>1e-2)
          cout<<i<<"th= "<<tmp<<endl;
        sum_ES+=tmp;
      }

    }
  }

  cout<<endl;
  cout<<"It took me "<<t<<" clicks ( "<<((float) t)/CLOCKS_PER_SEC<<" processing seconds.)"<<endl;
  cout<<"costed me "<<time<<" second in real time"<<endl;
  cout<<"norm after run = "<<norm<<endl;
  cout<<"max after run = "<<max<<endl;
  cout<<"location after run = "<<location<<endl;
  cout<<endl;
  cout<<"pobability of GS: "<<sum_GS<<endl;
  cout<<"pobability of non-GS: "<<sum_ES<<endl;
  cout<<"sum of two"<<sum_ES+sum_GS<<endl;

  cout<<"reutrn "<<0<<endl;
  return 0;
}


/* initialize the system
  N: number of spins
  T: end time
  tau: the size of the timestep
*/
void spin_system::initialize(int N_sys_user_defined, int N_env_user_defined, double T_user_defined, double tau_user_defined){
  N_sys = N_sys_user_defined;
  N_env = N_env_user_defined;
  N = N_sys+N_env;
  T = T_user_defined;
  tau = tau_user_defined;
  nofstates = (int) pow(2,N);// number of states construct by the number of spins

  /* initialize the coupling factor J for double spin Hamiltonian,
    and the factor h for single spin Hamiltonian.
    J should have size as N*(N-1) and the size of h is simply N.
  */
  J_x = new double [N_sys*N_sys];
  J_y = new double [N_sys*N_sys];
  J_z = new double [N_sys*N_sys];
  Jx_env = new double [N_env*N_env];
  Jy_env = new double [N_env*N_env];
  Jz_env = new double [N_env*N_env];
  Jx_se = new double [N_env*N_sys];
  Jy_se = new double [N_env*N_sys];
  Jz_se = new double [N_env*N_sys];
  h_x = new double [N];
  h_y = new double [N];
  h_z = new double [N];
  h_x_start=1;
  for (int i = 0; i < N_sys*N_sys; i++){
    J_x[i] = 0.;
    J_y[i] = 0.;
    J_z[i] = 0.;
  }
  for (int i = 0; i < N_env*N_env; i++){
    Jx_env[i] = 0.;
    Jy_env[i] = 0.;
    Jz_env[i] = 0.;
  }
  for (int i = 0; i < N_env*N_sys; i++){
    Jx_se[i] = 0.;
    Jy_se[i] = 0.;
    Jz_se[i] = 0.;
  }
  for (int i = 0; i < N; i++){
    h_x[i]=0.;
    h_y[i]=0.;
    h_z[i]=0.;
  }


  read(N_sys*N_sys,J_z,"J2.txt");
  read(N_sys*N_sys,J_x,"J2x.txt");
  read(N_sys*N_sys,J_y,"J2y.txt");
  read(N_env*N_env,Jx_env,"Jx_env.txt");
  read(N_env*N_env,Jy_env,"Jy_env.txt");
  read(N_env*N_env,Jz_env,"Jz_env.txt");
  read(N_sys*N_env,Jx_se,"Jx_se.txt");
  read(N_sys*N_env,Jy_se,"Jy_se.txt");
  read(N_sys*N_env,Jz_se,"Jz_se.txt");
  read(N,h_z,"h2.txt");

  // G=1.0;
  // Jenv_generate(N_env,G);
  G=0.0;
  Jse_generate(N_sys,N_env,G);

  // environment(N_env,0.00001);

  /* initialize the wave function in the ground state
  */

  psi_real = new double [nofstates];
  psi_imaginary = new double [nofstates];
  for (int i = 0; i < nofstates; i++) {
      psi_real[i]      = 0;
      psi_imaginary[i] = 0;
  }
  generate_initial_sys_state("allx");
  generate(N,psi_real,psi_imaginary,"psi_real.dat","psi_imagine.dat","env_real.dat","env_imag.dat");
  psi_tmp_real = new double [nofstates];
  psi_tmp_imaginary = new double [nofstates];
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0;
    psi_tmp_imaginary[i] = 0;
  }

  average_energy=0.;
  /* initialize the  matrix. We have two kind of Hamiltonian operator here.
    First is the single spin operator, ss_operator.
    Second is the double spin operator, ds_operator.
    Ref. Article: Computational Methods for Simulating Quantum Computers euqation (68) & equation (70).
  */
  ss_operator_real      = new double [4];
  ss_operator_imaginary = new double [4];
  ds_operator_real      = new double [8];
  ds_operator_imaginary = new double [8];
  for (int i = 0; i < 4; i++) {
    ss_operator_real[i]      = 0;
    ss_operator_imaginary[i] = 0;
  }
  for (int i = 0; i < 8; i++) {
    ds_operator_real[i]      = 0;
    ds_operator_imaginary[i] = 0;
  }



  /* initialize the Hamiltonian matrix for calculating the energy */
  // H_real = new double [nofstates*(nofstates+1)/2];
  // H_imaginary = new double [nofstates*(nofstates+1)/2];
  // for (int i = 0; i < nofstates*(nofstates+1)/2; i++) {
  //   H_real[i] = 0.;
  //   H_imaginary[i]= 0.;
  // }




  Gamma=0; //time evolution of the initial Hamiltonian. Goes from 1 to 0
  Delta=0; //time evolution of the desired Hamiltonian. Goes from 0 to 1

}


/*
  Set up the initial basis state
  Input:
    d: define which mode to set up
      "read": read "psi_real.dat" and "psi_imagine_dat"
      "1upRan": set spin_0 as up and others are random.
      "allx": as the ground state of sigma_x.
*/
void spin_system::generate_initial_sys_state(char const * d){


  if ("read"==d){
    read(nofstates,psi_real,"psi_real.dat");
    read(nofstates,psi_imaginary,"psi_imagine.dat");
    cout<<"RRR"<<endl;
  }
  else if ("1upRand"==d) {
    srand (time(NULL));
    double normalize_factor=0.;
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
      if(i%2==0){
        psi_real[i]      = (double) rand()/RAND_MAX;//pow(nofstates, -0.5);
        psi_imaginary[i] = (double) rand()/RAND_MAX;
        normalize_factor += psi_real[i]*psi_real[i] + psi_imaginary[i]*psi_imaginary[i];
      } else {
        psi_real[i]      = 0;
        psi_imaginary[i] = 0;
      }
    }
    ofstream Psi_r_out("psi_real.dat");
    ofstream Psi_i_out("psi_imagine.dat");
    normalize_factor = sqrt(normalize_factor);
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
      if(i%2==0){
        psi_real[i]      =psi_real[i]/normalize_factor;
        psi_imaginary[i] =psi_imaginary[i]/normalize_factor;
      }
      Psi_r_out<<psi_real[i]<<endl;
      Psi_i_out<<psi_imaginary[i]<<endl;
    }

  }
  else if ("allx"==d){
    ofstream Psi_r_out("psi_real.dat");
    ofstream Psi_i_out("psi_imagine.dat");
    for (int i = 0; i < (int) pow(2,N_sys); i++) {
        psi_real[i]      =pow((int) pow(2,N_sys), -0.5);
        psi_imaginary[i] =0;
        Psi_r_out<<psi_real[i]<<endl;
        Psi_i_out<<psi_imaginary[i]<<endl;
    }
  }

  else {
    cout<<endl<<"Wrong Parameter for function: generate_initial_state(char)"<<endl;
    cout<<"!!!BE AWARE of the correctness of the result!!!"<<endl<<endl;

  }

}
/*
  to operate sigma_x,_y,_z*h
  Input:
    t: current time
*/
void spin_system::single_spin_op(double t){

  // #pragma omp parallel for
  for (int k = 0; k < N; k++) {
    int i1=(int) pow(2,k);

    /* update the single spin Hamiltonian matrix with t.
      In principal, we can save some computing time here,
      because of reading convenience I didn't do so.
    */
    double norm=0;
    //set the initial transverse field
    if (k<N_sys) {
      norm=sqrt((Gamma*h_x_start+Delta*h_x[k])*(Gamma*h_x_start+Delta*h_x[k])+Delta*h_y[k]*Delta*h_y[k]+Delta*h_z[k]*Delta*h_z[k]);
    } else {
      norm=0.;
    }



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
      ss_operator_real[1]      = Delta*h_y[k]*b;//sin(tau*0.5*norm)/norm;
      ss_operator_real[2]      = -1*Delta*h_y[k]*b;//sin(tau*0.5*norm)/norm;
      ss_operator_real[3]      = a;//cos(tau*0.5*norm);
      ss_operator_imaginary[0] = Delta*h_z[k]*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[1] = (Gamma*h_x_start+Delta*h_x[k])*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[2] = (Gamma*h_x_start+Delta*h_x[k])*b;//sin(tau*0.5*norm)/norm;
      ss_operator_imaginary[3] = -1*Delta*h_z[k]*b;//sin(tau*0.5*norm)/norm;
    }
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



/*
  To operate sigma_x*sigma_x
  Input:
    t: the current time point
*/
void spin_system::double_spin_op_x(double t){
  // #pragma omp parallel for
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jx_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jx_se[k+(l-N_sys)*N_sys];
      } else {
        J=J_x[k+l*N_sys];
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

        double b_block=b*Delta*tau*0.5;
        double c_block=c*Delta*tau*0.5;

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
/*
  To operate sigma_y*sigma_y
  Input:
    t: the current time point
*/
void spin_system::double_spin_op_y(double t){
  // #pragma omp parallel for
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jy_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jy_se[k+(l-N_sys)*N_sys];
      } else {
        J=J_y[k+l*N_sys];
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

        double b_block=b*Delta*tau*0.5;
        double c_block=c*Delta*tau*0.5;

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


/*
  To operate sigma_z*sigma_z
  Input:
    t: the current time point
*/
void spin_system::double_spin_op_z(double t){

  // #pragma omp parallel for default(none)
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=0.;
      if (k>=N_sys) {
        J=Jz_env[(k-N_sys)+(l-N_sys)*N_env];
      } else if(l>=N_sys && k<N_sys) {
        J=Jz_se[k+(l-N_sys)*N_sys];
      } else {
        J=J_z[k+l*N_sys];
      }
  // for (int i = 0; i < count_z; i++) {
    // int k=Jz_k_marked[i];
    // int l=Jz_l_marked[i];
      if(abs(J)>1e-15){


        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=J/1.;//4.;
        double b=0;//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=0;//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double a_block=a*Delta*tau;
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
          // if(3==ch)
            // cout<<"k, l, n0, n1, n2, n3 = "<<k<<", "<<l<<", "<<n0<<", "<<n1<<", "<<n2<<", "<<n3<<endl;
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

/*
  calculate energy w/o construct matrix
  Input:
    t: the current time point

*/
void spin_system::energy(double t){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  average_energy=0.;
  double check_img=0.;

  double hx=-1*h_x_start;
  for (int k = 0; k < N_sys; k++) {
    int i1=(int) pow(2,k);
    double hz=-1*h_z[k];
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
}

/*
  calculate the average spin
  Input:
    d: 'x' for sigma_x, 'y' for sigma_y, 'z' for sigma_z
    which_spin: which spin we want to calculate (count from 0)
*/
double spin_system::spin(char d,int which_spin){

  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  average_spin=0.;
  double check_img=0.;

  for (int k = 0; k < N_sys; k++) {
    if (k==which_spin) {
      int i1=(int) pow(2,k);
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
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }

  if (abs(check_img)>1e-13)
    cout<<"Something went wrong in function spin()"<<check_img<<endl;
  return average_spin;


}

/*
  Read in the number from a txt into an array
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


/*
  set up an environment by diagonalize the H matrix.
  Input:
    int: the # of spins.
*/
void spin_system::environment(int N, double Temperature){


  //Construct the H_env matrix for the later use for solving energy of the heat bath
  double Boltzmann_const= 1;//we use 1 here insted of 0.000086173303 ev/K;
  int nofstates=(int) pow(2,N);
  double* H_env;
  H_env=new double[nofstates*(nofstates+1)/2];
  /* // Since I don't set single spin magnetization for environment, I comment this part now
  for (int k = 0; k < N; k++) {
    int i1=(int) pow(2,k);
    double h_x_init=0;
    double h_z_tmp =0;

    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l - i2 +i2/i1;
      int j = i+i1;
        H[i+i*(i+1)/2] += 1.*Delta*h_z_tmp;
        H[j+j*(j+1)/2] += -1.*Delta*h_z_tmp;
        H[i+j*(j+1)/2] += 1.*Gamma*h_x_init;
    }
  }*/
  for (int k = 0; k < N; k++) {
    for (int l = k+1; l < N; l++) {
      double Jx=-1*Jx_env[k+l*N];
      double Jy=-1*Jy_env[k+l*N];
      double Jz=-1*Jz_env[k+l*N];
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
  double w[nofstates];
  int n=nofstates;
  double vl,vu;
  int il,iu;
  double abstol = 2*DLAMCH("S");
  int m;
  complex<double> ap[nofstates*(nofstates+1)/2];
  for (int i = 0; i < nofstates*(nofstates+1)/2; i++)
    ap[i]=H_env[i];
  complex<double> z[nofstates*nofstates];
  int lda = nofstates;
  complex<double> work[2*nofstates];
  double rwork[7*nofstates];
  int iwork[5*nofstates];
  int ifail[nofstates];
  int info;
  zhpevx("Vector","All","Upper", &n, ap, &vl, &vu, &il, &iu, &abstol, &m, w, z, &lda, work, rwork, iwork, ifail, &info );
  double sum=0.;
  // cout<<w[0]<<endl;
  for (int i = 1; i < nofstates; i++) {
    // cout<<w[i]<<endl;
    w[i]=exp(-1*(w[i]-w[0])/(Temperature*Boltzmann_const));
    sum+=w[i];
  }
  w[0]=1;//exp(-1*(w[0]-w[0])/(Temperature*Boltzmann_const));
  sum+=w[0];
  cout<<"sum of z: "<<sum<<endl;
  for (int i = 0; i < nofstates; i++) {
    w[i]=w[i]/sum;
  }
  sum=0;
  for (int i = 0; i < nofstates; i++) {
    sum+=w[i];
  }
  cout<<"inside function Environment: "<<sum<<endl;
  ofstream envr_out("env_real.dat");
  ofstream envi_out("env_imag.dat");
  for (int i = 0; i < nofstates; i++) {
    envr_out<<w[i]<<endl;
    envi_out<<0<<endl;
  }

}


/*
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
  cout<<"sum inside generate functoin: "<<sum<<endl;
}
/*
  Make J_env txt file randomly between [-1,1]
  It is the coupling factor for environment spins

  Input:
    N: # of spins
    G: strength
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



/*
  Make J_se txt file randomly between [-1,1]
  and then multiply with a global interaction factor G,
  control how strong the system interact with the environment.
  Input:
    N: # of spins= N_env+N_sys
    G: stenght factor
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
*/
void spin_system::run(){
  // ofstream out_data("../Result/product_formula/gs_T_1e100_t_1e-1.dat");
  // ofstream E_out("../Result/Verify/spin_x4_pd.dat");
  ofstream output("output.dat");
  output<<"Time Energy ";
  for (int i = 0; i < N_sys; i++) {
    output<<"Sx_"<<i<<" "<<"Sy_"<<i<<" "<<"Sz_"<<i<<" ";
  }
  output<<endl;
  int total_steps=0;
  total_steps=(int) T/tau;

  /////test not go over whole J again and again/////
  //////////////////////////////////////////////////

  count_z=0;
  count_y=0;
  count_x=0;
  Jz_k_marked=new double [N_sys*(N_sys+1)/2];
  Jz_l_marked=new double [N_sys*(N_sys+1)/2];
  Jy_k_marked=new double [N_sys*(N_sys+1)/2];
  Jy_l_marked=new double [N_sys*(N_sys+1)/2];
  Jx_k_marked=new double [N_sys*(N_sys+1)/2];
  Jx_l_marked=new double [N_sys*(N_sys+1)/2];
  for (int k = 0; k <N_sys ; k++) {
    for (int l = k+1; l < N_sys; l++) {
      if (abs(J_z[k+l*N_sys])>1e-15){
        Jz_k_marked[count_z]=k;
        Jz_l_marked[count_z]=l;
        count_z+=1;
      }
      if (abs(J_y[k+l*N_sys])>1e-15){
        Jy_k_marked[count_y]=k;
        Jy_l_marked[count_y]=l;
        count_y+=1;
      }
      if (abs(J_x[k+l*N_sys])>1e-15){
        Jx_k_marked[count_x]=k;
        Jx_l_marked[count_x]=l;
        count_x+=1;
      }
    }
  }

  cout<<"Elements in Jx, Jy, and Jz: "<<count_x<<" "<<count_y<<" "<<count_z<<" "<<endl;
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  for (int step = 0; step < total_steps+1; step++){
    // for (int i = 0; i < nofstates*(nofstates+1)/2; i++) {
    //   H_real[i]=0.;
    //   H_imaginary[i]=0.;
    // }
    Delta=step*tau/T;
    Gamma=1-Delta;

    energy(step*tau);
    output<<step*tau<<" "<<average_energy<<" ";
    for (int s = 0; s < N_sys; s++) {
      output<<spin('x',s)<<" "<<spin('y',s)<<" "<<spin('z',s)<<" ";
    }
    output<<endl;
    single_spin_op(step*tau);
    double_spin_op_x(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_z(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_x(step*tau);
    single_spin_op(step*tau);


    double gs=0.;

    // E_out<<step*tau/T<<" "<<energy_Hmatrix<<" "<<average_energy<<endl;
    // E_out<<step*tau<<" "<<average_spin<<endl;


    // gs=psi_real[176]*psi_real[176]+psi_imaginary[176]*psi_imaginary[176];
    // out_data<<step*tau/T<<" "<<gs<<endl;
    // gs =0.;

  }
}
