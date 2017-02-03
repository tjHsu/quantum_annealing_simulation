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
#include <time.h>


using namespace std;
class spin_system {
private:
  int N;
  double T;
  double tau;
  int nofstates;

  double* ss_operator_real;
  double* ss_operator_imaginary;
  double* ds_operator_real;
  double* ds_operator_imaginary;
  double* J_x;
  double* J_y;
  double* J_z;
  double* h_x;
  double* h_y;
  double* h_z;
  double  h_x_start;
  double Gamma;
  double Delta;


public:
  void initialize(int ,double ,double);
  void single_spin_op(double);
  void double_spin_op(double);
  void double_spin_op_x(double);
  void double_spin_op_y(double);
  void double_spin_op_z(double);
  void generate_initial_state(char const *);
  void energy(double);
  void spin(char,int);
  void run();
  void read(int,double*,char const *);
  double* psi_real;
  double* psi_imaginary;
  double* psi_tmp_real;
  double* psi_tmp_imaginary;
  double* H_real;
  double* H_imaginary;
  double average_energy;
  double average_spin;


};

int main(int argc, char* argv[]){

  spin_system test;

  test.initialize(8,1000.,0.01);


  double norm=0.;
  for (int i = 0; i < (int) pow(2,8); i++) {
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
  for (int i = 0; i < (int) pow(2,8); i++) {
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

  ofstream Coefficient_out("coefficient_out_pd.dat");
  for (int i = 0; i < (int) pow(2,8); i++) {
    double tmp=test.psi_real[i]*test.psi_real[i]+test.psi_imaginary[i]*test.psi_imaginary[i];
    if(tmp<1e-3){
      Coefficient_out<<0<<endl;
    } else {
    Coefficient_out<<tmp<<endl;
    }
  }

  cout<<"reutrn "<<0<<endl;
  return 0;
}


/* initialize the system
  N: number of spins
  T: end time
  tau: the size of the timestep
*/
void spin_system::initialize(int N_user_defined, double T_user_defined, double tau_user_defined){
  N = N_user_defined;
  T = T_user_defined;
  tau = tau_user_defined;
  nofstates = (int) pow(2,N);// number of states construct by the number of spins

  /* initialize the wave function in the ground state
  */

  psi_real = new double [nofstates];
  psi_imaginary = new double [nofstates];


  for (int i = 0; i < nofstates; i++) {
      psi_real[i]      = 0;
      psi_imaginary[i] = 0;
  }


  generate_initial_state("allx");




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

  /* initialize the coupling factor J for double spin Hamiltonian,
    and the factor h for single spin Hamiltonian.
    J should have size as N*(N-1) and the size of h is simply N.
  */
  J_x = new double [N*N];
  J_y = new double [N*N];
  J_z = new double [N*N];
  h_x = new double [N];
  h_y = new double [N];
  h_z = new double [N];
  h_x_start=1;
  for (int i = 0; i < N*N; i++){
    J_x[i] = 0.;
    J_y[i] = 0.;
    J_z[i] = 0.;
  }

  for (int i = 0; i < N; i++){
    h_x[i]=0.;
    h_y[i]=0.;
    h_z[i]=0.;
  }

  /* initialize the Hamiltonian matrix for calculating the energy */
  H_real = new double [nofstates*(nofstates+1)/2];
  H_imaginary = new double [nofstates*(nofstates+1)/2];
  for (int i = 0; i < nofstates*(nofstates+1)/2; i++) {
    H_real[i] = 0.;
    H_imaginary[i]= 0.;
  }

  read(N*N,J_z,"J2.txt");
  read(N*N,J_x,"J2x.txt");
  read(N*N,J_y,"J2y.txt");
  read(N,h_z,"h2.txt");


  Gamma=1; //time evolution of the initial Hamiltonian. Goes from 1 to 0
  Delta=0; //time evolution of the desired Hamiltonian. Goes from 0 to 1

}

void spin_system::generate_initial_state(char const * d){


  if ("read"==d){
    read(nofstates,psi_real,"psi_real.dat");
    read(nofstates,psi_imaginary,"psi_imagine.dat");
    cout<<"RRR"<<endl;
  }
  else if ("1upRand"==d) {
    srand (time(NULL));
    double normalize_factor=0.;
    for (int i = 0; i < nofstates; i++) {
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
    for (int i = 0; i < nofstates; i++) {
      if(i%2==0){
        psi_real[i]      =psi_real[i]/normalize_factor;
        psi_imaginary[i] =psi_imaginary[i]/normalize_factor;
      }
      Psi_r_out<<psi_real[i]<<endl;
      Psi_i_out<<psi_imaginary[i]<<endl;
    }

  }
  else if ("allx"==d){
    for (int i = 0; i < nofstates; i++) {
        psi_real[i]      =pow(nofstates, -0.5);
        psi_imaginary[i] =0;
    }
  }
  else {
    cout<<endl<<"Wrong Parameter for function: generate_initial_state(char)"<<endl;
    cout<<"!!!BE AWARE of the correctness of the result!!!"<<endl<<endl;

  }

}
/* to operate sigma_x,_y,_z*h
*/
void spin_system::single_spin_op(double t){

  for (int k = 0; k < N; k++) {
    int i1=(int) pow(2,k);

    /* update the single spin Hamiltonian matrix with t.
      In principal, we can save some computing time here,
      because of reading convenience I didn't do so.
    */
    double norm=0;
    //set the initial transverse field

    norm=sqrt((Gamma*h_x_start+Delta*h_x[k])*(Gamma*h_x_start+Delta*h_x[k])+Delta*h_y[k]*Delta*h_y[k]+Delta*h_z[k]*Delta*h_z[k]);


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

  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=J_x[k+l*N];
      if(abs(J)>1e-15){

        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=0;//J_z[k+l*N]/1.;//4.;
        double b=J_x[k+l*N];//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=J_x[k+l*N];//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double b_block=b*Delta*tau*0.5;
        double c_block=c*Delta*tau*0.5;



        ds_operator_real[0]      = cos(b_block);//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[2]      = cos(c_block);//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[5]      = cos(c_block);//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[7]      = cos(b_block);//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[1] = sin(b_block);//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[3] = sin(c_block);//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[4] = sin(c_block);//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[6] = sin(b_block);//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);

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

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] + ds_operator_real[1]*psi_real[n3]
                                  -ds_operator_imaginary[0]*psi_imaginary[n0] - ds_operator_imaginary[1]*psi_imaginary[n3];
          psi_imaginary_temp_n0 = ds_operator_imaginary[0]*psi_real[n0] + ds_operator_imaginary[1]*psi_real[n3]
                                  +ds_operator_real[0]*psi_imaginary[n0] + ds_operator_real[1]*psi_imaginary[n3];

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] + ds_operator_real[3]*psi_real[n2]
                                  -ds_operator_imaginary[2]*psi_imaginary[n1] - ds_operator_imaginary[3]*psi_imaginary[n2];
          psi_imaginary_temp_n1 = ds_operator_imaginary[2]*psi_real[n1] + ds_operator_imaginary[3]*psi_real[n2]
                                  +ds_operator_real[2]*psi_imaginary[n1] + ds_operator_real[3]*psi_imaginary[n2];

          psi_real_temp_n2      = ds_operator_real[4]*psi_real[n1] + ds_operator_real[5]*psi_real[n2]
                                  -ds_operator_imaginary[4]*psi_imaginary[n1] - ds_operator_imaginary[5]*psi_imaginary[n2];
          psi_imaginary_temp_n2 = ds_operator_imaginary[4]*psi_real[n1] + ds_operator_imaginary[5]*psi_real[n2]
                                  +ds_operator_real[4]*psi_imaginary[n1] + ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = ds_operator_real[6]*psi_real[n0] + ds_operator_real[7]*psi_real[n3]
                                  -ds_operator_imaginary[6]*psi_imaginary[n0] - ds_operator_imaginary[7]*psi_imaginary[n3];
          psi_imaginary_temp_n3 = ds_operator_imaginary[6]*psi_real[n0] + ds_operator_imaginary[7]*psi_real[n3]
                                  +ds_operator_real[6]*psi_imaginary[n0] + ds_operator_real[7]*psi_imaginary[n3];

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

  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=J_y[k+l*N];
      if(abs(J)>1e-15){

        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=0;//J_z[k+l*N]/1.;//4.;
        double b=-J_y[k+l*N];//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=J_y[k+l*N];//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;

        double b_block=b*Delta*tau*0.5;
        double c_block=c*Delta*tau*0.5;



        ds_operator_real[0]      = cos(b_block);//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[2]      = cos(c_block);//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_real[5]      = cos(c_block);//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_real[7]      = cos(b_block);//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        ds_operator_imaginary[1] = sin(b_block);//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[3] = sin(c_block);//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[4] = sin(c_block);//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        ds_operator_imaginary[6] = sin(b_block);//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);

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

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] + ds_operator_real[1]*psi_real[n3]
                                  -ds_operator_imaginary[0]*psi_imaginary[n0] - ds_operator_imaginary[1]*psi_imaginary[n3];
          psi_imaginary_temp_n0 = ds_operator_imaginary[0]*psi_real[n0] + ds_operator_imaginary[1]*psi_real[n3]
                                  +ds_operator_real[0]*psi_imaginary[n0] + ds_operator_real[1]*psi_imaginary[n3];

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] + ds_operator_real[3]*psi_real[n2]
                                  -ds_operator_imaginary[2]*psi_imaginary[n1] - ds_operator_imaginary[3]*psi_imaginary[n2];
          psi_imaginary_temp_n1 = ds_operator_imaginary[2]*psi_real[n1] + ds_operator_imaginary[3]*psi_real[n2]
                                  +ds_operator_real[2]*psi_imaginary[n1] + ds_operator_real[3]*psi_imaginary[n2];

          psi_real_temp_n2      = ds_operator_real[4]*psi_real[n1] + ds_operator_real[5]*psi_real[n2]
                                  -ds_operator_imaginary[4]*psi_imaginary[n1] - ds_operator_imaginary[5]*psi_imaginary[n2];
          psi_imaginary_temp_n2 = ds_operator_imaginary[4]*psi_real[n1] + ds_operator_imaginary[5]*psi_real[n2]
                                  +ds_operator_real[4]*psi_imaginary[n1] + ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = ds_operator_real[6]*psi_real[n0] + ds_operator_real[7]*psi_real[n3]
                                  -ds_operator_imaginary[6]*psi_imaginary[n0] - ds_operator_imaginary[7]*psi_imaginary[n3];
          psi_imaginary_temp_n3 = ds_operator_imaginary[6]*psi_real[n0] + ds_operator_imaginary[7]*psi_real[n3]
                                  +ds_operator_real[6]*psi_imaginary[n0] + ds_operator_real[7]*psi_imaginary[n3];

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
  // Delta = t/T;
  // cout<<t<<endl;
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double J=J_z[k+l*N];
      // if((J==0)) cout<<" equal ZERO!! "<<k<<" "<<l<<endl;
      if(abs(J)>1e-15){


        /* update the double spin Hamiltonian matrix with t.
          In principal, we can save some computing time here,
          because of reading convenience I didn't do so.
        */

        double a=J_z[k+l*N]/1.;//4.;
        double b=0;//(J_x[k+l*N]-J_y[k+l*N])/1.;//4.;
        double c=0;//(J_x[k+l*N]+J_y[k+l*N])/1.;//4.;




        ds_operator_real[0]      = cos( a*Delta*tau);//cos( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_real[1]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_real[2]      = cos(-a*Delta*tau);//cos(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_real[3]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_real[4]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_real[5]      = cos(-a*Delta*tau);//cos(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_real[6]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_real[7]      =  cos(a*Delta*tau);//cos( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_imaginary[0] = sin( a*Delta*tau);//sin( a*Delta*tau)*cos(b*Delta*tau);
        ds_operator_imaginary[1] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_imaginary[2] = sin(-a*Delta*tau);//sin(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_imaginary[3] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_imaginary[4] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        ds_operator_imaginary[5] = sin(-a*Delta*tau);//sin(-a*Delta*tau)*cos(c*Delta*tau);
        ds_operator_imaginary[6] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        ds_operator_imaginary[7] = sin( a*Delta*tau);//sin( a*Delta*tau)*cos(b*Delta*tau);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);

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

          psi_real_temp_n0      = ds_operator_real[0]*psi_real[n0] + ds_operator_real[1]*psi_real[n3]
                                  -ds_operator_imaginary[0]*psi_imaginary[n0] - ds_operator_imaginary[1]*psi_imaginary[n3];
          psi_imaginary_temp_n0 = ds_operator_imaginary[0]*psi_real[n0] + ds_operator_imaginary[1]*psi_real[n3]
                                  +ds_operator_real[0]*psi_imaginary[n0] + ds_operator_real[1]*psi_imaginary[n3];

          psi_real_temp_n1      = ds_operator_real[2]*psi_real[n1] + ds_operator_real[3]*psi_real[n2]
                                  -ds_operator_imaginary[2]*psi_imaginary[n1] - ds_operator_imaginary[3]*psi_imaginary[n2];
          psi_imaginary_temp_n1 = ds_operator_imaginary[2]*psi_real[n1] + ds_operator_imaginary[3]*psi_real[n2]
                                  +ds_operator_real[2]*psi_imaginary[n1] + ds_operator_real[3]*psi_imaginary[n2];

          psi_real_temp_n2      = ds_operator_real[4]*psi_real[n1] + ds_operator_real[5]*psi_real[n2]
                                  -ds_operator_imaginary[4]*psi_imaginary[n1] - ds_operator_imaginary[5]*psi_imaginary[n2];
          psi_imaginary_temp_n2 = ds_operator_imaginary[4]*psi_real[n1] + ds_operator_imaginary[5]*psi_real[n2]
                                  +ds_operator_real[4]*psi_imaginary[n1] + ds_operator_real[5]*psi_imaginary[n2];

          psi_real_temp_n3      = ds_operator_real[6]*psi_real[n0] + ds_operator_real[7]*psi_real[n3]
                                  -ds_operator_imaginary[6]*psi_imaginary[n0] - ds_operator_imaginary[7]*psi_imaginary[n3];
          psi_imaginary_temp_n3 = ds_operator_imaginary[6]*psi_real[n0] + ds_operator_imaginary[7]*psi_real[n3]
                                  +ds_operator_real[6]*psi_imaginary[n0] + ds_operator_real[7]*psi_imaginary[n3];

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
  for (int k = 0; k < N; k++) {
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

  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      double Jx=-1*J_x[k+l*N];
      double Jy=-1*J_y[k+l*N];
      double Jz=-1*J_z[k+l*N];
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
  if (abs(check_img)>1e-15)
    cout<<"Something went wrong in functoin energy() "<<check_img<<endl;
}

/*
  calculate the average spin
  Input:
    d: 'x' for sigma_x, 'y' for sigma_y, 'z' for sigma_z
    which_spin: which spin we want to calculate (count from 0)
*/
void spin_system::spin(char d,int which_spin){

  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]      = 0.;
    psi_tmp_imaginary[i] = 0.;
  }
  average_spin=0.;
  double check_img=0.;

  for (int k = 0; k < N; k++) {
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
          cout<<"run input for calculating average spin"<<endl;
        }

      }
    }
  }

  for (int i = 0; i < nofstates; ++i) {
    // cout<<2++1*3<<endl;
    average_spin += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  if (abs(check_img)>1e-15)
    cout<<"Something went wrong in function spin()"<<check_img<<endl;



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


/* The main process to run the simulation
  16.12.2016: I haven't add the time evolution part. It should be added after.
  26.12.2016: time evolution part added.
  03.02.2017: can handle with sigma_(x,y,z) and sigma_(x,y,z)*sigma_(x,y,z) respectly.
*/
void spin_system::run(){
  ofstream out_data("../Result/product_formula/gs_T_1e100_t_1e-1.dat");
  ofstream E_out("../Result/Verify/spin_x4_pd.dat");

  int total_steps=0;
  total_steps=(int) T/tau;

  for (int step = 0; step < total_steps+1; step++){
    for (int i = 0; i < nofstates*(nofstates+1)/2; i++) {
      H_real[i]=0.;
      H_imaginary[i]=0.;
    }
    Delta=step*tau/T;
    Gamma=1-Delta;
    // spin('x',0);
    single_spin_op(step*tau);
    double_spin_op_x(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_z(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_x(step*tau);
    single_spin_op(step*tau);
    // energy(step*tau);

    double gs=0.;

//     E_out<<step*tau/T<<" "<<energy_Hmatrix<<" "<<average_energy<<endl;
    E_out<<step*tau<<" "<<average_spin<<endl;


    gs=psi_real[176]*psi_real[176]+psi_imaginary[176]*psi_imaginary[176];
    out_data<<step*tau/T<<" "<<gs<<endl;
    gs =0.;

  }
}
