//#include "spinsys.h" test
#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

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
  double Gamma;
  double Delta;

public:
  void initialize(int,double,double);
  void test();
  void single_spin_op(double, int);
  void double_spin_op(double, int);
  void run(int);
  void read(int,double*,char const *);
  double* psi_real;
  double* psi_imaginary;
};

/* Some Notes for myself: */
/* Don't use class complex<double>. It is terribly slow */

int main(int argc, char* argv[]){
spin_system test;
test.initialize(5,10,0.018);
for (int i = 0; i < (int) pow(2,5); i++) {
  cout<<i<<"th real"<<test.psi_real[i]<<endl;
  cout<<i<<"th imag"<<test.psi_imaginary[i]<<endl;
}

test.run(atoi(argv[1]));

for (int i = 0; i < (int) pow(2,5); i++) {
  cout<<i<<"th real after ope "<<test.psi_real[i]<<endl;
  cout<<i<<"th imag after ope "<<test.psi_imaginary[i]<<endl;
}
// test.test();
double* array;
array = new double [25];
test.read(25,array,"J.txt");
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
    psi_real[i]      = pow(nofstates,-0.5);
    psi_imaginary[i] = 0;
  }

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
  J_x = new double [N*(N-1)];
  J_y = new double [N*(N-1)];
  J_z = new double [N*(N-1)];
  h_x = new double [N];
  h_y = new double [N];
  h_z = new double [N];
  for (int i = 0; i < N*(N-1); i++){
    J_x[i] = 0.;
    J_y[i] = 0.;
    J_z[i] = 0.;
  }

  for (int i = 0; i < N; i++){
    h_x[i]=0.;
    h_y[i]=0.;
    h_z[i]=0.;
  }

  /* Other variables:
  */
    Gamma=1; //time evolution of the initial Hamiltonian. Goes from 1 to 0
    Delta=0; //time evolution of the desired Hamiltonian. Goes from 0 to 1

  /* The following part is in doubr
  */
  // hz=0;
  // E=0;
  // J=0;

}


/* to operate sigma_x,_y,_z*h
*/
void spin_system::single_spin_op(double t, int t_on){

  for (int k = 0; k < N; k++) {
    int i1=(int) pow(2,k);

    /* update the single spin Hamiltonian matrix with t.
      In principal, we can save some computing time here,
      because of reading convenience I didn't do so.
    */
    double norm=0;
    //set the initial transverse field
    double h_x_init=0.;
    if (!t_on) {
      //without t dependent
      norm=sqrt(h_x[i1]*h_x[i1]+h_y[i1]*h_y[i1]+h_z[i1]*h_z[i1]);
      // cout<<"t_on = "<<t_on<<"off"<<endl;
    }
    else {
      //with t dependent
      Delta = t/T;
      norm=sqrt(((1-Delta)*h_x_init+Delta*h_x[i1])*((1-Delta)*h_x_init+Delta*h_x[i1])+Delta*h_y[i1]*Delta*h_y[i1]+Delta*h_z[i1]*Delta*h_z[i1]);
      // cout<<"t_on = "<<t_on<<"on"<<endl;
    }


    if (norm-0<1e-15) {
      ss_operator_real[0]      = cos(tau*norm*0.5);
      ss_operator_real[1]      = 0;
      ss_operator_real[2]      = 0;
      ss_operator_real[3]      = cos(tau*norm*0.5);
      ss_operator_imaginary[0] = 0;
      ss_operator_imaginary[0] = 0;
      ss_operator_imaginary[0] = 0;
      ss_operator_imaginary[0] = 0;
    }

    else {

      ss_operator_real[0]      = cos(tau*norm*0.5);
      ss_operator_real[1]      = Delta*h_y[i1]*sin(tau*norm*0.5)/norm;
      ss_operator_real[2]      = -1*Delta*h_y[i1]*sin(tau*norm*0.5)/norm;
      ss_operator_real[3]      = cos(tau*norm*0.5);
      ss_operator_imaginary[0] = Delta*h_z[i1]*sin(tau*norm*0.5)/norm;
      ss_operator_imaginary[0] = ((1-Delta)*h_x_init+Delta*h_x[i1])*sin(tau*norm*0.5)/norm;
      ss_operator_imaginary[0] = ((1-Delta)*h_x_init+Delta*h_x[i1])*sin(tau*norm*0.5)/norm;
      ss_operator_imaginary[0] = -1*Delta*h_z[i1]*sin(tau*norm*0.5)/norm;
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


/* to operate J*S.*S
*/
void spin_system::double_spin_op(double t, int t_on){
  for (int k = 0; k <N ; k++) {
    for (int l = k; l < N; l++) {

      /* update the double spin Hamiltonian matrix with t.
        In principal, we can save some computing time here,
        because of reading convenience I didn't do so.
      */
      double a=J_z[k+l*N]/4;
      double b=(J_x[k+l*N]-J_y[k+l*N])/4;
      double c=(J_x[k+l*N]+J_y[k+l*N])/4;

      Delta = t/T;

      ds_operator_real[0]      = cos( a*Delta*tau)*cos(b*Delta*tau);
      ds_operator_real[1]      =-sin( a*Delta*tau)*sin(b*Delta*tau);
      ds_operator_real[2]      = cos(-a*Delta*tau)*cos(c*Delta*tau);
      ds_operator_real[3]      =-sin(-a*Delta*tau)*sin(c*Delta*tau);
      ds_operator_real[4]      =-sin(-a*Delta*tau)*sin(c*Delta*tau);
      ds_operator_real[5]      = cos(-a*Delta*tau)*cos(c*Delta*tau);
      ds_operator_real[6]      =-sin( a*Delta*tau)*sin(b*Delta*tau);
      ds_operator_real[7]      = cos( a*Delta*tau)*cos(b*Delta*tau);
      ds_operator_imaginary[0] = sin( a*Delta*tau)*cos(b*Delta*tau);
      ds_operator_imaginary[1] = cos( a*Delta*tau)*sin(b*Delta*tau);
      ds_operator_imaginary[2] = sin(-a*Delta*tau)*cos(c*Delta*tau);
      ds_operator_imaginary[3] = cos(-a*Delta*tau)*sin(c*Delta*tau);
      ds_operator_imaginary[4] = cos(-a*Delta*tau)*sin(c*Delta*tau);
      ds_operator_imaginary[5] = sin(-a*Delta*tau)*cos(c*Delta*tau);
      ds_operator_imaginary[6] = cos( a*Delta*tau)*sin(b*Delta*tau);
      ds_operator_imaginary[7] = sin( a*Delta*tau)*cos(b*Delta*tau);

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

        /* Following the similar manner in singli_spin_op,
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
    // cout<<Array[i]<<endl;
  }
}


/* The main process to run the simulation
  16.12.2016: I haven't add the time evolution part. It should be added after.
*/
void spin_system::run(int t_on){
  for(double t=0. ; t<T ; t+=tau){
    single_spin_op(t,t_on);
    double_spin_op(t,t_on);
    single_spin_op(t,t_on);
    // cout<<t<<endl;
  }
}


void spin_system::test(){
  cout<<nofstates<<endl;
  cout<<tau<<endl;
}
