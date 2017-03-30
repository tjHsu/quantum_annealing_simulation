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
#include <random>
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
  double Temperature;
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

  /* optimizing by skipping the zero term.*/
  // double* h_k_marked;
  // int count_h;
  double* Jz_k_marked;
  double* Jz_l_marked;
  double* Jz_J_marked;
  int count_z;
  double* Jy_k_marked;
  double* Jy_l_marked;
  double* Jy_J_marked;
  int count_y;
  double* Jx_k_marked;
  double* Jx_l_marked;
  double* Jx_J_marked;
  int count_x;
  double* Jz_k_eng_marked;
  double* Jz_l_eng_marked;
  double* Jz_J_eng_marked;
  int count_z_eng;
  double* J_k_combine_marked;
  double* J_l_combine_marked;
  double* Jx_J_combine_marked;
  double* Jy_J_combine_marked;
  double* Jz_J_combine_marked;
  int count_combine;

  int J_index;

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
  void exp_appr_op(double, int);
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
  void sumaverage(int, double*, double*, complex<double>*, double*, double*);

  void read(int,double*,char const *);



public:
  void skip_zeroterm();
  void initialize(int, int ,double ,double, double, double, int, int);
  void run();
  void random_wavef_run();

  double* spin_return;
  double* energy_sys_return;
  double* energy_env_return;
  double* energy_se_return;
  double* energy_all_return;
  double* coefficient_return;
  double  success_probability_return;
  int gs_sol;
};
