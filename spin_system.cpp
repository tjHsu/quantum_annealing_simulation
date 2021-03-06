#include "spin_system.h"

/* Initialize the annealing system with environment.
  Input:
    N_sys_user_defined: number of spins of system
    N_env_user_defined: number of spins of environment
    T_user_defined: The Annealing time
    tau_user_defined: the size of timestep of simulation
    G_user_defined: the coupling factor for Jse (system and environment)
    J_index_user_defined: If not specify, the program will use J4 and h4 now. otherwise pleasr define which set of J the program should use.
    env_on: if equal to 1, calculate the eigenvectors and eigenvalues of the environment spins using Lapack.
    random_approx_on: if equal to 1, prepare a random state which approximate exponential (-B*H_e/2m)^m
*/
void spin_system::initialize(int N_sys_user_defined, int N_env_user_defined, double T_user_defined, double tau_user_defined, double Temperature_user_defined, double G_user_defined, int env_on, int J_index_user_defined, int random_approx_on){
  Temperature = Temperature_user_defined;
  N_sys = N_sys_user_defined;
  N_env = N_env_user_defined;
  N = N_sys+N_env;
  nofstates = (int) pow(2,N);// number of states construct by the number of total spins
  T = T_user_defined;
  tau = tau_user_defined;
  J_index=J_index_user_defined;
  G = G_user_defined;

  /* initialize the coupling factor J for double spin Hamiltonian,
    and the factor h for single spin Hamiltonian.
    J is for system, J_env is for environment, and J_se is the interaction inbetween system and environment
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

  /*
  J_index ==0: for test case. Please chage target file here
  J_index !=0: for regular case. Please change from the for loop in main.cpp
  */
  if (0==J_index) {
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
    gs_sol=119;/*the solution ground state for this Hamiltonian.*/
    cout<<"For J_index = 0, solution state is set to : "<<gs_sol<<endl;
  } else {
    cout<<"For J_index != 0, read in file from the path: ";

    ostringstream index_strs;
    index_strs<<J_index;;
    string index_str=index_strs.str();
    cout<<index_str<<endl;

    string J_str("/home0/t.hsu/Documents/2016_WS_MT/8spin/J");
    string h_str("/home0/t.hsu/Documents/2016_WS_MT/8spin/h");
    string s_str("/home0/t.hsu/Documents/2016_WS_MT/8spin/s");
    // ostringstream strs;
    // strs<<"/home0/t.hsu/Documents/2016_WS_MT/8spin/J"<<(J_index)<<".txt";
    // string str=strs.str();

    string tmp_str;
    /*read in Jz*/
    tmp_str=J_str+index_str+".txt";
    const char *testChars = tmp_str.c_str();
    cout<<testChars<<endl;
    read(N_sys*N_sys,J_z,testChars);
    /*read in hz*/
    tmp_str=h_str+index_str+".txt";
    testChars=tmp_str.c_str();
    cout<<testChars<<endl;
    read(N_sys,h_z,testChars);

    // string teststr=to_string(1);
    /*read in ground state index*/
    tmp_str=s_str+index_str+".txt";
    testChars=tmp_str.c_str();
    cout<<testChars<<endl;
    ifstream myfile;
    myfile.open(testChars);
    if(!myfile.is_open()){
      cout<<"Unable to open the file: "<<testChars<<endl;
    }
    myfile>>gs_sol;
    myfile.close();

  }


  /*G is the global factor usually set between -1 and 1.*/
  // G=1.0;
  // Jenv_generate(N_env,G);/*randomly generate J_env*/
  // G=0.05;
  // Jse_generate(N_sys,N_env,G);/*randomly generate J_se*/

  if (1==env_on) {
    /* initialize the array for the enivironment's partition factor w[], and its eignevector in computational basis z[] */
    z= new complex<double> [(int)pow(2,N_env)*(int)pow(2,N_env)]();
    w= new double [(int)pow(2,N_env)]();
    environment(N_env,Temperature);/*get w[] and z[] with lapack diagonalization.*/
  }

  /* Change the coupling factor */
  for (int i = 0; i < N_sys*N_env; i++) {
    Jx_se[i]=Jx_se[i]*G;
    Jy_se[i]=Jy_se[i]*G;
    Jz_se[i]=Jz_se[i]*G;
  }




  /* initialize the wave function in the ground state */
  psi_real = new double [nofstates]();
  psi_imaginary = new double [nofstates]();
  psi_tmp_real = new double [nofstates]();
  psi_tmp_imaginary = new double [nofstates]();
  psi_sys_real = new double [(int)pow(2,N_sys)]();
  psi_sys_imaginary = new double [(int)pow(2,N_sys)]();

  /*this funtion can prepare system into some predefine state*/
  set_initial_sys_state("allx");

  /*initialize the psi function to the random wavefunction, if we use random wave function algorithm*/
  if (1==random_approx_on) {
    random_product(psi_real,psi_imaginary,psi_sys_real,psi_sys_imaginary);
    int M=10000;
    for (int i = 0; i < M; i++) {
      if (i%1000==0)
        cout<<".."<<i<<".."<<endl;
      exp_appr_op(M);
    }
    ofstream Psi_r_out("output_psi_appr_real.dat");
    ofstream Psi_i_out("output_psi_appr_imagine.dat");
    for (int i = 0; i < nofstates; i++) {
        Psi_r_out<<psi_real[i]<<endl;
        Psi_i_out<<psi_imaginary[i]<<endl;
    }

  } else if (0==random_approx_on) {
    read(nofstates,psi_real,"output_psi_appr_real.dat");
    read(nofstates,psi_imaginary,"output_psi_appr_imagine.dat");
    cout<<"Read for psi_real and psi_imaginary after exp_appr function"<<endl;
  } else {
    cout<<"Something go wrong with exp_appr_op() during Initialization";
  }

  /* initialize the  matrix. We have two kind of Hamiltonian operator here.
    First is the single spin operator, ss_operator.
    Second is the double spin operator, ds_operator.
    Ref. Article: Computational Methods for Simulating Quantum Computers euqation (68) & equation (70).
  */
  ss_operator_real      = new double [4]();
  ss_operator_imaginary = new double [4]();
  ds_operator_real      = new double [8]();
  ds_operator_imaginary = new double [8]();

  // /* uncommend if want to use spin_allinone(); */
  // psi_tmp_x_real=new double [nofstates];
  // psi_tmp_x_imaginary=new double [nofstates];
  // psi_tmp_y_real=new double [nofstates];
  // psi_tmp_y_imaginary=new double [nofstates];
  // psi_tmp_z_real=new double [nofstates];
  // psi_tmp_z_imaginary=new double [nofstates];


  /*set the return measuremt*/
  int total_steps=(int) (T/tau);

  /* prepare array for some mesurable output*/
  spin_return = new double [N*3*(total_steps+1)]();//3 state for x,y,z
  energy_sys_return = new double [total_steps+1]();
  energy_env_return = new double [total_steps+1]();
  energy_se_return  = new double [total_steps+1]();
  energy_all_return = new double [total_steps+1]();
  success_probability_return = 0;
  coefficient_return = new double [nofstates]();

  /*set array for skip_zeroterm()*/
  // count_h=0;
  // h_k_marked=new double [N_sys+N_env]();
  int count=0;
  for (int k = 0; k <N ; k++) {
    for (int l = k+1; l < N; l++) {
      if(k>=N_sys){
        if (abs(Jz_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15||abs(Jy_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15||abs(Jx_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15)
          count+=1;
      } else if (l>=N_sys && k<N_sys) {
        if (abs(Jx_se[k+(l-N_sys)*N_sys])>1e-15||abs(Jy_se[k+(l-N_sys)*N_sys])>1e-15||abs(Jz_se[k+(l-N_sys)*N_sys])>1e-15)
          count+=1;
      } else {
        if (abs(J_x[k+l*N_sys])>1e-15||abs(J_y[k+l*N_sys])>1e-15||abs(J_z[k+l*N_sys])>1e-15)
          count+=1;
      }
    }
  }

  count_z_eng=0;
  Jz_k_eng_marked=new double [count]();//[(N_sys)*((N_sys)+1)/2]();
  Jz_l_eng_marked=new double [count]();//[(N_sys)*((N_sys)+1)/2]();
  Jz_J_eng_marked=new double [count]();//[(N_sys)*((N_sys)+1)/2]();
  count_z=0;
  count_y=0;
  count_x=0;
  Jz_k_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jz_l_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jz_J_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jy_k_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jy_l_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jy_J_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jx_k_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jx_l_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jx_J_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  J_k_combine_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  J_l_combine_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jx_J_combine_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jy_J_combine_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  Jz_J_combine_marked=new double [count]();//[(N_sys+N_env)*((N_sys+N_env)+1)/2]();
  count_combine=0;
  // cout<<"count="<<count<<endl;


  Gamma=1.; //time evolution of the initial Hamiltonian. Goes from 1 to 0
  Delta=1.; //time evolution of the desired Hamiltonian. Goes from 0 to 1


}

/* Set up the initial basis state
  Input:
    d: define which mode to set up
      "read": read "output_psi_real.dat" and "psi_imagine_dat"
      "1upRan": set spin_0 as up in x and others are random.
      "allx": as the ground state of sigma_x.
      "allRand": all random.
*/
void spin_system::set_initial_sys_state(char const * d){

  if ("read"==d){
    read((int)pow(2,N_sys),psi_sys_real,"output_psi_real.dat");
    read((int)pow(2,N_sys),psi_sys_imaginary,"output_psi_imagine.dat");
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
    ofstream Psi_r_out("output_psi_real.dat");
    ofstream Psi_i_out("output_psi_imagine.dat");
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
    ofstream Psi_r_out("output_psi_real.dat");
    ofstream Psi_i_out("output_psi_imagine.dat");
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
    ofstream Psi_r_out("output_psi_real.dat");
    ofstream Psi_i_out("output_psi_imagine.dat");
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
  // for (int i = 0; i < count_h; i++) {

    double hx=0.;
    double hy=0.;
    double hz=0.;
    // int k=h_k_marked[i];
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


    if (norm<1e-14) {
      continue;
      // ss_operator_real[0]      = 1;//cos(tau*norm*0.5);
      // ss_operator_real[1]      = 0;
      // ss_operator_real[2]      = 0;
      // ss_operator_real[3]      = 1;//cos(tau*norm*0.5);
      // ss_operator_imaginary[0] = 0;
      // ss_operator_imaginary[1] = 0;
      // ss_operator_imaginary[2] = 0;
      // ss_operator_imaginary[3] = 0;
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
  // for (int k = 0; k <N ; k++) {
    // for (int l = k+1; l < N; l++) {
  for (int i = 0; i<count_x; i++) {
    int k=Jx_k_marked[i];
    int l=Jx_l_marked[i];
    double J=Jx_J_marked[i];
    if(l<N_sys){
      J=Delta*J;
    }
      // double J=0.;
      // if (k>=N_sys) {
      //   J=Jx_env[(k-N_sys)+(l-N_sys)*N_env];
      // } else if(l>=N_sys && k<N_sys) {
      //   J=Jx_se[k+(l-N_sys)*N_sys];
      // } else {
      //   J=Delta*J_x[k+l*N_sys];
      // }
  // for (int i = 0; i < count_x; i++) { //optimize use
    // int k=Jx_k_marked[i];
    // int l=Jx_l_marked[i];
      // if(abs(J)>1e-15){

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

        // ds_operator_real[0]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_real[2]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_real[5]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_real[7]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_imaginary[1] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_imaginary[3] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_imaginary[4] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_imaginary[6] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,cos_b_block,cos_c_block,sin_b_block,sin_c_block)
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

          psi_real_temp_n0      = // ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  // /*-ds_operator_imaginary[0]*psi_imaginary[n0]*/ - ds_operator_imaginary[1]*psi_imaginary[n3];
                                  cos_b_block*psi_real[n0]- sin_b_block*psi_imaginary[n3];
          psi_imaginary_temp_n0 = // /*ds_operator_imaginary[0]*psi_real[n0] +*/ ds_operator_imaginary[1]*psi_real[n3]
                                  // +ds_operator_real[0]*psi_imaginary[n0] /*+ ds_operator_real[1]*psi_imaginary[n3]*/;
                                  sin_b_block*psi_real[n3]+cos_b_block*psi_imaginary[n0];

          psi_real_temp_n1      = // ds_operator_real[2]*psi_real[n1] /*+ ds_operator_real[3]*psi_real[n2]*/
                                  // /*-ds_operator_imaginary[2]*psi_imaginary[n1]*/ - ds_operator_imaginary[3]*psi_imaginary[n2];
                                  cos_c_block*psi_real[n1]- sin_c_block*psi_imaginary[n2];
          psi_imaginary_temp_n1 = // /*ds_operator_imaginary[2]*psi_real[n1] +*/ ds_operator_imaginary[3]*psi_real[n2]
                                  // +ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;
                                  sin_c_block*psi_real[n2]+cos_c_block*psi_imaginary[n1];

          psi_real_temp_n2      = // /*ds_operator_real[4]*psi_real[n1] + */ds_operator_real[5]*psi_real[n2]
                                  // -ds_operator_imaginary[4]*psi_imaginary[n1] /*- ds_operator_imaginary[5]*psi_imaginary[n2]*/;
                                  cos_c_block*psi_real[n2]-sin_c_block*psi_imaginary[n1];
          psi_imaginary_temp_n2 = // ds_operator_imaginary[4]*psi_real[n1] /*+ ds_operator_imaginary[5]*psi_real[n2]*/
                                  // /*+ds_operator_real[4]*psi_imaginary[n1] */+ ds_operator_real[5]*psi_imaginary[n2];
                                  sin_c_block*psi_real[n1] +cos_c_block*psi_imaginary[n2];

          psi_real_temp_n3      = // /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  // -ds_operator_imaginary[6]*psi_imaginary[n0] /*- ds_operator_imaginary[7]*psi_imaginary[n3]*/;
                                  cos_b_block*psi_real[n3]-sin_b_block*psi_imaginary[n0];
          psi_imaginary_temp_n3 = // ds_operator_imaginary[6]*psi_real[n0] /*+ ds_operator_imaginary[7]*psi_real[n3]*/
                                  // /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];
                                  sin_b_block*psi_real[n0] + cos_b_block*psi_imaginary[n3];

          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      // ;}
    // ;}
  // ;}
  }
}

/* Operating sigma_y*sigma_y
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[],ds_operator_real[],ds_operator_imaginary[]
*/
void spin_system::double_spin_op_y(double t){
  // for (int k = 0; k <N ; k++) {
    // for (int l = k+1; l < N; l++) {
  for (int i = 0; i < count_y; i++) {
    int k=Jy_k_marked[i];
    int l=Jy_l_marked[i];
    double J=Jy_J_marked[i];
    if(l<N_sys){
      J=Delta*J;
    }
      // double J=0.;
      // if (k>=N_sys) {
      //   J=Jy_env[(k-N_sys)+(l-N_sys)*N_env];
      // } else if(l>=N_sys && k<N_sys) {
      //   J=Jy_se[k+(l-N_sys)*N_sys];
      // } else {
      //   J=Delta*J_y[k+l*N_sys];
      // }

  // for (int i = 0; i < count_y; i++) {
    // int k=Jy_k_marked[i];
    // int l=Jy_l_marked[i];
      // if(abs(J)>1e-15){

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



        // ds_operator_real[0]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_real[1]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_real[2]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_real[3]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_real[4]      = 0;//-sin(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_real[5]      = cos_c_block;//cos(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_real[6]      = 0;//-sin( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_real[7]      = cos_b_block;//cos( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_imaginary[0] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);
        // ds_operator_imaginary[1] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_imaginary[2] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_imaginary[3] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_imaginary[4] = sin_c_block;//cos(-a*Delta*tau*0.5)*sin(c*Delta*tau*0.5);
        // ds_operator_imaginary[5] = 0;//sin(-a*Delta*tau*0.5)*cos(c*Delta*tau*0.5);
        // ds_operator_imaginary[6] = sin_b_block;//cos( a*Delta*tau*0.5)*sin(b*Delta*tau*0.5);
        // ds_operator_imaginary[7] = 0;//sin( a*Delta*tau*0.5)*cos(b*Delta*tau*0.5);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,cos_b_block,cos_c_block,sin_b_block,sin_c_block)
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

          psi_real_temp_n0      = // ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  // /*-ds_operator_imaginary[0]*psi_imaginary[n0]*/ - ds_operator_imaginary[1]*psi_imaginary[n3];
                                  cos_b_block*psi_real[n0]- sin_b_block*psi_imaginary[n3];
          psi_imaginary_temp_n0 = // /*ds_operator_imaginary[0]*psi_real[n0] +*/ ds_operator_imaginary[1]*psi_real[n3]
                                  // +ds_operator_real[0]*psi_imaginary[n0] /*+ ds_operator_real[1]*psi_imaginary[n3]*/;
                                  sin_b_block*psi_real[n3]+cos_b_block*psi_imaginary[n0];
          psi_real_temp_n1      = // ds_operator_real[2]*psi_real[n1] /*+ ds_operator_real[3]*psi_real[n2]*/
                                  // /*-ds_operator_imaginary[2]*psi_imaginary[n1]*/ - ds_operator_imaginary[3]*psi_imaginary[n2];
                                  cos_c_block*psi_real[n1]- sin_c_block*psi_imaginary[n2];
          psi_imaginary_temp_n1 = // /*ds_operator_imaginary[2]*psi_real[n1] +*/ ds_operator_imaginary[3]*psi_real[n2]
                                  // +ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;
                                  sin_c_block*psi_real[n2] +cos_c_block*psi_imaginary[n1];
          psi_real_temp_n2      = // /*ds_operator_real[4]*psi_real[n1] + */ds_operator_real[5]*psi_real[n2]
                                  // -ds_operator_imaginary[4]*psi_imaginary[n1] /*- ds_operator_imaginary[5]*psi_imaginary[n2]*/;
                                  cos_c_block*psi_real[n2] -sin_c_block*psi_imaginary[n1];
          psi_imaginary_temp_n2 = // ds_operator_imaginary[4]*psi_real[n1] /*+ ds_operator_imaginary[5]*psi_real[n2]*/
                                  // /*+ds_operator_real[4]*psi_imaginary[n1] */+ ds_operator_real[5]*psi_imaginary[n2];
                                  sin_c_block*psi_real[n1]+ cos_c_block*psi_imaginary[n2];
          psi_real_temp_n3      = // /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  // -ds_operator_imaginary[6]*psi_imaginary[n0] /*- ds_operator_imaginary[7]*psi_imaginary[n3]*/;
                                  cos_b_block*psi_real[n3] -sin_b_block*psi_imaginary[n0];
          psi_imaginary_temp_n3 = // ds_operator_imaginary[6]*psi_real[n0] /*+ ds_operator_imaginary[7]*psi_real[n3]*/
                                  // /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];
                                  sin_b_block*psi_real[n0] + cos_b_block*psi_imaginary[n3];
          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      // ;}
    // ;}
  // ;}
  }
}

/* Operating sigma_z*sigma_z
  Input:
    t: current annealing time point
  Side effect/Changed:
    psi_real[], psi_imaginary[],ds_operator_real[],ds_operator_imaginary[]
*/
void spin_system::double_spin_op_z(double t){

  // for (int k = 0; k <N ; k++) {
    // for (int l = k+1; l < N; l++) {
  for (int i = 0; i < count_z; i++) {
    int k=Jz_k_marked[i];
    int l=Jz_l_marked[i];
    double J=Jz_J_marked[i];
    if(l<N_sys){
      J=Delta*J;
    }
      //
      // double J=0.;
      // if (k>=N_sys) {
      //   J=Jz_env[(k-N_sys)+(l-N_sys)*N_env];
      // } else if(l>=N_sys && k<N_sys) {
      //   J=Jz_se[k+(l-N_sys)*N_sys];
      // } else {
      //   J=Delta*J_z[k+l*N_sys];
      // }
  // for (int i = 0; i < count_z; i++) {
    // int k=Jz_k_marked[i];
    // int l=Jz_l_marked[i];
      // if(abs(J)>1e-15){


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
        // ds_operator_real[0]      = cos_a_block;//cos( a*Delta*tau)*cos(b*Delta*tau);
        // ds_operator_real[1]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        // ds_operator_real[2]      = cos_a_block;//cos(-a*Delta*tau)*cos(c*Delta*tau);
        // ds_operator_real[3]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        // ds_operator_real[4]      =0;//-sin(-a*Delta*tau)*sin(c*Delta*tau);
        // ds_operator_real[5]      = cos_a_block;//cos(-a*Delta*tau)*cos(c*Delta*tau);
        // ds_operator_real[6]      =0;//-sin( a*Delta*tau)*sin(b*Delta*tau);
        // ds_operator_real[7]      = cos_a_block;//cos( a*Delta*tau)*cos(b*Delta*tau);
        // ds_operator_imaginary[0] = sin_a_block;//sin( a*Delta*tau)*cos(b*Delta*tau);
        // ds_operator_imaginary[1] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        // ds_operator_imaginary[2] = -sin_a_block;//sin(-a*Delta*tau)*cos(c*Delta*tau);
        // ds_operator_imaginary[3] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        // ds_operator_imaginary[4] = 0;//cos(-a*Delta*tau)*sin(c*Delta*tau);
        // ds_operator_imaginary[5] = -sin_a_block;//sin(-a*Delta*tau)*cos(c*Delta*tau);
        // ds_operator_imaginary[6] = 0;//cos( a*Delta*tau)*sin(b*Delta*tau);
        // ds_operator_imaginary[7] = sin_a_block;//sin( a*Delta*tau)*cos(b*Delta*tau);

        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,cos_a_block,sin_a_block)
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

          psi_real_temp_n0      = //ds_operator_real[0]*psi_real[n0] /*+ ds_operator_real[1]*psi_real[n3]*/
                                  //-ds_operator_imaginary[0]*psi_imaginary[n0] /*- ds_operator_imaginary[1]*psi_imaginary[n3]*/;
                                  cos_a_block*psi_real[n0] -sin_a_block*psi_imaginary[n0];
          psi_imaginary_temp_n0 = //ds_operator_imaginary[0]*psi_real[n0] /*+ ds_operator_imaginary[1]*psi_real[n3]*/
                                  //+ds_operator_real[0]*psi_imaginary[n0]/* + ds_operator_real[1]*psi_imaginary[n3]*/;
                                  sin_a_block*psi_real[n0] +cos_a_block*psi_imaginary[n0];

          psi_real_temp_n1      = //ds_operator_real[2]*psi_real[n1]           /*+ ds_operator_real[3]*psi_real[n2]*/
                                  //-ds_operator_imaginary[2]*psi_imaginary[n1] /*- ds_operator_imaginary[3]*psi_imaginary[n2]*/;
                                  cos_a_block*psi_real[n1] -(-sin_a_block)*psi_imaginary[n1];
          psi_imaginary_temp_n1 = //ds_operator_imaginary[2]*psi_real[n1] /*+ ds_operator_imaginary[3]*psi_real[n2]*/
                                  //+ds_operator_real[2]*psi_imaginary[n1] /*+ ds_operator_real[3]*psi_imaginary[n2]*/;
                                  (-sin_a_block)*psi_real[n1] +cos_a_block*psi_imaginary[n1];
          psi_real_temp_n2      = // /*ds_operator_real[4]*psi_real[n1] +*/ ds_operator_real[5]*psi_real[n2]
                                  // /*-ds_operator_imaginary[4]*psi_imaginary[n1]*/ - ds_operator_imaginary[5]*psi_imaginary[n2];
                                  cos_a_block*psi_real[n2]- (-sin_a_block)*psi_imaginary[n2];
          psi_imaginary_temp_n2 = // /*ds_operator_imaginary[4]*psi_real[n1] */+ ds_operator_imaginary[5]*psi_real[n2]
                                  // /*+ds_operator_real[4]*psi_imaginary[n1]*/ + ds_operator_real[5]*psi_imaginary[n2];
                                  (-sin_a_block)*psi_real[n2]+ cos_a_block*psi_imaginary[n2];
          psi_real_temp_n3      = // /*ds_operator_real[6]*psi_real[n0] +*/ ds_operator_real[7]*psi_real[n3]
                                  // /*-ds_operator_imaginary[6]*psi_imaginary[n0]*/ - ds_operator_imaginary[7]*psi_imaginary[n3];
                                  cos_a_block*psi_real[n3]- sin_a_block*psi_imaginary[n3];
          psi_imaginary_temp_n3 = // /*ds_operator_imaginary[6]*psi_real[n0] */+ ds_operator_imaginary[7]*psi_real[n3]
                                  // /*+ds_operator_real[6]*psi_imaginary[n0]*/ + ds_operator_real[7]*psi_imaginary[n3];
                                  sin_a_block*psi_real[n3]+ cos_a_block*psi_imaginary[n3];
          psi_real[n0]      = psi_real_temp_n0;
          psi_imaginary[n0] = psi_imaginary_temp_n0;
          psi_real[n1]      = psi_real_temp_n1;
          psi_imaginary[n1] = psi_imaginary_temp_n1;
          psi_real[n2]      = psi_real_temp_n2;
          psi_imaginary[n2] = psi_imaginary_temp_n2;
          psi_real[n3]      = psi_real_temp_n3;
          psi_imaginary[n3] = psi_imaginary_temp_n3;
        }
      // ;}
    // ;}
  // ;}
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
  #: since now the system doesn't contain Jy,Jx, i skip this. If needed should uncomment the complete code
*/
double spin_system::energy_sys(double t){
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
    double Gammahx=Gamma*hx;
    double Deltahz=Delta*hz;
    #pragma omp parallel for default(none) shared(i1,hx,hz,Gammahx,Deltahz)
    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l -i2+i2/i1;
      int j = i+i1;
      // /*sigma_x*/
      // psi_tmp_real[i]     += Gamma*hx*psi_real[j];
      // psi_tmp_imaginary[i]+= Gamma*hx*psi_imaginary[j];
      // psi_tmp_real[j]     += Gamma*hx*psi_real[i];
      // psi_tmp_imaginary[j]+= Gamma*hx*psi_imaginary[i];
      //
      // /*sigma_z*/
      // psi_tmp_real[i]     += Delta*hz*psi_real[i];
      // psi_tmp_imaginary[i]+= Delta*hz*psi_imaginary[i];
      // psi_tmp_real[j]     += -Delta*hz*psi_real[j];
      // psi_tmp_imaginary[j]+= -Delta*hz*psi_imaginary[j];

      //collapse the same term to one equation, above are the reference.
      psi_tmp_real[i]     += Gammahx*psi_real[j] +Deltahz*psi_real[i];
      psi_tmp_imaginary[i]+= Gammahx*psi_imaginary[j] +Deltahz*psi_imaginary[i];
      psi_tmp_real[j]     += Gammahx*psi_real[i] -Deltahz*psi_real[j];
      psi_tmp_imaginary[j]+= Gammahx*psi_imaginary[i] -Deltahz*psi_imaginary[j];
    }
  }

  // for (int k = 0; k <N_sys ; k++) {
  //   for (int l = k+1; l < N_sys; l++) {
  //     double Jx=0;//-1*J_x[k+l*N_sys];
  //     double Jy=0;//-1*J_y[k+l*N_sys];
  //     double Jz=-1*J_z[k+l*N_sys];
  //     if(abs(Jz)>1e-15||abs(Jy)>1e-15||abs(Jx)>1e-15){
  for (int i = 0; i < count_z_eng; i++) {
    int k=Jz_k_eng_marked[i];
    int l=Jz_l_eng_marked[i];
    // if (k>=N_sys||l>=N_sys)
    //    continue;
      double Jx=0;
      double Jy=0;
      double Jz=-1*Jz_J_eng_marked[i];
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        double DelJz=Delta*Jz;
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz,DelJz)
        for (int m = 0; m < nofstates; m+=4) {
          int n3 = m & njj;
          int n2 = m-n3+(n3+n3)/njj;
          int n1 = n2 & nii;
          int n0 = n2-n1+n1/nii;
          n1=n0+nii;
          n2=n0+njj;
          n3=n1+njj;

          // /*sigma_x*sigma_x*/
          // psi_tmp_real[n0]      += Delta*Jx*psi_real[n3];
          // psi_tmp_imaginary[n0] += Delta*Jx*psi_imaginary[n3];
          // psi_tmp_real[n1]      += Delta*Jx*psi_real[n2];
          // psi_tmp_imaginary[n1] += Delta*Jx*psi_imaginary[n2];
          // psi_tmp_real[n2]      += Delta*Jx*psi_real[n1];
          // psi_tmp_imaginary[n2] += Delta*Jx*psi_imaginary[n1];
          // psi_tmp_real[n3]      += Delta*Jx*psi_real[n0];
          // psi_tmp_imaginary[n3] += Delta*Jx*psi_imaginary[n0];
          // /*sigma_y*sigma_y*/
          // psi_tmp_real[n0]      += -Delta*Jy*psi_real[n3];
          // psi_tmp_imaginary[n0] += -Delta*Jy*psi_imaginary[n3];
          // psi_tmp_real[n1]      += Delta*Jy*psi_real[n2];
          // psi_tmp_imaginary[n1] += Delta*Jy*psi_imaginary[n2];
          // psi_tmp_real[n2]      += Delta*Jy*psi_real[n1];
          // psi_tmp_imaginary[n2] += Delta*Jy*psi_imaginary[n1];
          // psi_tmp_real[n3]      += -Delta*Jy*psi_real[n0];
          // psi_tmp_imaginary[n3] += -Delta*Jy*psi_imaginary[n0];
          // /*sigma_z*sigma_z*/
          psi_tmp_real[n0]      += DelJz*psi_real[n0];//Delta*Jz*psi_real[n0];
          psi_tmp_imaginary[n0] += DelJz*psi_imaginary[n0];//Delta*Jz*psi_imaginary[n0];
          psi_tmp_real[n1]      += -DelJz*psi_real[n1];//-Delta*Jz*psi_real[n1];
          psi_tmp_imaginary[n1] += -DelJz*psi_imaginary[n1];//-Delta*Jz*psi_imaginary[n1];
          psi_tmp_real[n2]      += -DelJz*psi_real[n2];//-Delta*Jz*psi_real[n2];
          psi_tmp_imaginary[n2] += -DelJz*psi_imaginary[n2];//-Delta*Jz*psi_imaginary[n2];
          psi_tmp_real[n3]      += DelJz*psi_real[n3];//Delta*Jz*psi_real[n3];
          psi_tmp_imaginary[n3] += DelJz*psi_imaginary[n3];//Delta*Jz*psi_imaginary[n3];
          //collapse the same term to one equation, above are the reference.
          // psi_tmp_real[n0]      += Delta* /*((Jx-Jy)*psi_real[n3]*/+Jz*psi_real[n0];
          // psi_tmp_imaginary[n0] += Delta* /*((Jx-Jy)*psi_imaginary[n3]*/+Jz*psi_imaginary[n0];
          // psi_tmp_real[n1]      += Delta* /*((Jx+Jy)*psi_real[n2]*/-Jz*psi_real[n1];
          // psi_tmp_imaginary[n1] += Delta* /*((Jx+Jy)*psi_imaginary[n2]*/-Jz*psi_imaginary[n1];
          // psi_tmp_real[n2]      += Delta* /*((Jx+Jy)*psi_real[n1]*/-Jz*psi_real[n2];
          // psi_tmp_imaginary[n2] += Delta* /*((Jx+Jy)*psi_imaginary[n1]*/-Jz*psi_imaginary[n2];
          // psi_tmp_real[n3]      += Delta* /*((Jx-Jy)*psi_real[n0]*/+Jz*psi_real[n3];
          // psi_tmp_imaginary[n3] += Delta* /*((Jx-Jy)*psi_imaginary[n0]*/+Jz*psi_imaginary[n3];


        }
      // }
    // }
  }
  for (int i = 0; i < nofstates; ++i) {
    average_energy += psi_real[i]*psi_tmp_real[i] - -1*psi_imaginary[i]*psi_tmp_imaginary[i];
    // check_img += psi_real[i]*psi_tmp_imaginary[i] + -1*psi_imaginary[i]*psi_tmp_real[i];
  }
  // if (abs(check_img)>1e-13)
  //   cout<<"Something went wrong in functoin energy()   "<<check_img<<endl;

  return average_energy;
}


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
  complex<double> ImgNum(0,1);
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

    // cout<<"hx hy hz ="<<hx<<" "<<hy<<" "<<hz<<" "<<endl;

    for (int l = 0; l < nofstates; l+=2) {
      int i2= l & i1;
      int i = l - i2 +i2/i1;
      int j = i+i1;
        H_env[i+i*(i+1)/2] += hz;
        H_env[j+j*(j+1)/2] += -hz;
        H_env[i+j*(j+1)/2] += -hy*ImgNum;
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

  int lda = nofstates;
  complex<double> work[2*nofstates];
  double rwork[7*nofstates];
  int iwork[5*nofstates];
  int ifail[nofstates];
  int info;

  zhpevx("Vector","All","Upper", &n, H_env, &vl, &vu, &il, &iu, &abstol, &m, w, z, &lda, work, rwork, iwork, ifail, &info );
  if(info!=0){
    cout<<"info = "<<info<<endl;
  }
  ofstream eng_out("output_env_energy_spectrum.dat");
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

  ofstream envr_out("output_env_partition_factor.dat");
  ofstream env_basis_out("output_env_basis.dat");
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

// /* !!!!!!!!WILL BE REMOVE in the future, since it might be a wrong implementation!!!!!!!
//   Try to read the initial basis state from the system and the Enivironment
//   And then combine then together into a new state.
//   Input:
//     int: the # of total spin.
//     double*: the real part
//     double*: the imaginary part
//     char const*: the basis state of the system
//     char const*: the basis state of the environment
// */
// void spin_system::generate(int N, double* array_real, double* array_imagine, char const* filename_sysr, char const* filename_sysi, char const* filename_envr, char const* filename_envi ){
//
//   double *sys_real, *sys_imag, *env_real, *env_imag;
//   // int N_half=(int) pow(2,N/2);
//   int nos_sys=(int) pow(2,N_sys);
//   int nos_env=(int) pow(2,N_env);
//   sys_real=new double[nos_sys];
//   sys_imag=new double[nos_sys];
//   env_real=new double[nos_env];
//   env_imag=new double[nos_env];
//   read(nos_sys, sys_real, filename_sysr);
//   read(nos_sys, sys_imag, filename_sysi);
//   read(nos_env, env_real, filename_envr);
//   read(nos_env, env_imag, filename_envi);
//
//   ofstream state_out("output_state_complete.dat");
//   for (int i = 0; i < nos_env; i++) {
//     for (int j = 0; j < nos_sys; j++) {
//       array_real[i*nos_sys+j]=env_real[i]*sys_real[j]-env_imag[i]*sys_imag[j];
//       array_imagine[i*nos_sys+j]=env_imag[i]*sys_real[j]+env_real[i]*sys_imag[j];
//       state_out<<array_real[i*nos_sys+j]<<" "<<array_imagine[i*nos_sys+j]<<endl;
//     }
//   }
//   double sum=0;
//   int nos=(int) pow(2,N_sys+N_env);
//   for (int i = 0; i < nos; i++) {
//     sum+=array_real[i]*array_real[i]+array_imagine[i]*array_imagine[i];
//   }
//   // cout<<"sum inside generate functoin: "<<sum<<endl;
// }

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

/* Skin zero term and put non-zero term into a 1D array.
  Input:
    none
  Side Effect:
  count_z, count_y, count_x
  Jz_k_marked, Jz_l_marked, Jz_J_marked, Jy_k_marked, Jy_l_marked, Jy_J_marked, Jx_k_marked, Jx_l_marked, Jx_J_marked,
*/
void spin_system::skip_zeroterm(){
  /*skip zero term in h*/
  // for (int k = 0; k <(N_sys+N_env) ; k++) {
  //   if(k>=N_sys){
  //     h_k_marked[count_h]=k;
  //     count_h+=1;
  //   } else {
  //     if((h_x[k]*h_x[k]+h_y[k]*h_y[k]+h_z[k]*h_z[k])>1e-15)
  //       h_k_marked[count_h]=k;
  //       count_h+=1;
  //   }
  // }
  /*skip zero term in J*/
  for (int k = 0; k <(N_sys+N_env) ; k++) {
    for (int l = k+1; l < (N_sys+N_env); l++) {
      if(k>=N_sys){
        if (abs(Jx_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15){
          Jx_k_marked[count_x]=k;
          Jx_l_marked[count_x]=l;
          Jx_J_marked[count_x]=Jx_env[(k-N_sys)+(l-N_sys)*N_env];
          count_x+=1;
        }
        if (abs(Jy_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15){
          Jy_k_marked[count_y]=k;
          Jy_l_marked[count_y]=l;
          Jy_J_marked[count_y]=Jy_env[(k-N_sys)+(l-N_sys)*N_env];
          count_y+=1;
        }
        if (abs(Jz_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15){
          Jz_k_marked[count_z]=k;
          Jz_l_marked[count_z]=l;
          Jz_J_marked[count_z]=Jz_env[(k-N_sys)+(l-N_sys)*N_env];
          count_z+=1;
        }
        if (abs(Jz_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15&&abs(Jy_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15&&abs(Jx_env[(k-N_sys)+(l-N_sys)*N_env])>1e-15) {
          J_k_combine_marked[count_combine]=k;
          J_l_combine_marked[count_combine]=l;
          Jx_J_combine_marked[count_combine]=Jx_env[(k-N_sys)+(l-N_sys)*N_env];
          Jy_J_combine_marked[count_combine]=Jy_env[(k-N_sys)+(l-N_sys)*N_env];
          Jz_J_combine_marked[count_combine]=Jz_env[(k-N_sys)+(l-N_sys)*N_env];
          count_combine+=1;
        }

      } else if (l>=N_sys && k<N_sys) {
        if (abs(Jx_se[k+(l-N_sys)*N_sys])>1e-15){
          Jx_k_marked[count_x]=k;
          Jx_l_marked[count_x]=l;
          Jx_J_marked[count_x]=Jx_se[k+(l-N_sys)*N_sys];
          count_x+=1;
        }
        if (abs(Jy_se[k+(l-N_sys)*N_sys])>1e-15){
          Jy_k_marked[count_y]=k;
          Jy_l_marked[count_y]=l;
          Jy_J_marked[count_y]=Jy_se[k+(l-N_sys)*N_sys];
          count_y+=1;
        }
        if (abs(Jz_se[k+(l-N_sys)*N_sys])>1e-15){
          Jz_k_marked[count_z]=k;
          Jz_l_marked[count_z]=l;
          Jz_J_marked[count_z]=Jz_se[k+(l-N_sys)*N_sys];
          count_z+=1;
        }
        if (abs(Jx_se[k+(l-N_sys)*N_sys])>1e-15&&abs(Jy_se[k+(l-N_sys)*N_sys])>1e-15&&abs(Jz_se[k+(l-N_sys)*N_sys])>1e-15) {
          J_k_combine_marked[count_combine]=k;
          J_l_combine_marked[count_combine]=l;
          Jx_J_combine_marked[count_combine]=Jx_se[k+(l-N_sys)*N_sys];
          Jy_J_combine_marked[count_combine]=Jy_se[k+(l-N_sys)*N_sys];
          Jz_J_combine_marked[count_combine]=Jz_se[k+(l-N_sys)*N_sys];
          count_combine+=1;
        }

      } else {
        if (abs(J_x[k+l*N_sys])>1e-15){
          Jx_k_marked[count_x]=k;
          Jx_l_marked[count_x]=l;
          Jx_J_marked[count_x]=J_x[k+l*N_sys];
          count_x+=1;
        }
        if (abs(J_y[k+l*N_sys])>1e-15){
          Jy_k_marked[count_y]=k;
          Jy_l_marked[count_y]=l;
          Jy_J_marked[count_y]=J_y[k+l*N_sys];
          count_y+=1;
        }
        if (abs(J_z[k+l*N_sys])>1e-15){
          Jz_k_marked[count_z]=k;
          Jz_l_marked[count_z]=l;
          Jz_J_marked[count_z]=J_z[k+l*N_sys];
          count_z+=1;
          Jz_k_eng_marked[count_z_eng]=k;
          Jz_l_eng_marked[count_z_eng]=l;
          Jz_J_eng_marked[count_z_eng]=J_z[k+l*N_sys];
          count_z_eng+=1;
        }
        if (abs(J_x[k+l*N_sys])>1e-15&&abs(J_y[k+l*N_sys])>1e-15&&abs(J_z[k+l*N_sys])>1e-15) {
          J_k_combine_marked[count_combine]=k;
          J_l_combine_marked[count_combine]=l;
          Jx_J_combine_marked[count_combine]=J_x[k+l*N_sys];
          Jy_J_combine_marked[count_combine]=J_y[k+l*N_sys];
          Jz_J_combine_marked[count_combine]=J_z[k+l*N_sys];
          count_combine+=1;
        }
      }
    }
  }
 cout<<"# of elements in Jx, Jy, and Jz: "<<count_x<<" "<<count_y<<" "<<count_z<<" "<<endl;

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

  double norm=0.;
  for (int E_i = 0; E_i < (int) pow(2,N_env); E_i++) {
    direct_product(E_i,psi_real,psi_imaginary,z,psi_sys_real,psi_sys_imaginary);
    for (int i = 0; i < nofstates; i++)
      coefficient_return[i]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);
  }
  for (int i = 0; i < nofstates; i++)
    norm+=coefficient_return[i];
  cout<<"norm before run: "<<norm<<endl;
  for (int i = 0; i < nofstates; i++)
    coefficient_return[i]=0.;

  cout<<"START RUNNING!"<<endl;
  for (int E_i = 0; E_i < (int) pow(2,N_env); E_i++) {
    if (abs(w[E_i]-0)<1e-8)
      continue;
    direct_product(E_i,psi_real,psi_imaginary,z,psi_sys_real,psi_sys_imaginary);

    for (int step = 0; step < total_steps+1; step++){ //+1 because i count the 0 point and the last poing as well.

      Delta=step*tau/T;
      Gamma=1-Delta;

      if (step%2500==0)
        cout<<"E_i= "<<E_i<<", w[]= "<<w[E_i]<<", step: "<<step<<endl;

      energy_sys_return[step]+=w[E_i]*energy_sys(step*tau);
      energy_env_return[step]+=w[E_i]*energy_env(step*tau);
      // energy_se_return[step]+=w[E_i]*energy_se(step*tau);
      // energy_all_return[step]+=w[E_i]*energy_all(step*tau);
      // for (int s = 0; s < N; s++) {
      //   int index=step*N*3+s*3;
      //   spin_return[index]  +=w[E_i]*spin('x',s);
      //   spin_return[index+1]+=w[E_i]*spin('y',s);
      //   spin_return[index+2]+=w[E_i]*spin('z',s);
      // }

      for (int i = gs_sol; i < nofstates; i+=256)
        frequency[step]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);

      // frequency[step]+=w[E_i]*(psi_real[gs_sol]*psi_real[gs_sol]+psi_imaginary[gs_sol]*psi_imaginary[gs_sol]);
      single_spin_op(step*tau);
      double_spin_op_x(step*tau);
      double_spin_op_y(step*tau);
      double_spin_op_z(step*tau);
      double_spin_op_y(step*tau);
      double_spin_op_x(step*tau);
      single_spin_op(step*tau);

    }
    for (int i = 0; i < nofstates; i++)
      coefficient_return[i]+=w[E_i]*(psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i]);

  }
  cout<<"END RUNNING"<<endl;
  norm=0.;
  for (int i = 0; i < nofstates; i++)
    norm+=coefficient_return[i];

  cout<<"norm after run: "<<norm<<endl;
  // output the return value: coefficient, energy expectation value, and spin expectation value.
  ofstream Coefficient_out("output_coefficient.dat");
  for (int i = 0; i < nofstates; i++)
    Coefficient_out<<coefficient_return[i]<<endl;

  ostringstream strs;
  strs <<"G"<<(G/10.)<<"_"<<"Ts"<<(T*10)<<".dat";//output name for single output
  //strs <<"H"<<J_index<<"_"<<"Ts"<<(T*10)<<".dat";//output name for multiple output for landau ziener comparison
  string str = strs.str();
  //string strmain="/home/zam/t.hsu/Documents/2016_WS_MT/Product_formula/gap_result/output_general_";
  //strmain.append(str);
  const char *testChars = str.c_str();
  ofstream output(testChars);
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
  success_probability_return=frequency[total_steps];
}



/* Produce a random state with pow(2,N_env)
  Input:
    N_env: the # of environment spins
  Side Effect:
    RANDOM_FACTOR[]
*/
void spin_system::set_random(double N_env){
  int nos_env=(int) pow(2,N_env);
  RANDOM_FACTOR=new complex<double> [nos_env]();
  complex<double> ImgNum(0,1);
  double rand_tmp=0;
  srand(time(NULL));
  double mean=0.0;
  double std=1./3;
  double u,v;
  double Norm=0;
  for (int i = 0; i < nos_env; i++) {
    u = rand() / (double)RAND_MAX;
    v = rand() / (double)RAND_MAX;
    rand_tmp= sqrt(-2*log(u))*cos(2*M_PI*v)*std+mean;
    RANDOM_FACTOR[i]=rand_tmp;

    u = rand() / (double)RAND_MAX;
    v = rand() / (double)RAND_MAX;
    rand_tmp = sqrt(-2*log(u))*sin(2*M_PI*v)*std+mean;
    RANDOM_FACTOR[i]+=rand_tmp*ImgNum;
    Norm+=norm(RANDOM_FACTOR[i]);
  }

  Norm=1./sqrt(Norm);
  for (int i = 0; i < nos_env; i++)
    RANDOM_FACTOR[i]=RANDOM_FACTOR[i]*Norm;
}


/* direct product a random factor with the sub_sys wavefunction
  Intput:
    array_real: prepare the output real part
    array_imagine: prepare the output imaginary pary
    sys_real: the input sys real
    sys_imag: the input sys imag
  Side effecht/Changed:
    psi_real[],psi_imaginary[]
  *also take RANDOM_FACTOR as input, so please call set_random() before use.
*/
void spin_system::random_product(double* array_real, double* array_imagine, double* sys_real, double* sys_imag){
  int nos_sys=(int) pow(2,N_sys);
  int nos_env=(int) pow(2,N_env);

  for (int i = 0; i < nos_env; i++) {
    for (int j = 0; j < nos_sys; j++) {
      array_real[i*nos_sys+j]=RANDOM_FACTOR[i].real()*sys_real[j]-RANDOM_FACTOR[i].imag()*sys_imag[j];
      array_imagine[i*nos_sys+j]=RANDOM_FACTOR[i].imag()*sys_real[j]+RANDOM_FACTOR[i].real()*sys_imag[j];
    }
  }
}

/* operation to approach e^-BH/2m with [1+(-BH/2m)]
  Intput:
    M: how precise we want to approximate the e^-BH/2m
  Side effecht/Changed:
    psi_real[],psi_imaginary[]
*/
void spin_system::exp_appr_op(int M){
  for (int i = 0; i < nofstates; i++) {
    psi_tmp_real[i]=0;
    psi_tmp_imaginary[i]=0;
  }

  for (int k = N-N_env; k <N ; k++) {
    double hx=-1*h_x[k];
    double hy=-1*h_y[k];
    double hz=-1*h_z[k];

    if(abs(hx)>1e-15||abs(hy)>1e-15||abs(hz)>1e-15){
      int i1=(int) pow(2,k);
      #pragma omp parallel for default(none) shared(i1,hx,hy,hz)
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

        // double psi_real_temp_i      = 0;
        // double psi_imaginary_temp_i = 0;
        // double psi_real_temp_j      = 0;
        // double psi_imaginary_temp_j = 0;
        // // combine to one
        // psi_real_temp_i      = hx*psi_tmp_real[j]+hy*psi_tmp_imaginary[j]+hz*psi_tmp_real[i];
        // psi_imaginary_temp_i = hx*psi_tmp_imaginary[j]-hy*psi_tmp_real[j]+hz*psi_tmp_imaginary[i];
        // psi_real_temp_j      = hx*psi_tmp_real[i]-hy*psi_tmp_imaginary[i]-hz*psi_tmp_real[j];
        // psi_imaginary_temp_j = hx*psi_tmp_imaginary[i]+hy*psi_tmp_real[i]-hz*psi_tmp_imaginary[j];
        // psi_tmp_real[i]      += psi_real_temp_i;
        // psi_tmp_imaginary[i] += psi_imaginary_temp_i;
        // psi_tmp_real[j]      += psi_real_temp_j;
        // psi_tmp_imaginary[j] += psi_imaginary_temp_j;

        //added 03.04.2017
        // psi_tmp_real[i]      += hx*psi_real[j]+hy*psi_imaginary[j]+hz*psi_real[i];
        // psi_tmp_imaginary[i] += hx*psi_imaginary[j]-hy*psi_real[j]+hz*psi_imaginary[i];
        // psi_tmp_real[j]      += hx*psi_real[i]-hy*psi_imaginary[i]-hz*psi_real[j];
        // psi_tmp_imaginary[j] += hx*psi_imaginary[i]+hy*psi_real[i]-hz*psi_imaginary[j];

      }
    }
  }
  // for (int k = 0; k <N ; k++) {
  //   for (int l = k+1; l < N; l++) {
  //     double Jx=0.;
  //     double Jy=0.;
  //     double Jz=0.;
  //     if(k>=N_sys){
  //       Jx=-1*Jx_env[(k-N_sys)+(l-N_sys)*N_env];
  //       Jy=-1*Jy_env[(k-N_sys)+(l-N_sys)*N_env];
  //       Jz=-1*Jz_env[(k-N_sys)+(l-N_sys)*N_env];
  //     }
  //     // else if(l>=N_sys && k<N_sys){
  //     //   Jx=-1*Jx_se[k+(l-N_sys)*N_sys];
  //     //   Jy=-1*Jy_se[k+(l-N_sys)*N_sys];
  //     //   Jz=-1*Jz_se[k+(l-N_sys)*N_sys];
  //     // } else {
  //     //   Jx=-1*J_x[k+l*N_sys]*Delta;
  //     //   Jy=-1*J_y[k+l*N_sys]*Delta;
  //     //   Jz=-1*J_z[k+l*N_sys]*Delta;
  //     // }





  // int test=0;
  // for (int i = 0; i<count_combine; i++) {
  //   int k=J_k_combine_marked[i];
  //   int l=J_l_combine_marked[i];
  //   double Jx=-1*Jx_J_combine_marked[i];
  //   double Jy=-1*Jy_J_combine_marked[i];
  //   double Jz=-1*Jz_J_combine_marked[i];
  //
  //   if(k<N_sys){
  //     continue;
  //     // Jx=0;//Delta*Jx;
  //     // Jy=0;//Delta*Jy;
  //     // Jz=0;//Delta*Jz;
  //   }
    // cout<<Jx<<" "<<Jy<<" "<<Jz<<endl;
    for (int k = N-N_env; k <N ; k++) {
      for (int l = k+1; l < N; l++) {
        double Jx=-1*Jx_env[(k-N_sys)+(l-N_sys)*N_env];
        double Jy=-1*Jy_env[(k-N_sys)+(l-N_sys)*N_env];
        double Jz=-1*Jz_env[(k-N_sys)+(l-N_sys)*N_env];

      if(abs(Jx)>1e-15||abs(Jy)>1e-15||abs(Jz)>1e-15){
        int nii=(int) pow(2,k);
        int njj=(int) pow(2,l);
        #pragma omp parallel for default(none) shared(nii,njj,Jx,Jy,Jz,k)
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

          // double psi_real_temp_n0      = 0;
          // double psi_imaginary_temp_n0 = 0;
          // double psi_real_temp_n1      = 0;
          // double psi_imaginary_temp_n1 = 0;
          // double psi_real_temp_n2      = 0;
          // double psi_imaginary_temp_n2 = 0;
          // double psi_real_temp_n3      = 0;
          // double psi_imaginary_temp_n3 = 0;
          //
          // psi_real_temp_n0      =Jx*psi_tmp_real[n3]-Jy*psi_tmp_real[n3]+Jz*psi_tmp_real[n0];
          // psi_imaginary_temp_n0 =Jx*psi_tmp_imaginary[n3]-Jy*psi_tmp_imaginary[n3]+Jz*psi_tmp_imaginary[n0];
          // psi_real_temp_n1      =Jx*psi_tmp_real[n2]+Jy*psi_tmp_real[n2]-Jz*psi_tmp_real[n1];
          // psi_imaginary_temp_n1 =Jx*psi_tmp_imaginary[n2]+Jy*psi_tmp_imaginary[n2]-Jz*psi_tmp_imaginary[n1];
          // psi_real_temp_n2      =Jx*psi_tmp_real[n1]+Jy*psi_tmp_real[n1]-Jz*psi_tmp_real[n2];
          // psi_imaginary_temp_n2 =Jx*psi_tmp_imaginary[n1]+Jy*psi_tmp_imaginary[n1]-Jz*psi_tmp_imaginary[n2];
          // psi_real_temp_n3      =Jx*psi_tmp_real[n0]-Jy*psi_tmp_real[n0]+Jz*psi_tmp_real[n3];
          // psi_imaginary_temp_n3 =Jx*psi_tmp_imaginary[n0]-Jy*psi_tmp_imaginary[n0]+Jz*psi_tmp_imaginary[n3];
          //
          // psi_tmp_real[n0]      += psi_real_temp_n0;
          // psi_tmp_imaginary[n0] += psi_imaginary_temp_n0;
          // psi_tmp_real[n1]      += psi_real_temp_n1;
          // psi_tmp_imaginary[n1] += psi_imaginary_temp_n1;
          // psi_tmp_real[n2]      += psi_real_temp_n2;
          // psi_tmp_imaginary[n2] += psi_imaginary_temp_n2;
          // psi_tmp_real[n3]      += psi_real_temp_n3;
          // psi_tmp_imaginary[n3] += psi_imaginary_temp_n3;

          // //added 03.04.2017
          // psi_tmp_real[n0]      +=Jx*psi_real[n3]-Jy*psi_real[n3]+Jz*psi_real[n0];
          // psi_tmp_imaginary[n0] +=Jx*psi_imaginary[n3]-Jy*psi_imaginary[n3]+Jz*psi_imaginary[n0];
          // psi_tmp_real[n1]      +=Jx*psi_real[n2]+Jy*psi_real[n2]-Jz*psi_real[n1];
          // psi_tmp_imaginary[n1] +=Jx*psi_imaginary[n2]+Jy*psi_imaginary[n2]-Jz*psi_imaginary[n1];
          // psi_tmp_real[n2]      +=Jx*psi_real[n1]+Jy*psi_real[n1]-Jz*psi_real[n2];
          // psi_tmp_imaginary[n2] +=Jx*psi_imaginary[n1]+Jy*psi_imaginary[n1]-Jz*psi_imaginary[n2];
          // psi_tmp_real[n3]      +=Jx*psi_real[n0]-Jy*psi_real[n0]+Jz*psi_real[n3];
          // psi_tmp_imaginary[n3] +=Jx*psi_imaginary[n0]-Jy*psi_imaginary[n0]+Jz*psi_imaginary[n3];



          // test=k;
        }
      }
    }
  }
  for (int i = 0; i < nofstates; i++) {
    psi_real[i]=1*psi_real[i]+(-(1./Temperature)*psi_tmp_real[i]/(2.*M));
    psi_imaginary[i]=1*psi_imaginary[i]+(-(1./Temperature)*psi_tmp_imaginary[i]/(2.*M));
  }
  double norm=0;
  for (int i = 0; i < nofstates; i++)
    norm+=psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i];

  norm=1./sqrt(norm);
  for (int i = 0; i < nofstates; i++) {
    psi_real[i]=norm*psi_real[i];
    psi_imaginary[i]=norm*psi_imaginary[i];
  }

}


/* Run the simulation of TDSE with random wave function.
  Input:

  Side Effect:
    array_real[],array_imagine[]
  #: usually the side effect is on psi_real[] and psi_imaginary[]
*/
void spin_system::random_wavef_run(){
  int total_steps=(int) T/tau;
  double* frequency= new double [total_steps+1]();

  double norm=0;

  // random_product(psi_real,psi_imaginary,psi_sys_real,psi_sys_imaginary);
  for (int i = 0; i < nofstates; i++)
    coefficient_return[i]+=psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i];

  for (int i = 0; i < nofstates; i++)
    norm+=coefficient_return[i];

  cout<<"norm before run: "<<norm<<endl;
  cout<<"START RUNNING!"<<endl;

  //output
  ostringstream strs;
  strs <<"G"<<(G/10.)<<"_"<<"Ts"<<(T*10)<<".dat";//output name for single output

  string str = strs.str();
  string strmain="output_general_";
  strmain.append(str);
  const char *testChars = strmain.c_str();
  ofstream output(testChars);
  output<<"Time Energy_sys Energy_env Energy_se Energy_all Frequency ";
  for (int i = 0; i < N; i++)
    output<<"Sx_"<<i<<" "<<"Sy_"<<i<<" "<<"Sz_"<<i<<" ";
  output<<endl;
  //output set END

  for (int step = 0; step < total_steps+1; step++){ //+1 because i count the 0 point and the last poing as well.
    Delta=step*tau/T;
    Gamma=1-Delta;
    if (step%2500==0)
      cout<<"step:"<<step<<endl;

    energy_sys_return[step]=energy_sys(step*tau);
    energy_env_return[step]=energy_env(step*tau);
    energy_se_return[step]=energy_se(step*tau);
    energy_all_return[step]=energy_all(step*tau);
    // for (int s = 0; s < N; s++) {
    //   int index=step*N*3+s*3;
    //   spin_return[index]  =spin('x',s);
    //   spin_return[index+1]=spin('y',s);
    //   spin_return[index+2]=spin('z',s);
    // }

    // output start
    output<<step*tau<<" ";
    output<<energy_sys_return[step]<<" ";
    output<<energy_env_return[step]<<" ";
    output<<energy_se_return[step]<<" ";
    output<<energy_all_return[step]<<" ";
    output<<frequency[step]<<" ";
    for (int i = 0; i < 3*N; i++)
      output<<spin_return[step*3*N+i]<<" ";
    output<<endl;

    for (int i = gs_sol; i < nofstates; i+=256)
      frequency[step]+=psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i];

    single_spin_op(step*tau);
    double_spin_op_x(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_z(step*tau);
    double_spin_op_y(step*tau);
    double_spin_op_x(step*tau);
    single_spin_op(step*tau);
  }
  for (int i = 0; i < nofstates; i++)
    coefficient_return[i]=psi_real[i]*psi_real[i]+psi_imaginary[i]*psi_imaginary[i];

  cout<<"END RUNNING"<<endl;
  norm=0.;
  for (int i = 0; i < nofstates; i++)
    norm+=coefficient_return[i];
  cout<<"norm after run: "<<norm<<endl;
  // output the return value: coefficient, energy expectation value, and spin expectation value.
  ofstream Coefficient_out("output_coefficient.dat");
  for (int i = 0; i < nofstates; i++)
    Coefficient_out<<coefficient_return[i]<<endl;

  // for (int step = 0; step < total_steps+1; step++){
  //   output<<step*tau<<" ";
  //   output<<energy_sys_return[step]<<" ";
  //   output<<energy_env_return[step]<<" ";
  //   output<<energy_se_return[step]<<" ";
  //   output<<energy_all_return[step]<<" ";
  //   output<<frequency[step]<<" ";
  //   for (int i = 0; i < 3*N; i++) {
  //     output<<spin_return[step*3*N+i]<<" ";
  //   }
  //   output<<endl;
  // }
  success_probability_return=frequency[total_steps];

}
