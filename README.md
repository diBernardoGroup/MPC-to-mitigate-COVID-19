# MPC-to-mitigate-COVID-19
This repository contains the code used in the paper: 

> M. Coraggio, S. Xie, F. De Lellis, G. Russo, M. di Bernardo, "Intermittent non-pharmaceutical strategies to mitigate the COVID-19 epidemic in a network model of Italy via constrained optimization"

for more detail please read the paper at arxiv.org/abs/2103.16502


**SYSTEM REQUIREMENTS**

All the scripts listed were developed using MATLAB R2020b on Windows Server 2019 (64-bit) OS.
The hardware used for simulations is an Intel Xeon Gold 6226 chip and 32Gb RAM DDR4.
With this system, each of the simulations takes approximately one hour and a half. 


**REPOSITORY ORGANIZATION** 

The repository is organized as follows:
1. The folders 'National_Identification' and 'Regional_Identification' contain the scripts for the regional and national identification procedures;
2. Both the above identification folders contain a .txt that can be used to read the data off-line (read_national_data and read_regional_data provide the option to also read the latest data from 'Protezione Civile' GitHub repository-see the corresponding files for more details);
3. Both the above identification folders contain a .mat that can be used to run directly the stage3.m (this .mat file is generated as an output from stage2.m therefore it is not needed if one wants to run the identification from stage1.m);
4. The folder 'Regression' contains the scripts and the data files to run the regression algorithm;
5. The scripts related to the network model are located in the folder 'Simulator'.

**READ ME CONTENTS ORGANIZATION**

This READ ME is organized as follows:
1. Section 'Simulator Description' includes the information about the simulator used for the model and its functionality;
2. Section 'Identification Description' includes the scripts implementing the identification procedure of the paper and their functionality;
3. Section 'Regression Description' includes the instruction to repeat the regression procedure.
----------------------------------------------------------------------------------------------------------------------------------------------------
# Simulator Description

The program is run simply by running main.m

Parameters to change the scenario to be run are at the beginning of main.m and are described by comments in the code; the name of variables is taken from the corresponding variables in the paper.

To run scenario 1 in the paper, set the following:
1. contraction_method = 'global';
2. sigma_pool         = 0;
3. epsilon_H          = 0; 
4. epsilon_R          = 0.   

To run scenario 2 in the paper, set the following:
1. contraction_method = 'local';
2. sigma_pool         = 0;
3. epsilon_H          = 0.3; 
4. epsilon_R          = 1.3.    

To run scenario 3 in the paper, set the following:
1. contraction_method = 'local';
2. sigma_pool         = [0 0.5 1];
3. epsilon_H          = 0.3; 
5. epsilon_R          = 1.3.    

The results of these simulations are already contained in the "simulations" folder.

----------------------------------------------------------------------------------------------------------------------------------------------------
# Identification Description

In the following, we give a brief explanation of what each file of the identification does.

'stage1_r.m' performs the identification of the time windows and of \rho*\beta, \tau, I_0 values in each window using
regional data. The identification is performed using an ad hoc nonlinear identification procedure.


'stage2_r.m' identifies \rho*\beta, \tau and I_0 values, given the time windows using an ad hoc nonlinear 
identification algorithm with regional data.


'stage3_r.m' identifies the remaining parameters using an ordinary least squares algorithm on regional data.
 
'stage2_r.m', 'stage1_r.m' you need to select:
1. the initial guess for each parameter

Moreover, in 'stage2_r.m' you need to select:

2. the time windows specified as an array whose elements are the starting points of each time window

'Regional_stage1.mat' contains the parameters obtained after stage 1 and stage 2 of the identification procedure described in the SI. They embed a numerical matrix organized as follows: Each row corresponds to a time window and each column represents a parameter. Specifically, the columns correspond, in order, to :
1. v (product between rho and beta);
2. tau;
3. gamma;
4. I0;
5. If.

'Regional_stage2.mat' contains the parameters obtained after stage 3 of the identification procedure described in the SI. They embed a numerical matrix organized as follows: Each column corresponds to a time window and each row represents a parameter. Specifically, the columns correspond, in order, to:
1. etaq;
1. etah;
2. zeta;
3. alpha;
4. psi;
5. kappaq;
6. kappa.

read_regional_data.m: 
Reads the regional data from the Protezione Civile GitHub repository. Here you need to select:
1. the code of the region you want to analyze;
2. the length of the averaging filter.


id_and_sim_r.m: 
Identifies and compares model predictions with data collected


INPUT:     
1. tab_data        (data vector to fit the model);
2. ti              (initial time instant);
3. te              (final time instant);
4. initial_guess    (initial guess for the identification);
5. N               (Number of residents);
6. total_active    (Function of currently active people).
OUTPUT: 
1. pars            (parameters of the model (rho*beta, tau, g, I0);
2. y               (model prediction);
3. If              (Final number of infected people).

Idendtify_model.m: 
Identifies the parameters of the model (nonlinear part)

INPUT:  
1. data            (data vector to fit the model);
2. lim_inf         (Inferior limit for the parameters);
3. lim_sup         (Superior limit for the parameters);
4. times           (Time instants for the identification);
5. initial_guess    (initial guess for the identification);
6. N               (Number of residents);
7. total_active    (Function of currently active people).

OUTPUT: pars            (Parameters identified from the algorithm)

Find_Change: 
finds if there is a breakpoint in the window and where it happened

INPUT:  
1. data            (data vector to fit the model);
2. fit1            (model prediction in the first half of the window);
3. fit2            (model prediction in the second half of the window);
4. fitTot          (model prediction in the entire window);
5. N               (Number of residents);
6. total_active    (Function of currently active people);
7. pr              (Parameters estimated in the entire window).
 

OUTPUT: 
1. Where           (Point where the parameter change happened);
2. Change          (Boolean advising a change in parameters).

 

Least_Squares_id
Runs the constrained Least square identification (linear part)


INPUT:  
1. In              (Infected time series);
2. total_quar      (Quarantined time series);
3. total_hosp      (Hospitalized time series);
4. eta             (Eta identified);
5. total_dead      (Dead time series);
6. tspan           (Time frame for the identification);
7. tau             (tau identified at stage 2).

OUTPUT: pars            (Parameters identified)

For more detailed information please refer to the folder https://github.com/diBernardoGroup/Network-model-of-the-COVID-19

**IDENTIFICATION PROCEDURE**

Regional Identification: In this section, we provide instructions to run the identification procedure for a given region

1. Select the region code you want to analyze in 'read_regional_data.m' (movemeanK); 
2. Select the averaging filter length in 'read_regional_data.m'; 
3. Select the initial guess for the parameters in the script 'stage1_r.m';
4. Run 'stage1_r.m';
5. Merge the windows following the procedure described in section S2 of the SI;
6. Run 'stage2_r.m' given the merged time windows obtained at point 5;
7. Run 'stage3_r.m'.
