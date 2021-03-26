% PROGRAM INITIALIZATION --------------------------------------------------
clear
%clc;   close all
addpath('utilities');
addpath('parameters');
disp('Simulation started...');
rng default              % for reproducibility
comp_time_start = tic;   % start timer
n_x = 6;                 % don't change: dimension of state of single region (i.e. num. of compartments)


% USER INTERFACE ========================================================== 
% CONTROL PARAMETERS ------------------------------------------------------
load parameters_macroregions_2021-25-02;
% M: number of regions

x_1             = reshape([N(1)-(I_0(1)+Q_0(1)+H_0(1)+R_0(1)+D_0(1)), I_0(1), Q_0(1), H_0(1), R_0(1), D_0(1);  ...   % initial state for epidemic model [(S, I, Q, U, D, R) region 1, (S, I, Q, U, D, R) region 2 ...].
                           N(2)-(I_0(2)+Q_0(2)+H_0(2)+R_0(2)+D_0(2)), I_0(2), Q_0(2), H_0(2), R_0(2), D_0(2);  ...
                           N(3)-(I_0(3)+Q_0(3)+H_0(3)+R_0(3)+D_0(3)), I_0(3), Q_0(3), H_0(3), R_0(3), D_0(3)], ...
                               1, M*n_x);
contraction_method = 'local';        % 'global' or 'local'; how contraction is applied

rho_pool        = [0.3 0.4 0.5 0.6 0.7]; % values that rho can take
varphi_pool     = [sqrt(0.3) 1];     % values that varphi can take
sigma_pool      = [0 0.5 1];                 % values that sigma can take
epsilon_H       = 0.3;               % percentage of ICUs that triggers constraint on contraction
epsilon_R       = 1.3;               % value of R_t that triggers constraint on contraction
discount_factor = 0.9;               % discount in cost function. 0 < ... <= 1
norm_bound      = 0.99;              % < 1

T               = 365;               % overall simulation time
tau_min         = 14;                % control must be equal on tau_min days before it can be changed
T_p             = 29;                % prediction horizon for MPC
MPC_latency     = 5;                 % (>= 1) time instants MPC control is used if deemed fit by 1 MPC iteration (it's an approximation of a classical MPC)

is_plot_mpc_sessions   = false;      % true: plots each session in mpc
is_plot_ga_generations = false;      % true: plots genetic algortihm generations

    rho_pool_size    = length(rho_pool);    % do not change
    varphi_pool_size = length(varphi_pool); % do not change
    sigma_pool_size  = length(sigma_pool);  % do not change
    
% MODEL PARAMETERS --------------------------------------------------------
gener_interv = 4;      % (integer) generation interval (average days required for contagion after infection) [Kohlberg, Neyman, "Demystifying the math of the coronavirus", 2020.]
c_w          = 1 - 0.383; % fraction of people who cannot work with social distancing [Holgersen, Jia, Svenkerud 2020]
c_m          = 83.992;  % daily gdp per capita in euro [nominal GDP per capita 30657€ "World Economic Outlook Database, October 2020". IMF.org. International Monetary Fund. Retrieved 22 October 2020]

alpha_0     = 0.0671 * ones(M, 1);  % additional alpha due to increased testing 
alpha_tilde = (0.1477 - 0.0671) * ones(M, 1);  % additional alpha due to increased testing 
cost_of_1_test                                  = 59;    % in euros
tests_per_person_tested                         = 1.74;   % normally each person tested gets more than one test
people_tested_daily_over_popul_when_max_testing = 3.49e-3;

n_windows = 1 + ceil( (T_p-1) / tau_min);    % maximum possible number of control windows in T_p

% GENETIC ALGORITHM OPTIONS -----------------------------------------------
ga_population_size       = 4000;
ga_max_stall_generations = 50;                              % default: 50
ga_elite_count           = floor(0.2 * ga_population_size); % default: 0.05 * ga_population_size
ga_crossover_fraction    = 0.4;                             % default: 0.8
ga_constraint_tolerance  = 1e-4; 
ga_display               = 'off';                           % 'off' or 'final'
ga_use_parallel          = true;                            % default: false
%==========================================================================


% PREPARATION =============================================================
% CONSISTENCY CHECK -------------------------------------------------------
assert( all( size(x_1) == [1, n_x*M]), 'ERROR: Initial state has incorrect size.');

% COMPUTE DISCOUNT FACTORS ------------------------------------------------
discounts = compute_discounts(discount_factor, T_p);

% DECISION VARIABLES' STRUCTURE -------------------------------------------


size_win_len          = n_windows;           % win len is lengths of control windows
size_rho_level_win    = n_windows * M;       % rho_level_win is level of rho in a window (level is an index for rho_pool)
size_varphi_level_win = n_windows * M;
size_sigma_level_win  = n_windows * M;
% rho_level_win is [level_win1_region1, level_win2_region1, ..., level_win1_region2, level_win2_region2, ...]

idx_win_len          = 1 : size_win_len;     % indices of win_len in u_bar
idx_rho_level_win    = size_win_len + (1 : size_rho_level_win);
idx_varphi_level_win = size_win_len + size_rho_level_win + (1 : size_varphi_level_win);
idx_sigma_level_win  = size_win_len + size_rho_level_win + size_varphi_level_win + (1 : size_sigma_level_win);

% GENETIC ALGORITHM OPTIONS -----------------------------------------------
if is_plot_ga_generations == true,    ga_plot_fcn = @gaplotbestf;
else,    ga_plot_fcn = [];    end

ga_options = optimoptions( 'ga', ...
    'PopulationSize',      ga_population_size,    ...
    'MaxStallGenerations', ga_max_stall_generations, ...
    'EliteCount',          ga_elite_count, ...
    'CrossoverFraction',   ga_crossover_fraction, ...
    'ConstraintTolerance', ga_constraint_tolerance, ...
    'Display',             ga_display, ...
    'PlotFcn',             ga_plot_fcn, ...
    'UseParallel',         ga_use_parallel );

% INITIALIZATION OF VARIABLES ---------------------------------------------
N   = sum(reshape(x_1, M, n_x), 2);
c_alpha = people_tested_daily_over_popul_when_max_testing * ... % daily cost in euro per person to do maximum testing
          tests_per_person_tested * cost_of_1_test;

ctrl_size_rho    = size_rho_level_win / n_windows;    % number of rho variables in "control variable" (= M)
ctrl_size_varphi = size_varphi_level_win / n_windows;
ctrl_size_sigma  = size_sigma_level_win / n_windows;
ctrl_idx_rho     = 1 : ctrl_size_rho;                 % indices of rho variables in "control variable"
ctrl_idx_varphi  = ctrl_size_rho + (1 : ctrl_size_varphi);
ctrl_idx_sigma   = ctrl_size_rho + ctrl_size_varphi + (1 : ctrl_size_sigma);

x_mpc       = zeros(T, length(x_1));
new_inf     = zeros(T, M);   % new_inf_k = (N-S_k - (N-S_(k-1))
control_mpc = zeros(T, (ctrl_size_rho + ctrl_size_varphi + ctrl_size_sigma));
tau_mpc     = zeros(T, 1);   % tau(i) is the number of days before i that the control stayed the same
J_mpc       = zeros(T, 1);   % control cost computed by mpc
ga_flags    = [];            % flag given by genetic algorithm regarding success

new_inf_1 = -ones(2*gener_interv, M);  % up to time instant 1 included
rho_0     = -ones(1, M);   % value of rho on last day before MPC starts. all -1 means not existing
varphi_0  = -ones(1, M);   
sigma_0   = -ones(1, M);   
tau_1     = 0;             % number of days that had same value of control before MPC starts
MPC_lag   = 0;

x_mpc(1, :)   = x_1;
new_inf(1, :) = new_inf_1(end, :);
tau_mpc(1)    = tau_1;
control_0     = [rho_0, varphi_0, sigma_0];

% STRUCTS ASSEMBLY --------------------------------------------------------
MODEL_PARAM = struct;           
MODEL_PARAM.alpha            = alpha;
MODEL_PARAM.alpha_0          = alpha_0;
MODEL_PARAM.alpha_tilde      = alpha_tilde;
MODEL_PARAM.beta             = beta;
MODEL_PARAM.eta_H            = eta_H;
MODEL_PARAM.eta_Q            = eta_Q;
MODEL_PARAM.gamma            = gamma;
MODEL_PARAM.gener_interv     = gener_interv;
MODEL_PARAM.H_max            = H_max;
MODEL_PARAM.ICU_max          = ICU_max + ICU_att;
MODEL_PARAM.c_alpha          = c_alpha;
MODEL_PARAM.c_m              = c_m;
MODEL_PARAM.c_w              = c_w;
MODEL_PARAM.ctrl_idx_rho     = ctrl_idx_rho;
MODEL_PARAM.ctrl_idx_varphi  = ctrl_idx_varphi;
MODEL_PARAM.ctrl_idx_sigma   = ctrl_idx_sigma;
MODEL_PARAM.ctrl_size_rho    = ctrl_size_rho;
MODEL_PARAM.ctrl_size_varphi = ctrl_size_varphi;
MODEL_PARAM.ctrl_size_sigma  = ctrl_size_sigma;
MODEL_PARAM.kappa_H          = kappa_H;
MODEL_PARAM.kappa_Q          = kappa_Q;
MODEL_PARAM.N                = N;
MODEL_PARAM.n_x              = n_x;
MODEL_PARAM.M                = M;
MODEL_PARAM.Phi              = Phi;
MODEL_PARAM.psi              = psi;
MODEL_PARAM.zeta_baseline    = zeta_baseline;
MODEL_PARAM.zeta_slope       = zeta_slope;

CONTROL_PARAM = struct;
CONTROL_PARAM.discounts          = discounts;
CONTROL_PARAM.epsilon_H          = epsilon_H;
CONTROL_PARAM.epsilon_R          = epsilon_R;
CONTROL_PARAM.contraction_method = contraction_method;
CONTROL_PARAM.MPC_lag_max        = MPC_latency;
CONTROL_PARAM.norm_bound         = norm_bound; 
CONTROL_PARAM.rho_pool           = rho_pool;
CONTROL_PARAM.sigma_pool         = sigma_pool;
CONTROL_PARAM.T                  = T;
CONTROL_PARAM.T_p                = T_p;
CONTROL_PARAM.tau_min            = tau_min;
CONTROL_PARAM.varphi_pool        = varphi_pool;

DECISION_VAR = struct;
DECISION_VAR.n_windows             = n_windows;
DECISION_VAR.size_win_len          = size_win_len;
DECISION_VAR.size_rho_level_win    = size_rho_level_win;
DECISION_VAR.size_varphi_level_win = size_varphi_level_win;
DECISION_VAR.size_sigma_win        = size_sigma_level_win;
DECISION_VAR.idx_win_len           = idx_win_len;
DECISION_VAR.idx_rho_level_win     = idx_rho_level_win;
DECISION_VAR.idx_varphi_level_win  = idx_varphi_level_win;
DECISION_VAR.idx_sigma_win         = idx_sigma_level_win;
DECISION_VAR.rho_pool_size         = rho_pool_size;
DECISION_VAR.varphi_pool_size      = varphi_pool_size;
DECISION_VAR.sigma_pool_size       = sigma_pool_size;

PREV_STATE = struct;

MPC_OPT = struct;                  
MPC_OPT.ga_options           = ga_options;
MPC_OPT.is_plot_mpc_sessions = is_plot_mpc_sessions;
% =========================================================================


% EXECUTION ===============================================================
% MPC AT STEP 1 -----------------------------------------------------------
comp_time_mpc_start  = tic;
PREV_STATE.x_1       = x_mpc(1, :);
PREV_STATE.t_1       = 1;
PREV_STATE.tau_1     = tau_mpc(1);
PREV_STATE.control_0 = control_0;
PREV_STATE.MPC_lag   = MPC_lag;
PREV_STATE.new_inf_1 = new_inf_1;

[control_opt, MPC_lag, J_mpc(1), ga_flag, ga_outputs(1)] = run_mpc(MPC_OPT, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE); 

control_mpc(1, :) = control_opt(1, :);
x_mpc(2, :)       = epidemic_model(x_mpc(1, :), control_mpc(1, :), MODEL_PARAM);
tau_mpc(2)        = update_tau(tau_mpc(1), control_mpc(1, :), control_0, tau_min);
MPC_lag           = max(MPC_lag-1, 0);
new_inf(2, :)     = update_new_inf(x_mpc(2, 1:M)', x_mpc(1, 1:M)', MODEL_PARAM);

ga_flags = [ga_flags; 1, ga_flag];  % append ga_flag at time 1
if ga_flag > 0,    disp(['Step 1/', num2str(T), ' completed with MPC (feasible).']);
else,    disp(['Step 1/', num2str(T), ' completed with MPC (infeasible).']); end
comp_time_mpc_episode = toc(comp_time_mpc_start);
disp_computation_time(comp_time_mpc_episode);

% MPC AT STEPS >= 2 -------------------------------------------------------
for k = 2 : T   
    if tau_mpc(k) >= tau_min && MPC_lag < 1      % if control can change, run mpc
        comp_time_mpc_start  = tic;
        PREV_STATE.x_1       = x_mpc(k, :);
        PREV_STATE.t_1       = 1;
        PREV_STATE.tau_1     = tau_mpc(k);
        PREV_STATE.control_0 = control_mpc(k-1, :);
        PREV_STATE.MPC_lag   = MPC_lag;
        if k >= 2*gener_interv
            PREV_STATE.new_inf_1 = new_inf( (1 + k - 2*gener_interv : k), :);
        else
            PREV_STATE.new_inf_1 = [ new_inf_1( 1 + (end-1)-(2*gener_interv-k): (end-1), :) ; new_inf( 1:k, :) ];
        end
        
        [control_opt, MPC_lag, J_mpc(k), ga_flag, ga_outputs(k)] = run_mpc(MPC_OPT, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE); %#ok<SAGROW>
        
        control_mpc(k, :) = control_opt(1, :);
        
        ga_flags = [ga_flags; k, ga_flag]; %#ok<AGROW>
        if ga_flag > 0,    disp(['Step ', num2str(k), '/', num2str(T), ' completed with MPC (feasible).']);
        else,    disp(['Step ', num2str(k), '/', num2str(T), ' completed with MPC (infeasible).']); end
        comp_time_mpc_episode = toc(comp_time_mpc_start);
        disp_computation_time(comp_time_mpc_episode);
    else          % if control cannot change use previous value for control
        control_mpc(k, :) = control_mpc(k-1, :);
        
        if tau_mpc(k) >= tau_min && MPC_lag >= 1
            disp(['Step ', num2str(k), '/', num2str(T), ' completed without MPC (latency).']);
        else
            disp(['Step ', num2str(k), '/', num2str(T), ' completed without MPC.']);    
        end
    end
    
    x_mpc(k + 1, :)    = epidemic_model(x_mpc(k, :), control_mpc(k, :), MODEL_PARAM);
    if k < T
        tau_mpc(k + 1) = update_tau(tau_mpc(k), control_mpc(k, :), control_mpc(k-1, :), tau_min);
        new_inf(k + 1, :) = update_new_inf( x_mpc(k+1, 1:M)', x_mpc(k, 1:M)', MODEL_PARAM );
    end
    MPC_lag = max(MPC_lag-1, 0);
end
R_t_mpc = compute_R_t(new_inf, new_inf_1, MODEL_PARAM);

% COMPUTE COST ------------------------------------------------------------
J_step = zeros(T, 1);     % cost from one time instant
for k = 1 : T
    J_step(k) = cost_function_one_step( x_mpc(k, :), control_mpc(k, :), MODEL_PARAM );
end
J_step(T) = cost_function_last_step( );
J = sum(J_step); 

comp_time = toc(comp_time_start);
%==========================================================================


% RESULTS =================================================================
% REPORT ------------------------------------------------------------------
disp('---------- REPORT ----------');
check_constraints(x_mpc, control_mpc, tau_mpc, ga_flags, MODEL_PARAM, CONTROL_PARAM);
disp(['Final cost is ', num2str(J/1e9), ' bln €']);
disp_computation_time(comp_time);
disp('----------------------------');

% PLOT --------------------------------------------------------------------   
PREV_STATE.t_1       = 1;
PREV_STATE.control_0 = control_0;
plot_epidemic_model(x_mpc, control_mpc, tau_mpc, R_t_mpc, MODEL_PARAM, CONTROL_PARAM, PREV_STATE);
plot_costs(J_step);
%==========================================================================

rmpath('utilities');
 