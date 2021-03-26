function [control_bar_opt, MPC_lag, J_opt, ga_flag, ga_output] = run_mpc(MPC_OPT, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE)

% INPUT EXTRACTION --------------------------------------------------------
M = MODEL_PARAM.M;

MPC_lag_max = CONTROL_PARAM.MPC_lag_max;
T_p         = CONTROL_PARAM.T_p;
tau_min     = CONTROL_PARAM.tau_min;

idx_rho_level_win     = DECISION_VAR.idx_rho_level_win;
idx_sigma_level_win   = DECISION_VAR.idx_sigma_win;
idx_varphi_level_win  = DECISION_VAR.idx_varphi_level_win;
idx_win_len           = DECISION_VAR.idx_win_len;
rho_pool_size         = DECISION_VAR.rho_pool_size;
sigma_pool_size       = DECISION_VAR.sigma_pool_size;
size_rho_level_win    = DECISION_VAR.size_rho_level_win;
size_sigma_level_win  = DECISION_VAR.size_sigma_win;
size_varphi_level_win = DECISION_VAR.size_varphi_level_win;
size_win_len          = DECISION_VAR.size_win_len;
varphi_pool_size      = DECISION_VAR.varphi_pool_size;

control_0 = PREV_STATE.control_0;
new_inf_1 = PREV_STATE.new_inf_1;

is_plot_mpc_sessions = MPC_OPT.is_plot_mpc_sessions;
ga_options           = MPC_OPT.ga_options;

% LINEAR CONSTRAINT -------------------------------------------------------
n_vars           = size_win_len + size_rho_level_win + size_varphi_level_win + size_sigma_level_win;  % number of decision variable for genetic algorithm
int_vars         = [idx_win_len idx_rho_level_win idx_varphi_level_win idx_sigma_level_win];     % which decision variables are integer
lin_con_ineq_mat = [];    % matrix for linear inequality constraints: A in A*u_bar <= b
lin_con_ineq_vec = [];    % vector for linear inequality constraints: b in A*u_bar <= b
lin_con_eq_mat   = [];    % matrix for linear equality constraints: A in A*u_bar = b
lin_con_eq_vec   = [];    % matrix for linear equality constraints: b in A*u_bar = b
low_bound        = [1; tau_min*ones(size_win_len-1, 1); ...
                    ones(size_rho_level_win, 1); ...
                    ones(size_varphi_level_win, 1); ...
                    ones(size_sigma_level_win, 1) ];  % LB in u_bar >= LB
up_bound         = [T_p*ones(size_win_len, 1); ...
                    rho_pool_size*ones(size_rho_level_win, 1); ...
                    varphi_pool_size*ones(size_varphi_level_win, 1); ...
                    sigma_pool_size*ones(size_sigma_level_win, 1) ];     % UB in u_bar <= UB

% NONLINEAR CONSTRAINT FUNCTION -------------------------------------------
nonlin_con_handle = @(u_bar) nonlin_con(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE);    % nonlinear inequality and equality constraints

% COST FUNCTION -----------------------------------------------------------
cost_fun_handle = @(u_bar) cost_function(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE);    % cost function

% MPC EXECUTION -----------------------------------------------------------
[u_bar_opt, J_opt, ga_flag, ga_output] = ga(cost_fun_handle, n_vars, lin_con_ineq_mat, lin_con_ineq_vec, ...
       lin_con_eq_mat, lin_con_eq_vec, low_bound, up_bound, nonlin_con_handle, int_vars, ga_options );

% COMPUTE OUTPUT ----------------------------------------------------------
[control_bar_opt] = u_bar_2_control(u_bar_opt, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR);

first_win_len = u_bar_opt(1);
MPC_lag = min(first_win_len, MPC_lag_max);  % MPC latency is length of first window up to MPC_lag_max

% PLOT ONE MPC STEP -------------------------------------------------------
if is_plot_mpc_sessions
    x_bar_opt   = predict_epidemic_model(control_bar_opt, MODEL_PARAM, PREV_STATE);
    tau_1       = tau_min;   % if mpc was run, we can assume that tau_1 was equal to tau_min (except for first step, where tau_1 = and rho_last = -1)
    tau_bar_opt = compute_tau(control_bar_opt, tau_1, tau_min, control_0);
    new_inf_opt = compute_new_inf(x_bar_opt(:, 1:M), new_inf_1, MODEL_PARAM);
    R_t_opt     = compute_R_t(new_inf_opt, new_inf_1, MODEL_PARAM);
    plot_epidemic_model(x_bar_opt, control_bar_opt, tau_bar_opt, R_t_opt, MODEL_PARAM, CONTROL_PARAM, PREV_STATE);    
    input('Type Enter to continue...');
end

end

