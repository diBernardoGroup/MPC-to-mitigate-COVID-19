function [nonlin_con_ineq, nonlin_con_eq] = nonlin_con(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE)
% Note that x_bar and rho_bar have time on along column, regions along row

% INPUT EXTRACTION --------------------------------------------------------
ICU_max = MODEL_PARAM.ICU_max;
M       = MODEL_PARAM.M;

epsilon_H          = CONTROL_PARAM.epsilon_H;
epsilon_R          = CONTROL_PARAM.epsilon_R;
contraction_method = CONTROL_PARAM.contraction_method;
T_p                = CONTROL_PARAM.T_p;
tau_min            = CONTROL_PARAM.tau_min;

idx_win_len = DECISION_VAR.idx_win_len;

control_0 = PREV_STATE.control_0;
new_inf_1 = PREV_STATE.new_inf_1;

% COMPUTE CONTROL ---------------------------------------------------------
win_len       = u_bar(idx_win_len);
[control_bar] = u_bar_2_control(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR);

% SIMULATE EPIDEMIC MODEL -------------------------------------------------
x_bar   = predict_epidemic_model(control_bar, MODEL_PARAM, PREV_STATE);
ICU_bar = 0.1 * x_bar( 1:T_p, ( 3*M + (1:M) ) );

S_bar   = x_bar( 1:T_p, 1:M );
new_inf = compute_new_inf(S_bar, new_inf_1, MODEL_PARAM);
R_t_bar = compute_R_t(new_inf, new_inf_1, MODEL_PARAM);

% COMPUTE JACOBIAN MATRIX NORM --------------------------------------------
is_local_trigger_met  = (ICU_bar >= epsilon_H * ICU_max)  |  (R_t_bar >= epsilon_R);
is_global_trigger_met = any(is_local_trigger_met, 2);

switch contraction_method
    case 'global'
        contraction_constraints = compute_contraction_constraints_global(...
            x_bar, control_bar, T_p, is_global_trigger_met, MODEL_PARAM, CONTROL_PARAM);
    case 'local'
        contraction_constraints = compute_contraction_constraints_local(...
            x_bar, control_bar, T_p, is_local_trigger_met, MODEL_PARAM, CONTROL_PARAM);
    otherwise
        error('contraction_method is invalid.')
end
contr_con_size = length(contraction_constraints);

% ASSEMBLE CONSTRAINTS ----------------------------------------------------
nonlin_con_eq = [];  % cannot give eq constraints with integer constraints
nonlin_con_ineq = zeros(contr_con_size + 1, 1);

nonlin_con_ineq(1 : contr_con_size) = contraction_constraints;

if ~isequal( control_bar(1, :), control_0) || isequal( control_0, -ones(1, length(control_0)))
    nonlin_con_ineq(end) = tau_min - win_len(1);  % first window must be at least tau_min
else
    nonlin_con_ineq(end) = 0;  % first window can be any length
    % NB: this is ok because MPC is only called when tau_1 == tau_min, then the
    % first window is effectively a continuation of the previous one (before MPC)
end

end