function contraction_constraints = compute_contraction_constraints_global(...
    x_bar, control_bar, T_p, is_global_trigger_met, MODEL_PARAM, CONTROL_PARAM)
% returns T_p elements in column

% INPUT EXTRACTION --------------------------------------------------------
M = MODEL_PARAM.M;
norm_bound = CONTROL_PARAM.norm_bound;

% COMPUTATION -------------------------------------------------------------
contraction_constraints = zeros(T_p, 1);
for k = 1 : T_p
    if is_global_trigger_met(k)
        Psi = compute_Psi(x_bar(k,:), control_bar(k, :), MODEL_PARAM);
        dynamics_I = eye(M) + Psi;
        norm_rows = sum( abs(dynamics_I), 2);
        
        contraction_constraints(k) = max(norm_rows) - norm_bound;
    else
        contraction_constraints(k) = 0;
    end
end

end

