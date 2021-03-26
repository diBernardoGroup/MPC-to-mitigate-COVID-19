function J = cost_function( u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR, PREV_STATE )

% EXTRACT INPUT -----------------------------------------------------------
discounts = CONTROL_PARAM.discounts;
T_p       = CONTROL_PARAM.T_p;

% COMPUTE CONTROL ---------------------------------------------------------
control_bar = u_bar_2_control(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR);
x_bar       = predict_epidemic_model(control_bar, MODEL_PARAM, PREV_STATE);

% COMPUTE COST FUNCTION ---------------------------------------------------
J_step = zeros(T_p, 1);     % cost from one time instant
for k = 1 : T_p-1
    J_step(k) = cost_function_one_step( x_bar(k, :), control_bar(k, :), MODEL_PARAM );
end
J_step(T_p)   = cost_function_last_step(  );

J = sum(discounts .* J_step);            % cost from all prediction horizon

end