function [is_mpc_always_feasible, is_disease_under_threshold, is_control_too_quick] = ...
    check_constraints( x_mpc, control, tau_mpc, ga_flags, MODEL_PARAM, CONTROL_PARAM )

% INPUT EXTRACTION --------------------------------------------------------
H_max = MODEL_PARAM.H_max;
M     = MODEL_PARAM.M;

T       = CONTROL_PARAM.T;
tau_min = CONTROL_PARAM.tau_min;

% CHECK FEASIBILITY OF MPC ------------------------------------------------
infeasible_instants = find(ga_flags(:, 2) <= 0);
if ~isempty(infeasible_instants)
    is_mpc_always_feasible = true;
    disp([ 'WARNING: Could not solve MPC ' num2str(length(infeasible_instants)), ...
        '/', num2str(size(ga_flags, 1)), ' times.']);
else
    is_mpc_always_feasible = false;
    disp('OK: could always solve MPC');
end

% CHECK DIESEASE CONTROL --------------------------------------------------
is_disease_under_threshold = false(1, M);
for i = 1 : M
    H_i = x_mpc( :, 3*M + i );
    if max(H_i) <= H_max(i)
        is_disease_under_threshold(i) = true;
    else
        is_disease_under_threshold(i) = false;
    end 
end
if all(is_disease_under_threshold)
    disp('OK: Hospitalized below threshold for all regions.');
else
    disp( ['PROBLEM: Infected over threshold in regions ', num2str(find(~is_disease_under_threshold)), '.' ] ); 
end

% CHECK QUICKNESS OF CONTROL ----------------------------------------------
is_control_too_quick = false;
control_too_quick_instants = [];
for k = 2 : T
    if ~isequal(control(k, :), control(k-1, :) )
       if tau_mpc(k) < tau_min
           is_control_too_quick       = true; 
           control_too_quick_instants = [control_too_quick_instants; k ]; %#ok<AGROW>
       end
    end
end
if is_control_too_quick
    disp('PROBLEM: Control changed to quickly.');
else
    disp('OK: Control changed slowly enough.');
end

end