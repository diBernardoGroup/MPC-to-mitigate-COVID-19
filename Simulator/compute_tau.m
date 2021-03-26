function tau = compute_tau( control, tau_1, tau_min, control_0 )

T_p = size(control, 1);
tau = zeros(T_p, 1);

tau(1) = tau_1;
if T_p >= 2
    tau(2) = update_tau(tau(1), control(1, :), control_0, tau_min);
end
for k = 3 : T_p
    tau(k) = update_tau(tau(k-1), control(k-1, :), control(k-2, :), tau_min);
end
