function tau_k = update_tau( tau_k_m_1, control_k_m_1, control_k_m_2, tau_min )
% k_m_1 means k minus 1, k_m_2 means k minus 2

if isequal(control_k_m_1, control_k_m_2)
    tau_k = min( tau_k_m_1 + 1, tau_min);
else
    tau_k = 1;
end

end