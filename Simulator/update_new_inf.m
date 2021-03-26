function new_inf_k = update_new_inf(S_k, S_k_m_1, MODEL_PARAM)

N = MODEL_PARAM.N;

infected_k = N - S_k;
infected_k_m_1 = N - S_k_m_1;

new_inf_k = infected_k - infected_k_m_1;

end