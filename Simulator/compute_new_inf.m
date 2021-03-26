function new_inf = compute_new_inf(S_bar, new_inf_1, MODEL_PARAM)

M = MODEL_PARAM.M;

T = size(S_bar, 1);

new_inf = zeros(T, M);
new_inf(1, :) = new_inf_1(end, :);
for k = 2 : T
    S_k = S_bar(k, :)';
    S_k_m_1 = S_bar(k-1, :)';
    new_inf(k, :) = update_new_inf(S_k, S_k_m_1, MODEL_PARAM);
end