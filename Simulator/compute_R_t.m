function R_t_bar = compute_R_t(new_inf, new_inf_1, MODEL_PARAM)

gener_interv = MODEL_PARAM.gener_interv;
M            = MODEL_PARAM.M;

new_inf_shift = [new_inf_1; new_inf(2:end, :)];
T = size(new_inf, 1);

R_t_bar = -ones(T, M);

k_start = 0;
is_found_start_point = 0;
while ~is_found_start_point && k_start < T
    k_start = k_start + 1;    
    if all(new_inf_shift(k_start, :) ~= -1)
        is_found_start_point = 1;
    end
end
    
for k = k_start : T
    h = k + 2*gener_interv - 1;
    last_gi_days = h-gener_interv+1 : h;
    previous_gi_days = 1 + h - 2*gener_interv : h-gener_interv;
    %R_t(k, :) = new_inf_shift(last_gi_days(end), :) ./ new_inf_shift(previous_gi_days(end), :);     % V1
    R_t_bar(k, :) = mean(new_inf_shift(last_gi_days, :)) ./ mean(new_inf_shift(previous_gi_days, :));   % V2
    R_t_bar(k, isnan(R_t_bar(k, :))) = -1;
    R_t_bar(k, isinf(R_t_bar(k, :))) = -1;
end


% length of new_inf_1 changes according to whether you use V1 or V2