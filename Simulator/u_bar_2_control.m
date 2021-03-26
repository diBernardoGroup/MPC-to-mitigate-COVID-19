function [control_bar] = u_bar_2_control(u_bar, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR)

% INPUT EXTRACTION --------------------------------------------------------
ctrl_idx_rho     = MODEL_PARAM.ctrl_idx_rho;
ctrl_idx_varphi  = MODEL_PARAM.ctrl_idx_varphi;
ctrl_idx_sigma   = MODEL_PARAM.ctrl_idx_sigma;
ctrl_size_rho    = MODEL_PARAM.ctrl_size_rho;
ctrl_size_varphi = MODEL_PARAM.ctrl_size_varphi;
ctrl_size_sigma  = MODEL_PARAM.ctrl_size_sigma;

T_p = CONTROL_PARAM.T_p;

idx_rho_level_win    = DECISION_VAR.idx_rho_level_win;
idx_varphi_level_win = DECISION_VAR.idx_varphi_level_win;
idx_sigma_level_win  = DECISION_VAR.idx_sigma_win;
idx_win_len          = DECISION_VAR.idx_win_len;

% COMPUTATION -------------------------------------------------------------
win_len          = u_bar(idx_win_len);         % lengths of control windows
rho_level_win    = u_bar(idx_rho_level_win);   % level of rho in each window
varphi_level_win = u_bar(idx_varphi_level_win);   % level of varphi in each window
sigma_level_win  = u_bar(idx_sigma_level_win);   % level of varphi in each window

[rho_bar, varphi_bar, sigma_bar] = windows_2_control(win_len, rho_level_win, varphi_level_win, sigma_level_win, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR);

control_bar = zeros(T_p, ctrl_size_rho + ctrl_size_varphi + ctrl_size_sigma );
control_bar(:, ctrl_idx_rho)    = rho_bar;
control_bar(:, ctrl_idx_varphi) = varphi_bar;
control_bar(:, ctrl_idx_sigma)  = sigma_bar;
end


function [rho_bar, varphi_bar, sigma_bar] = windows_2_control(win_len, rho_level_win, varphi_level_win, sigma_level_win, MODEL_PARAM, CONTROL_PARAM, DECISION_VAR)

M = MODEL_PARAM.M;

rho_pool    = CONTROL_PARAM.rho_pool;
sigma_pool  = CONTROL_PARAM.sigma_pool;
T_p         = CONTROL_PARAM.T_p;
varphi_pool = CONTROL_PARAM.varphi_pool;

n_windows = DECISION_VAR.n_windows;

rho_bar    = zeros(T_p, M);     
varphi_bar = zeros(T_p, M);     
sigma_bar  = zeros(T_p, M);     
for i = 1 : M
    active_win                = window_lengths_2_active_windows(n_windows, win_len);
    rho_level_win_regional    = rho_level_win(    (1 : n_windows) + (i-1)*n_windows );
    varphi_level_win_regional = varphi_level_win( (1 : n_windows) + (i-1)*n_windows );
    sigma_level_win_regional  = sigma_level_win(  (1 : n_windows) + (i-1)*n_windows );
    for k = 1 : T_p
        rho_bar(k, i)    = rho_pool(    rho_level_win_regional(    active_win(k) ) );
        varphi_bar(k, i) = varphi_pool( varphi_level_win_regional( active_win(k) ) );
        sigma_bar(k, i)  = sigma_pool(  sigma_level_win_regional(  active_win(k) ) );
    end
end

end


function active_win = window_lengths_2_active_windows(n_windows, win_len)

windows_begin = zeros(n_windows, 1);    % time instants where windows start
windows_end   = zeros(n_windows, 1);    % time instants where windows end
active_win    = zeros(sum(win_len), 1); % active window for each time instant. NB: sum(win_len) >= T_p

windows_begin(1) = 1;
windows_end(1)   = win_len(1);
active_win(windows_begin(1) : windows_end(1) ) = 1;
for i = 2 : n_windows
    windows_begin(i) = windows_end(i-1) + 1;
    windows_end(i)   = windows_begin(i) + win_len(i) - 1;
    active_win(windows_begin(i) : windows_end(i) ) = i;
end

end