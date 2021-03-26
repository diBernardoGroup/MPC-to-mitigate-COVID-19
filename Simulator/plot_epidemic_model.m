function fig = plot_epidemic_model(x_bar, control_bar, tau_bar, R_t, MODEL_PARAM, CONTROL_PARAM, PREV_STATE)

addpath('utilities');
% INPUT EXTRACTION --------------------------------------------------------
ctrl_idx_rho    = MODEL_PARAM.ctrl_idx_rho;
ctrl_idx_varphi = MODEL_PARAM.ctrl_idx_varphi;
ctrl_idx_sigma  = MODEL_PARAM.ctrl_idx_sigma;
ICU_max         = MODEL_PARAM.ICU_max;
M               = MODEL_PARAM.M;

epsilon_H = CONTROL_PARAM.epsilon_H;
epsilon_R = CONTROL_PARAM.epsilon_R;
tau_min = CONTROL_PARAM.tau_min;
norm_bound = CONTROL_PARAM.norm_bound;

rho_0    = PREV_STATE.control_0(ctrl_idx_rho);
varphi_0 = PREV_STATE.control_0(ctrl_idx_varphi);
sigma_0  = PREV_STATE.control_0(ctrl_idx_sigma);
t_1      = PREV_STATE.t_1;

% INITIALIZATION ----------------------------------------------------------
rho_bar    = control_bar(:, ctrl_idx_rho);  
varphi_bar = control_bar(:, ctrl_idx_varphi);  
sigma_bar  = control_bar(:, ctrl_idx_sigma);  
H_to_ICU_scaling = 0.1;

fig = figure;    set(gcf, 'color', 'w')
set(fig, 'defaultAxesColorOrder', [[0 0 0]; [1 0 0]] )
plot_rows = 3;
T_x       = size(x_bar, 1);
T_c       = size(control_bar, 1);
T_max     = max(T_x, T_c);
time      = t_1 : t_1 + T_x - 1;
shadings_height = 1e10;
scatter_size = 15;

% COLORS ------------------------------------------------------------------
colors       = lines(7);
color_blue   = colors(1, :);
color_orange = colors(2, :);
color_gold   = colors(3, :);
% color_purple = colors(4, :);
% color_green  = colors(5, :);
% color_l_blue = colors(6, :);
% color_red    = colors(7, :);

% PLOT STATES -------------------------------------------------------------
shadings_states = zeros(time(end-1), M);
for k = time(1 : end-1)
    Psi = compute_Psi(x_bar(k,:), control_bar(k, :), MODEL_PARAM);
    dynamics_I = eye(M) + Psi;
    norm_rows = sum( abs(dynamics_I), 2);
    shadings_states(k, :) = norm_rows' < norm_bound;
    prova(k, :) = norm_rows';
end

for i = 1 : M
    subplot(plot_rows, M, i);    hold on; grid on;
    yyaxis right
    
    plot_shadings(shadings_states(:, i), time(1:end-1), shadings_height, 'red', 'green');
    line_ICU_max   = yline(ICU_max(i), 'r--');    
    line_ICU_thres = yline(epsilon_H * ICU_max(i), 'r-.');  
    H_i = x_bar(:, i + 3*M);
    line_ICU = plot(time, H_to_ICU_scaling * H_i,'r-','linewidth', 1);
    ylim([ 0, 1.05 * max(ICU_max(i), H_to_ICU_scaling*max(H_i)) ]);
    
    yyaxis left
    I_i = x_bar(:, i + 1*M);
    Q_i = x_bar(:, i + 2*M);
    line_I = plot(time, I_i, '-' , 'Color', color_blue, 'linewidth', 1);
    line_Q = plot(time, Q_i, '-' , 'Color', 'black', 'linewidth', 1);
    
    xlim([t_1 t_1 + T_max]);
    
%     if i == 1
%         ylabel('states', 'Interpreter', 'latex');
%     end
%     if i == M
%         legend([line_I, line_Q, line_ICU, line_ICU_max, line_ICU_thres], ...
%                 'I', 'Q', 'ICU', 'ICU_{max}', 'ICU thres.'  );
%     end
    
    box on
    set(gca,'FontSize',13.5);
end

subplot(plot_rows, M, 1);  title('North');
subplot(plot_rows, M, 2);  title('Center');
subplot(plot_rows, M, 3);  title('South');

% PLOT CONTROL ------------------------------------------------------------
shadings_control = tau_bar < (tau_min .* ones(time(end-1), M));

for i = 1 : M
    subplot(plot_rows, M, M+i);    hold on; grid on;
    plot_shadings(shadings_control(:, i), time(1:end-1), shadings_height, 'black', 0);   
    if sigma_0(i) ~= -1,  scatter(t_1, sigma_0(i),  scatter_size, color_gold,   'filled' ); end
    if varphi_0(i) ~= -1, scatter(t_1, varphi_0(i), scatter_size, color_orange, 'filled');  end
    if rho_0(i) ~= -1,    scatter(t_1, rho_0(i),    scatter_size, color_blue,   'filled');  end
%   plot(t_1:+t_1+T_c-1, tau_bar/tau_min,        'color', [color_green, 0.3], 'linewidth', 2.5);
    line_sigma  = plot(t_1:+t_1+T_c-1, sigma_bar(:, i),        'color', color_gold,   'linewidth', 1);    
    line_varphi = plot(t_1:+t_1+T_c-1, varphi_bar(:, i), '-.', 'color', color_orange, 'linewidth', 1);    
    line_rho    = plot(t_1:+t_1+T_c-1, rho_bar(:, i),    '--', 'color', color_blue,   'linewidth', 1);      
    
    xlim([t_1 t_1 + T_max]);
    ylim([0 1.05]);    
%     if i == 1
%         ylabel('control', 'Interpreter', 'latex');
%     end
%     if i == M
%         legend([line_rho, line_varphi, line_sigma], ...
%             '\rho_i', 'varphi_i', '\sigma_i');
%         %set(my_legend, 'Interpreter', 'latex');
%     end
    box on
    set(gca,'FontSize',13.5);
end

% PLOT R_t ----------------------------------------------------------------
max_R = max( max(max(R_t)), epsilon_R);
for i = 1 : M
    subplot(plot_rows, M, 2*M+i);    hold on;    grid on;
%     yyaxis left
%     line_S = plot(time, x_bar(:, i), 'linewidth', 1);
%     if i == 1
%         ylabel('repr. num.', 'Interpreter', 'latex');
%     end
%     axes = gca;
%     axes.XGrid = 'on';
%     axes.YGrid = 'off';
%     
%     yyaxis right
    idx_plot = find(R_t(:, i) ~= -1);
    line_R_t_thres = yline(epsilon_R, 'r-.');
    line_R_t       = plot(t_1 + idx_plot-1, R_t(idx_plot, i), 'linewidth', 1);
    ylim([0 max(1, 1.05*max_R)]);
    xlim([t_1 t_1 + T_max]);
%     if i == M
%         legend([line_S, line_R_t, line_R_t_thres], 'S', 'R_t','R_t thres.');
%     end
%     xlabel('$t$', 'Interpreter', 'latex');
    
    box on
    set(gca,'FontSize',13.5);
end
set( gcf, 'Position', [100 100 900 600])
saveas(gcf,'.\fig1.pdf')
end