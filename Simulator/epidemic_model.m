function x_k_p_1 = epidemic_model(x, control, MODEL_PARAM) 
% x and control are given as rows, x_k_p_1 must be a row

% INPUT EXTRACTION --------------------------------------------------------
alpha_0         = MODEL_PARAM.alpha_0;
alpha_tilde     = MODEL_PARAM.alpha_tilde;
beta            = MODEL_PARAM.beta;
ctrl_idx_rho    = MODEL_PARAM.ctrl_idx_rho;
ctrl_idx_varphi = MODEL_PARAM.ctrl_idx_varphi;
ctrl_idx_sigma  = MODEL_PARAM.ctrl_idx_sigma;
eta_H           = MODEL_PARAM.eta_H;
eta_Q           = MODEL_PARAM.eta_Q;
gamma           = MODEL_PARAM.gamma;
ICU_max         = MODEL_PARAM.ICU_max;
kappa_H         = MODEL_PARAM.kappa_H;
kappa_Q         = MODEL_PARAM.kappa_Q;
M               = MODEL_PARAM.M; 
Phi             = MODEL_PARAM.Phi; 
psi             = MODEL_PARAM.psi;
zeta_baseline   = MODEL_PARAM.zeta_baseline;
zeta_slope      = MODEL_PARAM.zeta_slope;

% COMPUTE CONTROL ---------------------------------------------------------
rho    = control(ctrl_idx_rho)';    % rho is now a column
varphi = control(ctrl_idx_varphi)'; % varphi is now a column
sigma  = control(ctrl_idx_sigma)';  % sigma is now a column

Phi   = compute_controlled_Phi(Phi, varphi);
alpha = alpha_0 + sigma .* alpha_tilde;

% UNSTACK x_bar -----------------------------------------------------------
S = x(        1:M  )';  % must be a column
I = x( M   + (1:M) )';
Q = x( 2*M + (1:M) )';
H = x( 3*M + (1:M) )';
%D = x( 4*M + (1:M) )';  % unused
R = x( 5*M + (1:M) )';
ICU = 0.1 * H;

% COMPUTE NEW VALUES ------------------------------------------------------
N_p          = Phi' * (S + I + R);
new_infected = S .* Phi * (rho .* beta ./ N_p .* (Phi' * I));
zeta         = zeta_baseline + zeta_slope * min(ICU ./ ICU_max', 1);

x_k_p_1 = x + [ - new_infected; ...
                new_infected - (gamma + alpha + psi) .* I; ...
                alpha .* I - (kappa_H + eta_Q) .* Q + kappa_Q .* H; ...
                psi .* I + kappa_H .* Q - (eta_H + zeta + kappa_Q) .* H; ...
                zeta .* H; ...
                gamma .* I + eta_Q .* Q + eta_H .* H]';  % output must be a row

end