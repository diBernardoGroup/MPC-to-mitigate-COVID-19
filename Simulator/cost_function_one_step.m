function J_step = cost_function_one_step( x, control, MODEL_PARAM )

% INPUT EXTRACTION --------------------------------------------------------
c_alpha         = MODEL_PARAM.c_alpha;
c_m             = MODEL_PARAM.c_m;
c_w             = MODEL_PARAM.c_w;
ctrl_idx_rho    = MODEL_PARAM.ctrl_idx_rho;
ctrl_idx_varphi = MODEL_PARAM.ctrl_idx_varphi;
ctrl_idx_sigma  = MODEL_PARAM.ctrl_idx_sigma;
Phi             = MODEL_PARAM.Phi;
M               = MODEL_PARAM.M;
N               = MODEL_PARAM.N;

% COMPUTE CONTROL ---------------------------------------------------------
rho    = control(ctrl_idx_rho)';  
varphi = control(ctrl_idx_varphi)';  
sigma  = control(ctrl_idx_sigma)';  

% UNSTACK x_bar -----------------------------------------------------------
S = x(        1:M  )';  % must be a column
I = x( M   + (1:M) )';
Q = x( 2*M + (1:M) )';
H = x( 3*M + (1:M) )';
D = x( 4*M + (1:M) )';  % unused
R = x( 5*M + (1:M) )';

% COMPTUTE COST -----------------------------------------------------------
Phi_orig_diag = diag(diag(Phi));

Phi_contr = compute_controlled_Phi(Phi, varphi);
idx_diag = 1 : M+1 : M^2;
Phi_contr(idx_diag) = zeros(M,1);

Phi_blocked = Phi - Phi_contr;
Phi_blocked(idx_diag) = zeros(M,1);

workers = S + I + R;
non_workers = c_w*Q + H + D;

J_step_regional = c_w * c_m * ( ...
    (1 - sqrt(rho)) .* ((Phi_orig_diag + Phi_contr)' * workers) ... % people who work in region i and suffer from rho(i)
    + sum(Phi_blocked, 2) .* workers ) ... % people who live in region i, but cannot leave towards other regions
    + c_m * non_workers ...                % people who live in region i but cannot work
    + c_alpha * N .* sigma;                % cost of increased testing

J_step = sum(J_step_regional);

end
