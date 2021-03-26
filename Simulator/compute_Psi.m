function Psi = compute_Psi(x, control, MODEL_PARAM)
% input are rows
% using infinite norm

% INPUT EXTRACTION --------------------------------------------------------
alpha_0         = MODEL_PARAM.alpha_0;
alpha_tilde     = MODEL_PARAM.alpha_tilde;
beta            = MODEL_PARAM.beta;
ctrl_idx_rho    = MODEL_PARAM.ctrl_idx_rho;
ctrl_idx_varphi = MODEL_PARAM.ctrl_idx_varphi;
ctrl_idx_sigma  = MODEL_PARAM.ctrl_idx_sigma;
gamma           = MODEL_PARAM.gamma;
M               = MODEL_PARAM.M; 
Phi             = MODEL_PARAM.Phi; 
psi             = MODEL_PARAM.psi;

% COMPUTE CONTROL ---------------------------------------------------------
rho    = control(ctrl_idx_rho)';    % rho is now a column
varphi = control(ctrl_idx_varphi)'; % varphi is now a column
sigma  = control(ctrl_idx_sigma)';  % sigma is now a column

Phi   = compute_controlled_Phi(Phi, varphi);
alpha = alpha_0 + sigma .* alpha_tilde;

% UNSTACK STATE -----------------------------------------------------------
S = x(        1:M  )';   % must be a column
I = x( M   + (1:M) )';
R = x( 5*M + (1:M) )';

% COMPUTE NORM ------------------------------------------------------------
N_p = Phi' * (S + I + R);
Psi = zeros(M);
temp = eye(M).*(alpha + psi + gamma);

for i = 1:M
    for j = 1:M
        Psi(i,j) = S(i) * beta(i) * (Phi(i,:) .* rho' ./ N_p' * Phi(j,:)') - temp(i,j);
    end
end

end