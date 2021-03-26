function controlled_Phi = compute_controlled_Phi(Phi, varphi)

M        = size(Phi, 1);
idx_diag = 1 : M+1 : M^2;

varphi_mat = varphi * varphi';
controlled_Phi = varphi_mat .* Phi;
controlled_Phi(idx_diag) = zeros(M, 1);
controlled_Phi(idx_diag) = ones(M, 1) - sum(controlled_Phi, 2);

end