function discounts = compute_discounts( discount_factor, T_p )

if discount_factor > 0  &&  discount_factor <= 1
    A = discount_factor * ones(T_p);    
    A = triu(A) - diag(diag(A));
    B = ones(T_p);
    B = tril(B);
    discount_mat = A + B;
    discounts = prod(discount_mat, 1)';
else
    discounts = ones(T_p, 1);
end

end