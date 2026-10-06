function [B, Tc, seB, seTc] = ols_line(x, y, direction)
    % tau = B*qd + Tc*direction 的最小平方解與標準誤（殘差 iid 假設）
    x = x(:); y = y(:);
    X = [x, direction * ones(size(x))];
    c = X \ y;
    B = c(1); Tc = c(2);
    r = y - X * c;
    dof = max(numel(y) - 2, 1);
    se = sqrt(diag((r' * r) / dof * inv(X' * X)));
    seB = se(1); seTc = se(2);
end
