function [Y,rk] = fn(m,n,k)
  % 创建秩为 k 的 m * n 矩阵
    P = orth(randn(m,k));
    Q = orth(randn(n,k))';
    Y = P*Q;
    rk = rank(Y);
end