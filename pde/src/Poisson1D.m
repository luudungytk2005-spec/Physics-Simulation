%% ========================================================================
%% SAMPLE: POISSON 1D - MATLAB
%% ========================================================================

%% ------------------------------------------------------------------------
%% Problem 1: Dirichlet boundary condition
%% ------------------------------------------------------------------------
% Problem: -u'' = sin(pi*x), x ∈ [0,1]
%          u(0) = 0, u(1) = 0
%  u(x) = sin(pi*x)/(pi^2)
a =0;b=1; % 2 biên
N = 50;h=(b-a)/(N+1); % điểm trong miền và bước lưới
x = linspace(a,b,N+2)'; %lưới đầy đủ
u_a =0; u_b =0; % điều kiện biên
u_exact = @(x) sin(pi*x)/(pi^2); % nghiệm giải tích
f = @(x) sin(pi*x); % hàm nguồn
e = ones(N, 1);
A = spdiags([e -2*e e], -1:1, N, N) / h^2;

M = zeros(N,N);
for i=1:N
    M(i,i) = -2;
    if i>1,M(i,i-1) = 1;end
    if i<N,M(i,i+1) = 1;end
end
M = M / h^2;
F = f(x(2:end-1));
F(1) = F(1) -u_a/h^2;
F(end) = F(end) - u_b/h^2;
u_inner = A \ F;
u_num = [u_a;u_inner;u_b];
u_ex = u_exact(x);
error = norm(u_num -u_ex,inf);
fprintf('Sai so cuc dai: %.6e\n',error)
figure('Name', 'Bài 1: Điều kiện Dirichlet');
plot(x, u_num, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
hold on;
plot(x, u_ex, 'r--', 'LineWidth', 2);
grid on;
xlabel('x'); ylabel('u(x)');
title('Phương trình Poisson với điều kiện Dirichlet');
legend('Nghiệm số', 'Nghiệm giải tích', 'Location', 'best');

fprintf('Bài 1 hoàn thành!\n\n');


