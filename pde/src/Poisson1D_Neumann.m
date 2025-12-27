%% ------------------------------------------------------------------------
%% BÀI 2: ĐIỀU KIỆN NEUMANN
%% ------------------------------------------------------------------------
% Bài toán: -u'' = 1, x ∈ [0,1]
%          u'(0) = 0, u'(1) = 0
% Nghiệm giải tích: u(x) = -x²/2 + x/2 + C (chọn C sao cho ∫u dx = 0)
a=0; b=1;
n = 50;
h =(b-a)/(n+1);
x = linspace(a,b,n+2)';
du_a =0; du_b=0;
f = @(x) ones(size(x));
A = spdiags(ones(n+2,1)*[1 -2 1],-1:1, n+2,n+2)/h^2;
A(1,1) =-2/h^2;A(1,2) = 2/h^2;
A(end,end) = -2/h^2;A(end,end-1)=2/h^2;
F = f(x);
F(1) = F(1) - 2*du_a/h;
F(end) = F(end) +2*du_b/h;

u_num = A \ F;
u_num = u_num - mean(u_num);

u_exact = -x.^2/2 + 1/6;
u_exact = u_exact - mean(u_exact);
% So sánh nghiệm

figure;
plot(x, u_num, x, u_exact);
xlabel('x');
ylabel('u');
legend('Nghiệm số', 'Nghiệm giải tích');
grid on;
