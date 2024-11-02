%% 清屏
clear ;
close all;
clc;
% 第一步，定义状态空间矩阵
% 定义状态矩阵 A, n x n 矩阵
A = [1 1; 0 1];
n= size (A,1);
% 定义输入矩阵 B, n x m 矩阵
B = [0 0.5; 1 0.5];
m = size(B,2);
C= [1 0];
p = 1;
%% 定义Q矩阵,n x n 矩阵
Q=[1 0;0 1];
% 定义R矩阵，p x p 矩阵
R=[1 0;0 1];
% 定义P矩阵，n x n 矩阵,Riccati equation solution
% K = [-0.0173 0.5367;0.5541 0.8225];
% Ak =A-B*K;
% Qk =Q+K'*R*K;
% P = dlyap(Ak',Qk);
P = dare(A, B, Q, R);
[Kmpc, ~, ~] = lqr(A, B, Q, R);
T = 100*P;
M_theta = [1 0 0 0;0 1 1 -2]';
N_theta = [1 0];

[H_rpi, h_rpi] = compute_rpi_git(A,B,Kmpc)

%% 绘图
figure;
hold on;

% 使用线段绘制约束区域的边界
x1 = linspace(-5, 5, 100);

for i = 1:size(H_rpi, 1)
    % 计算约束线
    if H_rpi(i, 2) ~= 0
        y1 = (h_rpi(i) - H_rpi(i, 1) * x1) / H_rpi(i, 2);
        plot(x1, y1, 'k--', 'LineWidth', 1.5);
    end
end

% 设置坐标轴和图形标题
xlabel('x_1');
ylabel('x_2');
title('由 H_rpi 和 h_rpi 定义的约束边界');
axis equal;
xlim([-6 6]);
ylim([-6 6]);
grid on;
hold off;