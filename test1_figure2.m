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

T = 100*P;
M_theta = [1 0 0 0;0 1 1 -2]';
N_theta = [1 0];
%% 定义预测区间K
N=3;
% 定义step数量k
k_steps=90;
% 定义矩阵 X_K， n x k 矩 阵
X_K = zeros(n,k_steps);
% 定义记录k步的人工稳态矩阵 X_s， n*k 矩阵
X_s =zeros(n,k_steps);
%% 初始状态变量值， n x 1 向量
X_K(:,1) =[-5;0];
X_s(:,1) =[-5;0];
% 定义记录k步的输入矩阵 U_K， m x k 矩阵
U_K=zeros(m,k_steps);
% % 定义输出矩阵 hat_Y_s， n*k 矩阵
% hat_Y_s =zeros(n,k_steps);
% 定义记录k步的人工稳态输出矩阵 Y_s， n*k 矩阵
Y_s =zeros(n,k_steps);

% 定义记录k步的人工稳态输入矩阵 U_s， n*k 矩阵
U_s =zeros(m,k_steps);
% 定义k步的实际输出矩阵 Y_k， n*k 矩阵
Y_k =zeros(n,k_steps);
%% 创建参考信号矩阵
% 使用 repmat 函数重复矩阵  
matrix1 = repmat([4.95; 0], 1, 30); % 2x30 矩阵，重复 [4.95; 0] 30 次  
matrix2 = repmat([-5.5; 0], 1, 30); % 2x30 矩阵，重复 [-5.5; 0] 30 次  
matrix3 = repmat([2; 0], 1, 30);    % 2x30 矩阵，重复 [2; 0] 30 次    
% 水平拼接子矩阵  
hat_X_s = [matrix1, matrix2, matrix3]; % 2x90 矩阵

%% 计算每一步的状态变量的值
for k = 1 : k_steps
    % 求得U_K(:,k)
    % Call MPC_Matrices 函数 求得 E,H矩阵 
    [E_1,E_2,H,yue_zuo,yue_you,M_x_bar1,M_u_bar1,Kmpc]=MPC_Matrices(A,B,Q,R,P,T,N,M_theta,X_K(:,k));
    [u_k,theta] = Prediction(X_K(:,k),E_1,E_2,H,N,m,hat_X_s(:,k),yue_zuo,yue_you);
    U_K(:,k) = u_k; % 已经是优化后得到的第一个元素，排列成k_steps列
    
    % 计算第k+1步时状态变量的值
    % observer update
    % X_s(:,k+1)=A*hat_X_s(:,k)+B*U_K(:,k);
    X_s(:,k) = M_x_bar1*theta;

    X_K(:,k+1)=(A*X_K(:,k)+B*U_K(:,k));
    Y_k(:,k+1)=C*X_K(:,k);
    
    U_s (:,k) = M_u_bar1*theta;
    Y_s (:,k) = C*X_s(:,k);
    
    %定义不变集
    
end
 
%% 绘制状态变量和输入的变化

% 创建一个 2x1 的子图
figure;

% 上半部分
subplot(2, 1, 1); % 2 行 1 列的子图，选择第 1 个
hold on; % 保持当前图形
plot(hat_X_s(1,:), 'k-.', 'DisplayName', 'the output target'); % 输出目标
plot(Y_s(1,:), 'r--', 'DisplayName', 'artificial reference'); % 人工参考
plot(Y_k(1,:), 'b-', 'DisplayName', 'Evolution of the output'); % 输出的演变
xlabel('列索引');
ylabel('值');
title('上半部分：输出相关曲线');
legend; % 显示图例
grid on; % 添加网格
hold off; % 结束绘制

% 下半部分
% subplot(2, 1, 2); % 选择第 2 个子图
%% 绘制状态肖像
% 绘制状态肖像
% figure;
% hold on; % 保持当前图形
% plot(X_K(1,:), X_K(2,:), '-o', 'MarkerSize', 3, 'LineWidth', 1.5); % 绘制状态肖像
% xlabel('State Variable x_1');  % x_1 代表第一个状态变量（如位置）
% ylabel('State Variable x_2');  % x_2 代表第二个状态变量（如速度）
% title('State Portrait');
% grid on;
% axis equal; % 确保x和y轴的比例相同
% 
% % 添加箭头以表示状态变化
% for i = 1:length(X_K)-1
%     % 在每两个状态之间添加箭头
%     quiver(X_K(1,i), X_K(2,i), X_K(1,i+1)-X_K(1,i), X_K(2,i+1)-X_K(2,i), ...
%            'AutoScale', 'off', 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 1);
% end
% 
% hold off; % 结束绘制

%% 绘制输入输出
% figure;
% hold on; % 保持当前图形
% plot(U_K(1,:), 'm-', 'DisplayName', '实际控制输入'); % 实际控制输入
% plot(U_s(1,:), 'g-', 'DisplayName', '人工稳态输入'); % 人工稳态输入
% xlabel('列索引');
% ylabel('值');
% title('下半部分：控制输入曲线');
% legend; % 显示图例
% grid on; % 添加网格
% hold off; % 结束绘制
%% 状态变化
% figure;
% hold on;
% plot(X_K(1,:)); 
% plot(X_K(2,:));



