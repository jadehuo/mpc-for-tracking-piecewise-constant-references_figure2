function  [E_1,E_2,H,yue_zuo,yue_you,M_x_bar1,M_u_bar1,Kmpc]=MPC_Matrices(A,B,Q,R,P,T,N,M_theta,x_k)
    n=size(A,1);   % A 是 n x n 矩阵, 得到 n
    m=size(B,2);   % B 是 n x m 矩阵, 得到 m
    %%%%%%%%%%%%
    M=[eye(n);zeros(N*n,n)]; % 初始化 M 矩阵. M 矩阵是 (N+1)n x n的， 
                             % 它上面是 n x n 个 "I", 这一步先把下半部
                             % 分写成 0 
    C_bar=zeros((N+1)*n,N*m); % 初始化 C 矩阵, 这一步令它有 (N+1)n x Nm 个 0
    % 定义M 和 C 
    tmp=eye(n);  %定义一个n x n 的 I 矩阵，存放A和幂次方
    %　更新Ｍ和C
    for i=1:N % 循环，i 从 1到 N
        rows =i*n+(1:n); %定义当前行数，从i x n开始，共n行 。如果 i = 2 且 n = 3，则 rows = 2*3 + (1:3)，即 rows = 6 + [1, 2, 3]，结果是 rows = [7, 8, 9]。
        C_bar(rows,:)=[tmp*B,C_bar(rows-n, 1:end-m)]; %将c矩阵填满，从前一个时间步的C矩阵中取出对应的列，并附加到当前的C矩阵
        tmp= A*tmp; %每一次将tmp左乘一次A
        M(rows,:)=tmp; %将M矩阵写满
    end
    %定义M_x, M_u
    M_1 = eye(n);
    M_2 = zeros(n);
    M_x =[M_1,M_2];
    M_u = [M_2,M_1];
    % 定义Q_bar和R_bar
    Q_bar1 = kron(eye(N),Q);
    Q_bar = blkdiag(Q_bar1,P);
    R_bar = kron(eye(N),R);
    M_x_bar1 = M_x*M_theta;
    M_x_bar = repmat(M_x_bar1, N+1, 1);%M_x_bar: （N+1）n x n
    M_u_bar1 = M_u*M_theta;
    M_u_bar = repmat(M_u_bar1, N, 1);%M_u_bar: Np x n
    %定义O_1,O_2
    O_11 = eye(N*m);
    O_12 = zeros(N*m,n);
    O_1 = [O_11,O_12];%O_1: Nm x (Nm+n)
    O_21 = zeros(n,N*m);
    O_22 = eye(n);
    O_2 = [O_21,O_22];%O_1: n x (Nm+n)
    % 计算G, E, H

    % G=M'*Q_bar*M; % G: n x n常数项
    E_1 = M'*Q_bar*(C_bar*O_1-M_x_bar*O_2);
    E_2 = T*M_x_bar1*O_2; % E:  nx (Nm+n)
    
    H = O_1'*(C_bar'*Q_bar*C_bar+R_bar)*O_1+O_2'*(M_x_bar'*Q_bar*M_x_bar+M_u_bar'*R_bar*M_u_bar+M_x_bar1'*T*M_x_bar1)*O_2-2*O_1'*C_bar'*Q_bar*M_x_bar*O_2-2*O_2'*M_u_bar'*R_bar*O_1; % H:(Nm+n) x (Nm+n)
    %% 约束条件
    % 计算 K,L
   
    [Kmpc, ~, ~] = lqr(A, B, Q, R);
    
    % F = [eye(N*m),zeros(N*m,1),zeros(N*m,1)]*M;%6*2矩阵
    % phi = [eye(N*m),zeros(N*m,1),zeros(N*m,1)]*C_bar;%6*6矩阵
    
    % Kmpc = [eye(m),zeros(m),zeros(m)]*inv(C_bar'*Q_bar*C_bar+R_bar)*C_bar'*Q_bar*M;
    
    % Kmpc = inv(R+B'*P*B)*B'*P*A;
    L = [-Kmpc, eye(m)] * M_theta;
    
    A_z_bar = [eye(n);-eye(n)];
    %左边的约束,
    yue1 = blkdiag(A_z_bar,A_z_bar,A_z_bar)*O_1;
    % yue211 = [eye(N*n),zeros(N*n,n)];
    % yue2 = blkdiag(A_z_bar,A_z_bar,A_z_bar)*yue211*C_bar*O_1;
   
    yue2 = blkdiag(A_z_bar,A_z_bar,A_z_bar,A_z_bar)*C_bar*O_1;
    yue4 = A_z_bar*(Kmpc*O_2*C_bar*O_1+L*O_2);

   
        lamda = 0.99;
        yue3 = blkdiag(A_z_bar,A_z_bar)*M_theta*O_2;
        yue311 = lamda*5*ones(2*n,1);
        yue312 = lamda*0.3*ones(2*m,1);
        yue31 = [yue311;yue312];
    
    %右边的约束
    
    yue11 = 0.3*ones(2*m*N,1);
    % yue21 = 5*ones(2*n*N,1)-blkdiag(A_z_bar,A_z_bar,A_z_bar)*yue211*M*x_k;
    yue21 = 5*ones(2*n*(N+1),1)-blkdiag(A_z_bar,A_z_bar,A_z_bar,A_z_bar)*M*x_k;
    yue41 = 0.3*ones(2*m,1)-A_z_bar*Kmpc*O_2*M*x_k;
    % 左右的约束整合起来
    yue_zuo = [yue1;yue2;yue3;yue4];
    yue_you = [yue11;yue21;yue31;yue41];
   
end