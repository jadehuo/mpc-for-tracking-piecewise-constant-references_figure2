function [u_k,theta]= Prediction(x_k,E_1,E_2,H,N,m,hat_x_s,yue_zuo,yue_you)
    n =length(x_k); 
    f = x_k'*E_1-hat_x_s'*E_2;% 线性项
    S = quadprog(H,f,yue_zuo,yue_you);
    % U_k = O_1*S;% Nm x 1,这个U_k和外面的U_k不一样
    u_k = S(1:m,1);
    theta =S(end-1:end,1);
end