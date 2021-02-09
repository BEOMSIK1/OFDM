function [Y_demod,S] = MU_BD(P,power_index,X_mod,H,SNR,modulation_order,Nt,Nr)
K=Nt/Nr;
%% Precoding Matrix
for i=1:K
    idx = Nr*i;
    H_tilde=H;
    H_tilde(idx-(Nr-1):idx,:) = [];                               % i 번째 채널을 제외한 H_slash_i 생성
    [~,~,V] = svd(H_tilde);                                       % singular value decomposition
    F(:,idx-(Nr-1):idx) = transpose(V(Nt-Nr+1:Nt,:));                        % precoding matrix (F_j)
end
%% Tx-Rx
X_power=1*(10^(P(power_index)/10));                               % rho 값 (modulation된 신호의 power, 프로베니우스 놈)
gamma=trace(F*F')./X_power;                                       % power 정규화 factor
precoded_X=(F*X_mod)./gamma;                                      % 사용자 신호 전처리
Y_received=awgn_noise(H*precoded_X,SNR);                          % H*전처리된 각 사용자 신호/gamma + noise
Y=Y_received.*gamma;                                              % 정규화 factor를 곱해주어 신호추출
H_BD=H*F;                                                         % BD matrix
Y_ZF=(inv(H_BD'*H_BD)*H_BD')*Y;                                   % ZF를 이용한 single user detection
Y_demod = base_demod(Y_ZF,modulation_order);
%% Sum rate
S=log2(det(eye(Nr*K)+X_power*(H_BD)*(H_BD)'));
