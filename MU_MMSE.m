function [Y_demod] = MU_MMSE(P,power_index,X_mod,H,SNR,modulation_order,K,Nt)

%% Precoding (MMSE)
tx_power=1*(10^(P(power_index)/10));                                 % rho 
z=K/tx_power;
gamma=trace(inv(H*H'+z*eye(Nt)))./tx_power;                          % normalize
precoded_matrix = (H'*inv(H*H'+z*eye(Nt)))./sqrt(gamma);             % 전처리 행렬
%% Tx-Rx
precoded_X = precoded_matrix*X_mod;              % 사용자 신호 전처리
Y_received = awgn_noise(H*precoded_X,SNR);       % H*전처리된 각 사용자 신호/sqrt(gamma) + noise
Y = Y_received.*sqrt(gamma);                     % 정규화 factor를 곱해주어 사용자 신호만 추출
Y_demod = base_demod(Y,modulation_order);