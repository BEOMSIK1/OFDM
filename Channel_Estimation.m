clc, clear, close all;
%% Parameters
FFT_Size = 1024;
GI_Size = FFT_Size/4;
SNR = [0 10 30];
Multipath = 7;

%% Channel genrate
h = rayleigh_channel(Multipath);
H = transpose(fft(h,FFT_Size));
amp_H = abs(H);
ang_H = wrapTo2Pi(angle(H));

for SNR_index=1:length(SNR)
    %% Pilot signal
    len= FFT_Size;
    if mod(len,2)==0
        N=len+1;
    else
        N=len;
    end
    tmp = zadoffChuSeq(1,N);
    pilot = transpose(tmp(1:len));
    X = diag(pilot);
    %% Tx-Rx
    [Y,N_0] = awgn_noise(X*H,SNR(SNR_index));
    %% Least Square
    H_hat_LS(:,SNR_index) = inv(X)*Y;
    amp_LS = abs(H_hat_LS);
    ang_LS = wrapTo2Pi(angle(H_hat_LS));
    %% MMSE
    H_hat_MMSE(:,SNR_index) =(H*H')*inv((H*H')+(N_0*eye(FFT_Size)))*inv(X)*Y;
    amp_MMSE = abs(H_hat_MMSE);
    ang_MMSE = wrapTo2Pi(angle(H_hat_MMSE));
    
end

plot(1:FFT_Size,amp_H,1:FFT_Size,amp_LS)
title('Channel Estimation (LS)'), xlabel('Subcarrier'),ylabel('Amplitude'),legend('perfect Channel','SNR=0','SNR=10','SNR=30'),axis([0 FFT_Size -inf inf])
grid on

figure
plot(1:FFT_Size,amp_H,1:FFT_Size,amp_MMSE)
title('Channel Estimation (MMSE)'), xlabel('Subcarrier'),ylabel('Amplitude'),legend('perfect Channel','SNR=0','SNR=10','SNR=30'), axis([0 FFT_Size -inf inf])
grid on

figure
plot(1:FFT_Size,ang_H,1:FFT_Size,ang_LS)
title('Channel Estimation (LS)'), xlabel('Subcarrier'),ylabel('Phase'),legend('perfect Channel','SNR=0','SNR=10','SNR=30'), axis([0 FFT_Size -inf inf])
grid on

figure
plot(1:FFT_Size,ang_H,1:FFT_Size,ang_MMSE)
title('Channel Estimation (MMSE)'), xlabel('Subcarrier'),ylabel('Phase'),legend('perfect Channel','SNR=0','SNR=10','SNR=30'), axis([0 FFT_Size -inf inf])
grid on
