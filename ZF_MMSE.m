clc, clear, close all;
%% Parameters

Modulation_Order=2;                   % 변조 방법
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Nr=4;
Nt=4;
SNR=0:3:30;
Iteration=100;
%%
for SNR_index=1:length(SNR)
    N_0=10^(-SNR(SNR_index)/10);
    for Iter=1:Iteration
        %% MIMO-OFDM
        Data = randi([0 1],[Nt Data_Size]);
        mod_data=base_mod(Data,Modulation_Order)./sqrt(Nt);                           % 변조 방식에 따른 데이터 변조
        IFFT_data=ifft(mod_data,FFT_Size,2)*sqrt(FFT_Size);                 % 변조된 데이터 IFFT연산(각 반송파에서 보내는 신호의 power를 1로 하기위해 sqrt(반송파 개수)를 곱해줌
        Add_CP_data=[IFFT_data(:,FFT_Size-GI_Size+1:end), IFFT_data];       % IFFT연산된 데이터에 CP삽입
        
        copied_cpdata=repmat(Add_CP_data,Nr,1);
        
        h=rayleigh_channel([Multi_path, Nt*Nr]);                             % time domain channel
        H=fft(h,FFT_Size,2);                                                % frequency domain channel
        for k=1:(Nt*Nr)
            hx(k,:)=conv(copied_cpdata(k,:),h(k,:));
        end
        
        for k=1:Nr
            hx_comb(k,:)=sum(hx([(Nt*(k-1)+1):Nt*k],:));
        end
        
        y=awgn_noise(hx_comb,SNR(SNR_index));                                    % 채널 통과된 데이터에 awgn추가 (y=h*x+n)
        y_remove_CP=y(:,GI_Size+1:GI_Size+FFT_Size);                        %remove CP
        Y=fft(y_remove_CP,FFT_Size,2)/sqrt(FFT_Size)*sqrt(Nt);                       % CP제거된 데이터 FFT연산(신호의 power를 1로 하기 위해 sqrt(반송파 개수)로 나눠줌)
        
        %% Zero Forcing
        H_rv=reshape(H,Nt,Nr,[]);
        for K=1:FFT_Size
            H_ch=transpose(H_rv(:,:,K));
            G_zf=inv(H_ch'*H_ch)*H_ch';
            X_hat(:,K)=G_zf*Y(:,K);
        end
        X_hat_demod_zf=base_demod(X_hat,Modulation_Order);
        num_error_zf(Iter,SNR_index)=biterr(Data,X_hat_demod_zf);
        %% Minimum Mean-Squared Error
        for K=1:FFT_Size
            H_ch=transpose(H_rv(:,:,K));
            H_bar=[H_ch;sqrt(N_0)*eye(Nt)];
            Y_bar=[Y;zeros(Nt,FFT_Size)];
            G_mmse=inv(H_bar'*H_bar)*H_bar';
            X_hat_mmse(:,K)=G_mmse*Y_bar(:,K);
        end
        X_hat_demod_mmse=base_demod(X_hat_mmse,Modulation_Order);
        num_error_mmse(Iter,SNR_index)=biterr(Data,X_hat_demod_mmse);
        
    end
    error_rate_zf=(sum(num_error_zf,1)/(Data_Size*Nt))/Iteration;
    error_rate_mmse=(sum(num_error_mmse,1)/(Data_Size*Nt))/Iteration;
end
semilogy(SNR,error_rate_zf,'-o')
hold on
semilogy(SNR,error_rate_mmse,'-*')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('ZF','MMSE'),axis([0 30 1e-4 1]),grid on