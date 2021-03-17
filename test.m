clc, clear, close all;
%% Parameters
Modulation_Order=2;                   % 변조 방법
FFT_Size=128;                         % 반송파 갯수
Data_Size =FFT_Size*Modulation_Order; % 데이터 크기
GI_Size=FFT_Size/4;                   % CP size
Multi_path=7;                         % 다중경로 갯수
Nr=2;
Nt=2;
SNR=0:3:30;
Iteration=1000;
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
        %% Zero Forcing-OSIC
        H_rv=reshape(H,Nt,Nr,[]);
        X_hat_zf_osic=ZF_OSIC(FFT_Size,Modulation_Order,Nt,H_rv,Y);
        X_hat_demod_zf_osic=base_demod(X_hat_zf_osic,Modulation_Order);
        num_error_zf_osic(Iter,SNR_index)=biterr(Data,X_hat_demod_zf_osic);
        %% Minimum Mean-Squared Error_OSIC
        X_hat_mmse_osic=MMSE_OSIC(FFT_Size,Modulation_Order,Nt,H_rv,Y,N_0);
        X_hat_demod_mmse_osic=base_demod(X_hat_mmse_osic,Modulation_Order);
        num_error_mmse_osic(Iter,SNR_index)=biterr(Data,X_hat_demod_mmse_osic);
        %% Zero Forcing_DFE
        X_hat_zf_dfe=ZF_DFE (FFT_Size,Modulation_Order,Nt,Nr,H_rv,Y);
        X_hat_demod_zf_dfe=base_demod(X_hat_zf_dfe,Modulation_Order);
        num_error_zf_dfe(Iter,SNR_index)=biterr(Data,X_hat_demod_zf_dfe);
        %% Minimum Mean-Squared Error_DFE
        X_hat_mmse_dfe=MMSE_DFE(FFT_Size,Modulation_Order,Nt,Nr,H_rv,Y,N_0);
        X_hat_demod_mmse_dfe=base_demod(X_hat_mmse_dfe,Modulation_Order);
        num_error_mmse_dfe(Iter,SNR_index)=biterr(Data,X_hat_demod_mmse_dfe);
    end
    error_rate_zf_dfe=(sum(num_error_zf_dfe,1)/(Data_Size*Nt))/Iteration;
    error_rate_mmse_dfe=(sum(num_error_mmse_dfe,1)/(Data_Size*Nt))/Iteration;
    error_rate_zf_osic=(sum(num_error_zf_osic,1)/(Data_Size*Nt))/Iteration;
    error_rate_mmse_osic=(sum(num_error_mmse_osic,1)/(Data_Size*Nt))/Iteration;
end
semilogy(SNR,error_rate_zf_osic,'-o')
hold on
semilogy(SNR,error_rate_mmse_osic,'-*')
hold on
semilogy(SNR,error_rate_zf_dfe,'-+')
hold on
semilogy(SNR,error_rate_mmse_dfe,'-v')

title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER')
legend('ZF-OSIC','MMSE-OSIC','ZF-DFE','MMSE-DFE'),grid on