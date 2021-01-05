clear,clc, close all
N = 36;
IFFT_N=512;
S=100000; % Number of samples
dB_scale=0:0.001:13-0.001;
%%OFDM(CCDF)
for k=1:S
real_value = randi([0,1], [1, N])*2 - 1;
imag_value = randi([0,1], [1, N])*2 - 1;
input_value = (real_value + imag_value*j)/sqrt(2);
time_samples = ifft(input_value,IFFT_N);
inst_power = abs(time_samples).^2;
peak_power = max(inst_power);
avg_power = mean(inst_power);
papr = peak_power/avg_power;
papr_dB = 10*log10(papr);

CCDF_OFDM_dB(k)=papr_dB;% Array(value of OFDM Papr_dB)

end
for k=1:length(dB_scale)
    M(k)=length(find(CCDF_OFDM_dB>dB_scale(k)));% Array(number of value that bigger than dB_scale(k))
    CCDF_OFDM_result(k)=M(k)/S; % result value
    
end

%%CCDF_SCFDMA
for k=1:S
real_value = randi([0,1], [1, N])*2 - 1;
imag_value = randi([0,1], [1, N])*2 - 1;
input_value = (real_value + imag_value*j)/sqrt(2);
dft_result=fft(input_value);
time_samples = ifft(dft_result,IFFT_N);
inst_power = abs(time_samples).^2;
peak_power = max(inst_power);
avg_power = mean(inst_power);
papr = peak_power/avg_power;
papr_dB = 10*log10(papr);

CCDF_scfdma_dB(k)=papr_dB;% Array(value of SC_FDMA Papr_dB)

end
for k=1:length(dB_scale)
    M(k)=length(find(CCDF_scfdma_dB>dB_scale(k)));% Array(number of value that bigger than dB_scale(k))
    CCDF_scfdma_result(k)=M(k)/S; % result value
    
end



semilogy(dB_scale,CCDF_OFDM_result)
xlabel('PAPR0 [dB]'), ylabel('CCDF'), axis([2 13 10^(-4) inf])
hold on
semilogy(dB_scale,CCDF_scfdma_result)
xlabel('PAPR0 [dB]'), ylabel('CCDF'), axis([2 13 10^(-4) inf])
legend('OFDM','SC-FDMA')
