clear,clc, close all
N = 36;
IFFT_N=512;

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

plot([1:IFFT_N],abs(inst_power))
avg_line=refline([0 avg_power]);
avg_line.Color='r';
peak_line=refline([0 peak_power]);
peak_line.Color='k';
axis([0 IFFT_N 0 (peak_power*1.1)])
ylabel('Instantaneous power(W)'),xlabel('Time(t)')
legend('inst power','avg power','peak power')