function z = OFDM_compare (N_samples, ifft_size)



f0 = 1e3;

T0 = 1/f0;



ifft_input_values = [0, randi([-1, 1], [1, N_samples])];



N_analog_samples = 1000;

x_analog = linspace(0, T0, N_analog_samples);

y_analog_real = zeros(size(x_analog));

y_analog_imag = zeros(size(x_analog));



for i=2:N_samples+1 

  y_analog_real = y_analog_real + ifft_input_values(i) * cos(2*pi*(i-1)*f0*x_analog);   

  y_analog_imag = y_analog_imag + ifft_input_values(i) * sin(2*pi*(i-1)*f0*x_analog);

end



x_sampled_values = linspace(0, T0*(1-1/ifft_size), ifft_size);

y_time_samples = ifft(ifft_input_values, ifft_size);



subplot(2,1,1);

plot(x_analog, y_analog_real, x_sampled_values, real(y_time_samples)*ifft_size, 'o');

subplot(2,1,2);

plot(x_analog, y_analog_imag, x_sampled_values, imag(y_time_samples)*ifft_size, 'o');

