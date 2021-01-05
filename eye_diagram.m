clear,clc, close all

%%mode setting
Mode=[9 10]; %Line coding(user select)
%1:Unipolar NRZ(Non) 2:Polar NRZ(Non)   3:Unipolar RZ(Non)
%4:Polar RZ(Non)     5:Bipolar RZ(Non)  6:Manchester(Non)
%7:Unipolar NRZ(LPF) 8:Polar NRZ(LPF)   9:Unipolar RZ(LPF)
%10:Polar RZ(LPF)    11:Bipolar RZ(LPF) 12:Manchester(LPF)

Legend_All={'Unipolar NRZ(Non)','Polar NRZ(Non)','Unipolar RZ(Non)'...
    'Polar RZ(Non)','Bipolar RZ(Non)','Manchester(Non)','Unipolar NRZ(LPF)'...
    'Polar NRZ(LPF)','Unipolar RZ(LPF)','Polar RZ(LPF)','Bipolar RZ(LPF)'...
    'Manchester(LPF)'};
Legend_Select=Legend_All(Mode);

%%Parameter setting
Rb=10000; %bitrate
Tb=1/Rb; %bit time
Fs=50000; %sampling frequency
Ts=1/Fs; %sampling period

N=Tb/Ts; %sample per bit
L=100; %number of bits
Data=randi([0 1],[1 L]); %bitstream
N_All=L*N; %totla number of sample
Index=1:N:N_All; %impulse index

t=0:Ts:L*Tb-Ts; %time

%%pulse parameter
ONE      = ones(1,N);                               %NRZ
ONE_ZERO = [ones(1,ceil(N/2)) zeros(1,floor(N/2))]; %RZ
ONE_ONE  = [ones(1,ceil(N/2)) -ones(1,floor(N/2))]; %Manchester

%%pulse, nosise parameter
f_Cut=1000; %LPF cutoff frequency
Noise_Power=0.01; %power of noise

%%eye diagram parameter
Tb_Eye=2; %number of bits presented eye diagram
N_Eye=Tb_Eye*N; %number of samples presented eye diagram
t_Eye=t(1:N_Eye+1);%time period presented eye diagram

Iteration=50; %number of repeat

%%result graph label setting
t_Ticks=(Tb)*[0:Tb_Eye];
t_Ticks_Label={'0' 'T_b'};
for k=2:Tb_Eye
    t_Ticks_Label=[t_Ticks_Label [num2str(k) 'T_b']];
end

%%main
m=1; %figure index
for M=Mode
    R=rem(M,6);
    s=zeros(1,N_All);
    
    if R==1 %Unipolar NRZ
        Pulse=Data;
        s(Index)=Pulse;
        s=conv(s,ONE);
    elseif R==2 %Polar NRZ
        Pulse=2*Data-1;
        s(Index)=Pulse;
        s=conv(s,ONE);
    elseif R==3 %Unipolar RZ
        Pulse=Data;
        s(Index)=Pulse;
        s=conv(s,ONE_ZERO);
    elseif R==4 %Polar RZ
        Pulse=2*Data-1;
        s(Index)=Pulse;
        s=conv(s,ONE_ZERO);
    elseif R==5 %Bipolar
        Pulse=Data;
        U_=cumsum([1 Pulse]);
        U=(-1).^U_(2:end);
        
        s(Index)=U.*Pulse; %(-1)^n 이전극성표현
        s=conv(s,ONE_ZERO);
    elseif R==0 %Machester
        Pulse=2*Data-1;
        s(Index)=Pulse;
        s=conv(s,ONE_ONE);
    end
    
    
    s=s(1:N_All+1);
    F=ceil(M/6); %filter mode : F=1(Non), F=2 (LPF)
    s_Filter=(F==1)*s+(F==2)*lowpass(s,f_Cut,Fs,'Steepness',0.5);
    
    for k=1:Iteration
        
        y=s_Filter+sqrt(Noise_Power)*randn(1,N*L+1);
        Index_Eye=N_Eye*(k-1)+1:N_Eye*k+1;
        y_Eye=y(Index_Eye);
        
        figure(m)
        if F==1&&Noise_Power<10^(-6)
            hold on, stairs(t_eye,y_Eye,'-b','LineWidth',2)
        else
            hold on, plot(t_Eye,y_Eye,'-b','LineWidth',2)
        end
        
        xlabel('t [sec]'), ylabel('y(t)'), grid on
        legend(Legend_Select(m)),set(gca,'Fontsize',15),xticks(t_Ticks)
        xticklabels(t_Ticks_Label),ylim([-2 2])
    end
    m=m+1;
end

