function [precoded_X,H_spc] = SPC_2bit(FFT_Size,mod_data,H,Nt)

%% Channel info.
ag=angle(H)-angle(H(1,:));
ag_idx=transpose(find(ag(2,:)<0));
ag_temp=zeros(Nt-1,FFT_Size);
ag_temp(ag_idx)=2*pi+ag(2,ag_idx);
ag(2,ag_idx)=0;
angle_H=ag(2,:)+ag_temp;
%% SPC 2-bit
angle_idx_1=find(angle_H<(pi/4)|angle_H>=(7*pi/4));             % state1 index
angle_idx_2=find(angle_H>=(pi/4)&angle_H<(3*pi/4));             % state2 index
angle_idx_3=find(angle_H>=(3*pi/4)&angle_H<(5*pi/4));           % state3 index
angle_idx_4=find(angle_H>=(5*pi/4)&angle_H<(7*pi/4));           % state4 index
pc_vector_2bit(1,1:FFT_Size)=1/sqrt(2);
pc_vector_2bit(2,angle_idx_1)=1/sqrt(2);
pc_vector_2bit(2,angle_idx_2)=exp(-j*pi/2)/sqrt(2);
pc_vector_2bit(2,angle_idx_3)=exp(-j*pi)/sqrt(2);
pc_vector_2bit(2,angle_idx_4)=exp(j*pi/2)/sqrt(2);
precoded_X=mod_data(1,:).*pc_vector_2bit;
H_spc=sum(H.*pc_vector_2bit);