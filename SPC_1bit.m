function [precoded_X,H_spc] = SPC_1bit(FFT_Size,mod_data,H,angle_diff)


angle_idx_1bit=angle_diff>(pi/2);               % 조건만족 -> state2, 불만족 -> state1
state1_idx=find(angle_idx_1bit==0);
state2_idx=find(angle_idx_1bit==1);
pc_vector_1bit(1,1:FFT_Size)=1/sqrt(2);
pc_vector_1bit(2,state1_idx)=1/sqrt(2);
pc_vector_1bit(2,state2_idx)=exp(-j*pi)/sqrt(2);
precoded_X=mod_data(1,:).*pc_vector_1bit;
H_spc=sum(H.*pc_vector_1bit);
