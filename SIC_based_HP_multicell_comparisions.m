% Multi-cell case - comparision with different parameters 
% Test_2 (Extra for rough work)
% -----------------------------------------------------------------------

% (1) comparision with Nr(number of Rx antennas)
clear
close all
%clc
tic
%SNR_dB=[-30:5:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
Nr_RF=4; % number of RF chains in each user(u) at Rx=no. of data streams at each u

Nt=128; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);

M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side

nr=64; % number of Rx antennas 

L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 
Nr_plot=[8:8:nr];

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_Nr=1:length(Nr_plot)
    i_Nr
    
    Nr=Nr_plot(i_Nr); % parameter (1) 
    Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
    ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
    ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
    Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
    
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_Nr)= real(temp4/N_iter);
    C3(i_Nr)= real(temp3/N_iter);     
    C1(i_Nr)= real(temp1/N_iter);
    C2(i_Nr)= real(temp2/N_iter);
end

plot(Nr_plot,C4,'r','Linewidth',1.5);
hold on
plot(Nr_plot,C3,'b','Linewidth',1.5);
hold on
plot(Nr_plot,C1,'g','Linewidth',1.5);
hold on
plot(Nr_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of Rx antennas at user u (Nr)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('Nt=128, NtRF=16, NrRF=4, L=3, users=4, cells=3, SNR=5 - [2.9 min]')
set(gcf,'color','white')
toc

% ----------------------------------------------------------------------
%{
% (2) comparision with Nt(number of Tx antennas)
clear
close all
%clc
tic
%SNR_dB=[-30:5:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
Nr_RF=4; % number of RF chains in each user(u) at Rx=no. of data streams at each u

nt=256; % number of Tx antennas
Nt_plot=[32:32:nt];

Nr=8; % number of Rx antennas 
Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
    
L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_Nt=1:length(Nt_plot)
    i_Nt
    
    Nt=Nt_plot(i_Nt); % parameter (1) 
    Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
    ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
    ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
    M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side
    
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_Nt)= real(temp4/N_iter);
    C3(i_Nt)= real(temp3/N_iter);     
    C1(i_Nt)= real(temp1/N_iter);
    C2(i_Nt)= real(temp2/N_iter);
end

plot(Nt_plot,C4,'r','Linewidth',1.5);
hold on
plot(Nt_plot,C3,'b','Linewidth',1.5);
hold on
plot(Nt_plot,C1,'g','Linewidth',1.5);
hold on
plot(Nt_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of BS antennas (Nt)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('Nr=8, NtRF=16, NrRF=4, L=3, users=4, cells=3, SNR=5 - [2.7 min]')
set(gcf,'color','white')
toc
%}
% ----------------------------------------------------------------------
%{
% (3) comparision with Num_users (number of users per cell)
clear
close all
%clc
tic
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
%Num_users=4; % Number of users
Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 

Nt=128; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side

Nr=16; % number of Rx antennas 
  
L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 

users_plot=[1,2,4,8,16];

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_users=1:length(users_plot)
    i_users
    
    Num_users=users_plot(i_users);
    Nr_RF=Nt_RF/Num_users; 
    
    Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
    ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
    ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
    Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
      
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_users)= real(temp4/N_iter);
    C3(i_users)= real(temp3/N_iter);     
    C1(i_users)= real(temp1/N_iter);
    C2(i_users)= real(temp2/N_iter);
end

plot(users_plot,C4,'r','Linewidth',1.5);
hold on
plot(users_plot,C3,'b','Linewidth',1.5);
hold on
plot(users_plot,C1,'g','Linewidth',1.5);
hold on
plot(users_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of users per cell (U)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=128x16, NtRF=16, NrRF=NtRF/Num_users, L=3,, cells=3, SNR=5 - [5.4 min]')
set(gcf,'color','white')
toc
%}
% ----------------------------------------------------------------------
%{
% (4) comparision with Num_cells (number of cells)
clear
close all
%clc
tic
N_iter=30; % earlier taken for 100 iterations 

max_cells=6; % Maximum number of cells 
Num_users=4; % Number of users
Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
Nr_RF=4; % number of RF chains per user at RX = number of data streams

Nt=128; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side

Nr=8; % number of Rx antennas 
Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
      
L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 

cells_plot=[1:max_cells];

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_cells=1:length(cells_plot)
    i_cells
    
    Num_cells=cells_plot(i_cells); 
    
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_cells)= real(temp4/N_iter);
    C3(i_cells)= real(temp3/N_iter);     
    C1(i_cells)= real(temp1/N_iter);
    C2(i_cells)= real(temp2/N_iter);
end

plot(cells_plot,C4,'r','Linewidth',1.5);
hold on
plot(cells_plot,C3,'b','Linewidth',1.5);
hold on
plot(cells_plot,C1,'g','Linewidth',1.5);
hold on
plot(cells_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of cells (C)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=128x8, NtRF=16, NrRF=4, L=3,, users=4, SNR=5 - [5.4 min]')
set(gcf,'color','white')
toc
%}

% ----------------------------------------------------------------------
%{
% (5) comparision with L (number of paths/rays)
clear
close all
%clc
tic
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
Nr_RF=4; % number of RF chains per user at RX = number of data streams
SNR=5; % it is in normal scale (not in dB) 

Nt=128; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side

Nr=8; % number of Rx antennas 
Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
      
paths=15; % number of rays(paths)

L_plot=[1:paths];

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_L=1:length(L_plot)
    i_L
    
    L=L_plot(i_L); 
    
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_L)= real(temp4/N_iter);
    C3(i_L)= real(temp3/N_iter);     
    C1(i_L)= real(temp1/N_iter);
    C2(i_L)= real(temp2/N_iter);
end

plot(L_plot,C4,'r','Linewidth',1.5);
hold on
plot(L_plot,C3,'b','Linewidth',1.5);
hold on
plot(L_plot,C1,'g','Linewidth',1.5);
hold on
plot(L_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of paths/rays (L)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=128x8, NtRF=16, NrRF=4, users=4, cells=3, SNR=5 - [5.4 min]')
set(gcf,'color','white')
toc
%}
% ----------------------------------------------------------------------
%{
% (6) comparision with Nt_RF(number of RF chains at Tx side)
clear
close all
%clc
tic
%SNR_dB=[-30:5:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
%Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
%Nr_RF=4; % number of RF chains in each user(u) at Rx=no. of data streams at each u

Nt=256; % number of Tx antennas

Nt_RF_plot=[4,8,16,32,64];

Nr=32; % number of Rx antennas 
 
L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_Nt_RF=1:length(Nt_RF_plot)
    i_Nt_RF
    
    Nt_RF=Nt_RF_plot(i_Nt_RF); % parameter (1) 
    Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
    ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
    ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
    M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side
    
    Nr_RF=Nt_RF/Num_users; 
    Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
    ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
    ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
    Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
       
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_Nt_RF)= real(temp4/N_iter);
    C3(i_Nt_RF)= real(temp3/N_iter);     
    C1(i_Nt_RF)= real(temp1/N_iter);
    C2(i_Nt_RF)= real(temp2/N_iter);
end

plot(Nt_RF_plot,C4,'r','Linewidth',1.5);
hold on
plot(Nt_RF_plot,C3,'b','Linewidth',1.5);
hold on
plot(Nt_RF_plot,C1,'g','Linewidth',1.5);
hold on
plot(Nt_RF_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of RF chains at Tx (NtRF) in one cell')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=256x32, L=3, users=4, NrRF=NtRF/U, cells=3, SNR=5-[42 min]')
set(gcf,'color','white')
toc
%}
% ----------------------------------------------------------------------
%{
% (7) comparision with Nr_RF(number of RF chains at Rx side)
clear
close all
%clc
tic
%SNR_dB=[-30:5:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
%Nt_RF=16; % number of RF chains at TX = number of data streams at Tx 
%Nr_RF=4; % number of RF chains in each user(u) at Rx=no. of data streams at each u

Nt=256; % number of Tx antennas

Nr=64; % number of Rx antennas 

Nr_RF_plot=[1,2,4,8,16];
 
L=3; % number of rays(paths)
SNR=5; % it is in normal scale (not in dB) 

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

for i_Nr_RF=1:length(Nr_RF_plot)
    i_Nr_RF
    
    Nr_RF=Nr_RF_plot(i_Nr_RF); % parameter (1) 
    Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
    ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
    ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);
    Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side
    
    Nt_RF=Nr_RF*Num_users; %number of RF chains at TX = number of data streams
    Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
    ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
    ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
    M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side
      
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        
        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%    
        [Precoder4,Combiner4]= optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder4(:,:,cell,cell))*Precoder4(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner4(:,:,i_bs,i_bs))*Combiner4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder4(:,:,i_bs,i_bs))*Precoder4(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner4(:,:,i_bs,i_bs)));
        
        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3]= SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder3(:,:,cell,cell))*Precoder3(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner3(:,:,i_bs,i_bs))*Combiner3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder3(:,:,i_bs,i_bs))*Precoder3(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner3(:,:,i_bs,i_bs)));
         
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%     
        [Precoder1,Combiner1]= optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder1(:,:,cell,cell))*Precoder1(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner1(:,:,i_bs,i_bs))*Combiner1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder1(:,:,i_bs,i_bs))*Precoder1(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner1(:,:,i_bs,i_bs)));
 
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%    
        [Precoder2,Combiner2]= SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR);
        if Num_cells==1
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*0*0'*H_eff(:,:,i_bs,cell)'), 1:Num_cells, 'UniformOutput', false));
        else
            interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs,cell)*(Precoder2(:,:,cell,cell))*Precoder2(:,:,cell,cell)'*H_eff(:,:,i_bs,cell)'), setdiff(1:Num_cells, i_bs), 'UniformOutput', false));
        end
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2(:,:,i_bs,i_bs)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*Combiner2(:,:,i_bs,i_bs))*Combiner2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*(Precoder2(:,:,i_bs,i_bs))*Precoder2(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*Combiner2(:,:,i_bs,i_bs)));
 
    end
    C4(i_Nr_RF)= real(temp4/N_iter);
    C3(i_Nr_RF)= real(temp3/N_iter);     
    C1(i_Nr_RF)= real(temp1/N_iter);
    C2(i_Nr_RF)= real(temp2/N_iter);
end

plot(Nr_RF_plot,C4,'r','Linewidth',1.5);
hold on
plot(Nr_RF_plot,C3,'b','Linewidth',1.5);
hold on
plot(Nr_RF_plot,C1,'g','Linewidth',1.5);
hold on
plot(Nr_RF_plot,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('number of RF chains at user u (NrRF) of Rx side in one cell')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=256x64, L=3, users=4, NtRF=NrRF*U, cells=3, SNR=5-[38 min]')
set(gcf,'color','white')
toc
%}
% ----------------------------------------------------------------------
% -----------------------------------------------------------------------

%mmWave channel
function [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada)
power=sqrt(Nr_w*Nr_h*Nt_w*Nt_h/L);
H=zeros(Nr_w*Nr_h, Nt_w*Nt_h, Num_users,Num_cells,Num_cells);  % One user channel
H_eff=zeros(Num_users*Nr_w*Nr_h,Nt_w*Nt_h,Num_cells,Num_cells); % effective channel 
d=lamada/2;
Nt=Nt_w*Nt_h;  
Nr=Nr_w*Nr_h;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);

for i_bs=1:1:Num_cells
for j_ms=1:1:Num_cells
for u=1:1:Num_users
AoD_az(u,:)=2*pi*rand(1,L); %uniformly distributed in [0,2pi]
AoD_el(u,:)=pi*rand(1,L)-pi/2; %uniformly distributed in [-pi/2,pi/2]
AoA_az(u,:)=2*pi*rand(1,L); %uniformly distributed in [0,2pi]
AoA_el(u,:)=pi*rand(1,L)-pi/2; %uniformly distributed in [-pi/2,pi/2]
alpha(u,:)=sqrt(1/2)*(randn(1,L)+1j*randn(1,L));
power_matrix=power*diag(alpha(u,:));

for l=1:L
    A_BS(:,l)=array_respones(AoD_az(u,l),AoD_el(u,l),Nt,ind_Nt_w,ind_Nt_h,d,lamada);
    A_MS(:,l)=array_respones(AoA_az(u,l),AoA_el(u,l),Nr,ind_Nr_w,ind_Nr_h,d,lamada);
end
H(:,:,u,j_ms,i_bs)=A_MS*power_matrix*A_BS';
end
H_new=permute(H(:,:,:,j_ms,i_bs),[1 3 2]);
H_new=reshape(H_new,[],size(H(:,:,:,j_ms,i_bs),2),1);
H_eff(:,:,j_ms,i_bs)=H_new;
end
end
end

%array_response (for UPA)
function a=array_respones(azimuth,elevation,N,ind_N_w,ind_N_h,d,lamada)
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1j*(2*pi/lamada)*d*([ind_N_w]*sin(azimuth(i))*sin(elevation(i))+[ind_N_h]*cos(elevation(i))))).'];
end
end

% -----------------------------------------------------------------------
% ----------------------------------------------------------------------
% optimal_precoding - fully-connected at both Tx and Rx sides. 
function [P,C]=optimal_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff)
P=zeros(Nt,Nt_RF,Num_cells,Num_cells);
C=zeros(Num_users*Nr,Num_users*Nr_RF,Num_cells,Num_cells);
comb_u=zeros(Nr,Nr_RF,Num_users,Num_cells,Num_cells);
%H_eff=zeros(Num_users*Nr,Nt,Num_cells,Num_cells);
channel=zeros(Nr,Nt);

for i_bs=1:1:Num_cells

for user=1:1:Num_users
    channel(:,:)= H(:,:,user,i_bs,i_bs);
    [Uu,~,~]=svd(channel);
    comb_u(:,:,user,i_bs,i_bs)=Uu(:,1:Nr_RF); % hybrid combiner matrix of user u.  
end

c_ms_cell=num2cell(comb_u(:,:,:,i_bs,i_bs),[1,2,4]);
C(:,:,i_bs,i_bs)=blkdiag(c_ms_cell{:});

[~,~,V]=svd(H_eff(:,:,i_bs,i_bs));
P(:,:,i_bs,i_bs)=V(:,1:Nt_RF); % P=precoder matrix at Tx side 

end
end

% --------------------------------------------------------------------------------

% SIC-based_hybrid_precoding using fully-connected structure at both Tx and Rx sides
function [P,C]=SIC_based_precoding_fully_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,H,H_eff)
P=zeros(Nt,Nt_RF,Num_cells,Num_cells); 
C=zeros(Num_users*Nr,Num_users*Nr_RF,Num_cells,Num_cells);
comb_u=zeros(Nr,Nr_RF,Num_users,Num_cells,Num_cells);
channel=zeros(Nr,Nt);

for i_bs=1:1:Num_cells

for user=1:1:Num_users
channel(:,:)= H(:,:,user,i_bs,i_bs);
[Uu,~,~]=svd(channel);
a_ms=exp(1i*angle(Uu(:,1:Nr_RF)))/sqrt(Nr); % phase=a_m %since a_m amplitude=1/sqrt(M)
d_ms=norm(Uu(:,1:Nr_RF),1)/sqrt(Nr);
comb_u(:,:,user,i_bs,i_bs)=a_ms*d_ms; % hybrid combiner matrix of user u 
end

c_ms_cell=num2cell(comb_u(:,:,:,i_bs,i_bs),[1,2,4]);
C(:,:,i_bs,i_bs)=blkdiag(c_ms_cell{:});

[~,~,V]=svd(H_eff(:,:,i_bs,i_bs));
a_bs=exp(1i*angle(V(:,1:Nt_RF)))/sqrt(Nt); 
d_bs=norm(V(:,1:Nt_RF),1)/sqrt(Nt);
P(:,:,i_bs,i_bs)=a_bs*d_bs; % P=precoder matrix at Tx side

end
end

% -------------------------------------------------------------------------------

% optimal_precoding - sub-connected at both Tx and Rx sides. 
function [P,C]=optimal_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR)
P=zeros(Nt,Nt_RF,Num_cells,Num_cells); % effective precoding matrix 
P_others=zeros(Nt,Nt_RF,Num_cells,Num_cells); % other Precoding matrices causing interference
C=zeros(Num_users*Nr,Num_users*Nr_RF,Num_cells,Num_cells); % effective combiner
comb_u=zeros(Nr,Nr_RF,Num_users,Num_cells,Num_cells); % Combiner of user u
%H_eff=zeros(Num_users*Nr,Nt,Num_cells,Num_cells);
channel=zeros(Nr,Nt);

for i_bs=1:1:Num_cells

for user=1:1:Num_users
channel(:,:)= H(:,:,user,i_bs,i_bs);
for j=1:Nr_RF
    Gu=channel*inv(eye(Nt)+(SNR/Nr_RF)*channel'*comb_u(:,1:(j-1),user,i_bs,i_bs)*comb_u(:,1:(j-1),user,i_bs,i_bs)'*channel)*channel';
    cc=zeros(Nr,1); %cc=combiner matrix at user u
    Su=Gu(Mu*(j-1)+1:Mu*(j-1)+Mu,Mu*(j-1)+1:Mu*(j-1)+Mu);  
    [Uu,~,~]=svd(Su);
    u1=Uu(:,1); 
    cc(Mu*(j-1)+1:Mu*(j-1)+Mu)=u1; 
    comb_u(:,j,user,i_bs,i_bs)=cc; % hybrid combiner matrix of user u.  
end
end

c_ms_cell=num2cell(comb_u(:,:,:,i_bs,i_bs),[1,2,4]);
C(:,:,i_bs,i_bs)=blkdiag(c_ms_cell{:});

for i=1:Nt_RF
    G=H_eff(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs)*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs))*C(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*P_others(:,1:(i-1),i_bs,i_bs)*P_others(:,1:(i-1),i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs))*C(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs);
    f=zeros(Nt,1); %f=precoding matrix of 1 RF chain at Tx side 
    S=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    a_bs=exp(1i*angle(v1))/sqrt(M); % phase=a_m %since a_bs amplitude=1/sqrt(M)
    d_bs=norm(v1,1)/sqrt(M);
    f(M*(i-1)+1:M*(i-1)+M)=a_bs*d_bs; 
    P_others(:,i,i_bs,i_bs)=f; % P=precoder matrix at Tx side 
end
end

for i_bs_n=1:1:Num_cells
for i_n=1:Nt_RF
    if Num_cells==1
        interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs_n,cell)*0*0'*H_eff(:,:,i_bs_n,cell)'), 1:Num_cells, 'UniformOutput', false));
    else
        interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs_n,cell)*(P_others(:,:,cell,cell))*P_others(:,:,cell,cell)'*H_eff(:,:,i_bs_n,cell)'), setdiff(1:Num_cells, i_bs_n), 'UniformOutput', false));
    end
    G=H_eff(:,:,i_bs_n,i_bs_n)'*C(:,:,i_bs_n,i_bs_n)*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C(:,:,i_bs_n,i_bs_n)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*C(:,:,i_bs_n,i_bs_n))*C(:,:,i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n)*P(:,1:(i_n-1),i_bs_n,i_bs_n)*P(:,1:(i_n-1),i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n)'*C(:,:,i_bs_n,i_bs_n))*C(:,:,i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n);
    f=zeros(Nt,1); %f=precoding matrix of 1 RF Rf chain at Tx side 
    S=G(M*(i_n-1)+1:M*(i_n-1)+M,M*(i_n-1)+1:M*(i_n-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    f(M*(i_n-1)+1:M*(i_n-1)+M)=v1; 
    P(:,i_n,i_bs_n,i_bs_n)=f; % P=precoder matrix of 1 RF chain at Tx side  
end
end
end

%----------------------------------------------------------------------------
% SIC-based_hybrid_precoding using sub-connected structure at both Tx and Rx sides
function [P,C]=SIC_based_precoding_sub_connected(Num_cells,Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,H_eff,SNR)
P=zeros(Nt,Nt_RF,Num_cells,Num_cells); % effective precoding matrix 
P_others=zeros(Nt,Nt_RF,Num_cells,Num_cells); % other Precoding matrices causing interference
C=zeros(Num_users*Nr,Num_users*Nr_RF,Num_cells,Num_cells); % effective combiner
comb_u=zeros(Nr,Nr_RF,Num_users,Num_cells,Num_cells); % Combiner of user u
channel=zeros(Nr,Nt);

for i_bs=1:1:Num_cells
%considering sub-connected structure at the Rx side 
for user=1:1:Num_users
channel(:,:)= H(:,:,user,i_bs,i_bs); 
for j=1:Nr_RF
    Gu=channel*inv(eye(Nt)+(SNR/Nr_RF)*channel'*comb_u(:,1:(j-1),user,i_bs,i_bs)*comb_u(:,1:(j-1),user,i_bs,i_bs)'*channel)*channel';
    cc=zeros(Nr,1); %cc=combiner matrix of 1 RF chain at Rx side 
    Su=Gu(Mu*(j-1)+1:Mu*(j-1)+Mu,Mu*(j-1)+1:Mu*(j-1)+Mu); %temp=S 
    [Uu,~,~]=svd(Su);
    u1=Uu(:,1); 
    
    a_ms=exp(1i*angle(u1))/sqrt(Mu); % phase=a_m %since a_m amplitude=1/sqrt(M)
    d_ms=norm(u1,1)/sqrt(Mu);
    cc(Mu*(j-1)+1:Mu*(j-1)+Mu)=a_ms*d_ms;
    comb_u(:,j,user,i_bs,i_bs)=cc; % F=P=combiner matrix of user u
    
end
end

c_ms_cell=num2cell(comb_u(:,:,:,i_bs,i_bs),[1,2,4]);
C(:,:,i_bs,i_bs)=blkdiag(c_ms_cell{:}); 

for i=1:Nt_RF
    G=H_eff(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs)*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs))*C(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)*P_others(:,1:(i-1),i_bs,i_bs)*P_others(:,1:(i-1),i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs)'*C(:,:,i_bs,i_bs))*C(:,:,i_bs,i_bs)'*H_eff(:,:,i_bs,i_bs);
    f=zeros(Nt,1); %f=precoding matrix of 1 RF chain at Tx side 
    S=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    a_bs=exp(1i*angle(v1))/sqrt(M); % phase=a_m %since a_bs amplitude=1/sqrt(M)
    d_bs=norm(v1,1)/sqrt(M);
    f(M*(i-1)+1:M*(i-1)+M)=a_bs*d_bs; 
    P_others(:,i,i_bs,i_bs)=f; % P=precoder matrix at Tx side 
end
end

for i_bs_n=1:1:Num_cells
for i_n=1:Nt_RF
    if Num_cells==1
        interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs_n,cell)*0*0'*H_eff(:,:,i_bs_n,cell)'), 1:Num_cells, 'UniformOutput', false));
    else
        interference=(arrayfun(@(cell) ((SNR/(Num_users*Nr))*H_eff(:,:,i_bs_n,cell)*(P_others(:,:,cell,cell))*P_others(:,:,cell,cell)'*H_eff(:,:,i_bs_n,cell)'), setdiff(1:Num_cells, i_bs_n), 'UniformOutput', false));
    end
    G=H_eff(:,:,i_bs_n,i_bs_n)'*C(:,:,i_bs_n,i_bs_n)*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C(:,:,i_bs_n,i_bs_n)'*(eye(Num_users*Nr)+(sum(cat(3,interference{:}),3)))*C(:,:,i_bs_n,i_bs_n))*C(:,:,i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n)*P(:,1:(i_n-1),i_bs_n,i_bs_n)*P(:,1:(i_n-1),i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n)'*C(:,:,i_bs_n,i_bs_n))*C(:,:,i_bs_n,i_bs_n)'*H_eff(:,:,i_bs_n,i_bs_n);
    f=zeros(Nt,1); %f=precoding matrix of 1 RF Rf chain at Tx side 
    S=G(M*(i_n-1)+1:M*(i_n-1)+M,M*(i_n-1)+1:M*(i_n-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    a_bs=exp(1i*angle(v1))/sqrt(M); % phase=a_m %since a_bs amplitude=1/sqrt(M)
    d_bs=norm(v1,1)/sqrt(M);
    f(M*(i_n-1)+1:M*(i_n-1)+M)=a_bs*d_bs; 
    P(:,i_n,i_bs_n,i_bs_n)=f; % P=precoder matrix of 1 RF chain at Tx side  
end
end
end
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------