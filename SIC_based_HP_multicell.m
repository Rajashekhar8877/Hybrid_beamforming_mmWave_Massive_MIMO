% Multi-cell case 
% Test_2 (Extra for rough work)

clear
close all
%clc
tic
SNR_dB=[-30:5:0];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_cells=3; % Number of cells 
Num_users=4; % Number of users
Nt_RF=8; % number of RF chains at TX = number of data streams at Tx 
Nr_RF=2; % number of RF chains in each user(u) at Rx=no. of data streams at each u

Nt=64; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);

Nr=8; % number of Rx antennas 
Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);

L=3; % number of rays(paths)
M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side
Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side

%cell=1;  
i_bs=1;  % cell for which capacity has to be found

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;

for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    
    temp1=0;temp2=0;temp3=0; temp4=0;
    
    for iter=1:N_iter
        [H,H_eff,A_BS,A_MS]=mmWave_channel(Num_cells,Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);
        %syms i_bs 
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
    C4(i_snr)= real(temp4/N_iter);
    C3(i_snr)= real(temp3/N_iter);     
    C1(i_snr)= real(temp1/N_iter);
    C2(i_snr)= real(temp2/N_iter);
end

plot(SNR_dB,C4,'r','Linewidth',1.5);
hold on
plot(SNR_dB,C3,'b','Linewidth',1.5);
hold on
plot(SNR_dB,C1,'g','Linewidth',1.5);
hold on
plot(SNR_dB,C2,'y','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)')
xlabel('SNR (dB)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=128x8, NtRF=16, NrRF=4, L=3, users=4, cells=3 - [2.7 min]')
%whitebg(figure,'white')
set(gcf,'color','white')
toc

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
%H_eff=zeros(Num_users*Nr,Nt,Num_cells,Num_cells);
channel=zeros(Nr,Nt);
%H_eff=zeros(Num_users*Nr,Nt);

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
%    a_bs=exp(1i*angle(v1))/sqrt(M); % phase=a_m %since a_bs amplitude=1/sqrt(M)
%    d_bs=norm(v1,1)/sqrt(M);
    f(M*(i-1)+1:M*(i-1)+M)=v1; 
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
P=zeros(Nt,Nt_RF,Num_cells,Num_cells); % effective precoding matrix (of cell k)
P_others=zeros(Nt,Nt_RF,Num_cells,Num_cells); % other Precoding matrices causing interference(e)
C=zeros(Num_users*Nr,Num_users*Nr_RF,Num_cells,Num_cells); % effective combiner (of cell k)
comb_u=zeros(Nr,Nr_RF,Num_users,Num_cells,Num_cells); % Combiner of user u (in cell k)
%H_eff=zeros(Num_users*Nr,Nt,Num_cells,Num_cells);
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