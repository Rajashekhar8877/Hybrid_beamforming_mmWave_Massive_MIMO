% GMD based hybrid precoding  

clear
close all
%clc
tic
SNR_dB=[-25:50:25];  SNR_linear=10.^(SNR_dB/10.);
N_iter=30; % earlier taken for 100 iterations 

Num_users=1; % Number of users
Nt_RF=8; %number of RF chains at TX = number of data streams at Tx 
Nr_RF=8; % number of RF chains in each user(u) at Rx=no. of data streams at each u

Nt=128; % number of Tx antennas
Nt_w=Nt_RF; Nt_h=Nt/Nt_RF;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);

Nr=128; % number of Rx antennas 
Nr_w=Nr_RF; Nr_h=Nr/Nr_RF;
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);

L=3; % number of rays(paths)
M=Nt/Nt_RF; % number of antennas connected to one RF chains at Tx side
Mu=Nr/Nr_RF; % number of antennas connected to one RF chains at Rx side

fc=28e9; % Frequencey 
lamada=3e8/fc; % wavelegenth;
for i_snr=1:length(SNR_linear)
    i_snr
    SNR=SNR_linear(i_snr);
    
    temp1=0;temp2=0;temp3=0; temp4=0; temp=5;
    
    for i=1:N_iter
        [H,A_BS,A_MS]=mmWave_channel(Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada);

        %%%%%% Optimal precoding - fully connected %%%%%%%%%%%%%%%%%%%%
        [Precoder4,Combiner4,H_eff]= optimal_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H);
        temp4=temp4+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner4'*Combiner4)*Combiner4'*H_eff*(Precoder4)*Precoder4'*H_eff'*Combiner4));

        %%%%%% SIC-based hybrid precoding - fully connected %%%%%%%%%%%%%%%
        [Precoder3,Combiner3,H_eff]= SIC_based_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H);
        temp3=temp3+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner3'*Combiner3)*Combiner3'*H_eff*(Precoder3)*Precoder3'*H_eff'*Combiner3));
        
        %%%%%% Optimal precoding - sub connected %%%%%%%%%%%%%%%%%%%%
        [Precoder1,Combiner1,H_eff]= optimal_precoding_sub_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,SNR);
        temp1=temp1+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner1'*Combiner1)*Combiner1'*H_eff*(Precoder1)*Precoder1'*H_eff'*Combiner1));
        
        %%%%%% SIC-based hybrid precoding - sub connected %%%%%%%%%%%%%%%
        [Precoder2,Combiner2,H_eff]= SIC_based_precoding_sub_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,SNR);
        temp2=temp2+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner2'*Combiner2)*Combiner2'*H_eff*(Precoder2)*Precoder2'*H_eff'*Combiner2));
                 
        %%%%%% GMD-based hybrid precoding - %%%%%%%%%%%%%%%
        [Precoder5,Combiner5,H_eff]= GMD_based_precoding(Num_users,Nt_RF,Nr_RF,Nt,Nr,H, A_BS, A_MS);
        temp5=temp5+log2(det(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(Combiner5'*Combiner5)*Combiner5'*H_eff*(Precoder5)*Precoder5'*H_eff'*Combiner5));
  
    end
    C4(i_snr)= real(temp4/N_iter);
    C3(i_snr)= real(temp3/N_iter);     
    C1(i_snr)= real(temp1/N_iter);
    C2(i_snr)= real(temp2/N_iter);
    C5(i_snr)= real(temp5/N_iter);
end

plot(SNR_dB,C4,'r','Linewidth',1.5);
hold on
plot(SNR_dB,C3,'b','Linewidth',1.5);
hold on
plot(SNR_dB,C1,'g','Linewidth',1.5);
hold on
plot(SNR_dB,C2,'y','Linewidth',1.5);
hold on
plot(SNR_dB,C5,'c','Linewidth',1.5);

legend('optimal precoding(fully-connected)','SIC-based hybrid precoding(fully-connected)','optimal precoding(sub-connected)','SIC-based hybrid precoding(sub-connected)','GMD based hybrid precoding')
xlabel('SNR (dB)')
ylabel('Achievable rate (bps/Hz)')
grid on 
title('NtxNr=128x8, NtRF=16, NrRF=4, L=3, users=4 - [2.7 min]')
%whitebg(figure,'white')
set(gcf,'color','white')
toc

%mmWave channel
function [H,A_BS,A_MS]=mmWave_channel(Num_users,Nr_w,Nr_h,Nt_w,Nt_h,L,lamada)
power=sqrt(Nr_w*Nr_h*Nt_w*Nt_h/L);
H=zeros(Nr_w*Nr_h, Nt_w*Nt_h, Num_users);  % One user channel
d=lamada/2;
Nt=Nt_w*Nt_h;  
Nr=Nr_w*Nr_h;
ind_Nt_w=reshape(repmat([0:1:Nt_w-1],Nt_h,1),1,Nt_w*Nt_h);
ind_Nt_h=repmat([0:1:Nt_h-1],1,Nt_w);
ind_Nr_w=reshape(repmat([0:1:Nr_w-1],Nr_h,1),1,Nr_w*Nr_h);
ind_Nr_h=repmat([0:1:Nr_h-1],1,Nr_w);

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
H(:,:,u)=A_MS*power_matrix*A_BS';
end
end

%array_response (for UPA)
function a=array_respones(azimuth,elevation,N,ind_N_w,ind_N_h,d,lamada)
a=[];
for i=1:length(azimuth)
    a=[a (sqrt(1/N)*exp(1j*(2*pi/lamada)*d*([ind_N_w]*sin(azimuth(i))*sin(elevation(i))+[ind_N_h]*cos(elevation(i))))).'];
end
end

% ---------------------------------------------------------------------------------

% GMD based hybrid precoding 
function [P,C,H_eff]=GMD_based_precoding(Num_users,Nt_RF,Nr_RF,Nt,Nr,H, A_BS, A_MS)

H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);

% hybrid precoder 
[U, S, V]=svd(H_eff);
U1=U(:,1:Nt_RF);
V1=V(:,1:Nt_RF);
[G, R, Q]=gmd(U, S, V, 0.0001);
Q1=Q(:,1:Nt_RF);
At=A_BS;

% OMP algorithm 
Qres=V1;
xee=algo_omp(Nt_RF, At, Q1);
%xee=algo_omp(Nt_RF, At, Qres);
%xee=algo_omp(Nt_RF, Q1, At);
%{
A1=[];
if abs(A1)~=Nt_RF
    [~,j]=max(abs(At'*(Qres)*Qres'*At));
    A1=[A1 j]; xee=A_BS(:,j);
    Y=pinv(xee)*Q1;
    Qres=(Q1-xee*Y)/norm((Q1-xee*Y),'fro'); 
end    
%}

QA=xee;
QD=pinv(QA)*Q1;
QD=sqrt(Nt_RF)*(QD)/norm((QA*QD), 'fro');

P=QA*QD;


% hybrid combiner  
channel=zeros(Nr,Nt);
for user=1:1:Num_users
    channel(:,:)= H(:,:,user);
    [Ur, Sr, Vr]=svd(channel);
U1r=Ur(:,1:Nr_RF);
V1r=Vr(:,1:Nr_RF);
[Gr, Rr, Qr]=gmd(Ur, Sr, Vr, 0.0001);
Q1r=Qr(:,1:Nr_RF);
Atr=A_MS;

% OMP algorithm 
Qresr=V1r;
xeer=algo_omp(Nr_RF, Atr, Q1r);
%xeer=algo_omp(Nr_RF, Atr, Qresr);
%xeer=algo_omp(Nr_RF, Q1r, Atr);
%{
A1r=[];
if abs(A1r)~=Nr_RF
    [~,jr]=max(abs(Atr'*(Qresr)*Qresr'*Atr));
    A1r=[A1r jr]; xeer=A_MS(:,j);
    Yr=pinv(xeer)*Q1r;
    Qresr=(Q1r-xeer*Yr)/norm((Q1r-xeer*Yr),'fro');
    
end    
%}

QAr=xeer;
QDr=pinv(QAr)*Q1r;
QDr=sqrt(Nr_RF)*(QDr)/norm((QAr*QDr), 'fro');

C=QAr*QDr;
end

%{
function  x = algo_omp(k, A, y)
    xbeg = zeros(size(A,2),1);
    support=[];
    temp=y;
    count = 1;
    while count < k+1
        ST = abs(A' * temp);
        [a, b] = max(ST);
        support = [support b];
        xfinal = A(:, support)\y;
        temp = y-A(:,support) * xfinal;
        count = count + 1;
    end
    x = xbeg;
    t = support';
    x(t) = xfinal;
end
%}
end

%{
% OMP algorithm from mathworks 
function [x] = algo_omp (K,y,A)
Res = y.' ;
[m,n] = size (A) ;
Q = zeros (m,K) ;
R = zeros (K,K) ;
Rinv = zeros (K,K) ;
w = zeros (m,K) ;
x = zeros (1,n) ;
for J = 1 : K
    
    %Index Search
    [V ,kkk] = max(abs(A'*Res)) ;
    kk (J) = kkk ;
    
    %Residual Update
    w (:,J) = A (:,kk (J)) ;
    for I = 1 : J-1
        if (J-1 ~= 0)
            R (I,J) = Q (:,I)' * w (:,J) ;
            w (:,J) = w (:,J) - R (I,J) * Q (:,I) ;
        end
    end
    R (J,J) = norm (w (:,J)) ;
    Q (:,J) = w (:,J) / R (J,J) ;
    Res = Res - (Q (:,J) * Q (:,J)' * Res) ;
  
end
%Least Squares
for J = 1 : K
    Rinv (J,J) = 1 / R (J,J) ;
    if (J-1 ~= 0)
        for I = 1 : J-1
            Rinv (I,J) = -Rinv (J,J) * (Rinv (I,1:J-1) * R (1:J-1,J)) ;
        end
    end
end
xx = Rinv * Q' * y.' ;
for I = 1 : K
    x (kk (I)) = xx (I) ;
end
end
%}


% OMP algorithm from Github 
function  x = algo_omp(k, A, y)
    xbeg = zeros(size(A,2),1);
    support=[];
    temp=y;
    count = 1;
    while count < k+1
        ST = abs(A' * temp);
        [a, b] = max(ST);
        support = [support b];
        xfinal = A(:, support)\y;
        temp = y-A(:,support) * xfinal;
        count = count + 1;
    end
    x = xbeg;
    t = support';
    x(t) = xfinal;
end



% Matlab implementation of the "Geometric Mean Decomposition"
% version of Hager, December 3, 2003
% slightly modified by Yi, April 19, 2004
% Copyright 2003, University of Florida, Gainesville, Florida
% 
%A = U*S*V' is the singular value decomposition of A
%           U, V unitary, S diagonal matrix with nonnegative
%           diagonal entries in decreasing order
%  = Q*R*P' is the geometric mean decomposition of A
%           P, Q unitary, R real upper triangular with r_ii =
%           geometric mean of the positive singular values of A,
%           1 <= i <= p, p = number of positive singular values
% All singular values smaller than tol treated as zero

function [Q, R, P] = gmd (U, S , V, tol)
if ( nargin < 4 )
%     tol = eps ;
    tol = 0.0001;
end
[m n] = size (S) ;
R = zeros (m, n) ;
P = V ;
Q = U ;
d = diag (S) ;
l = min (m, n) ;
for p = l : -1 : 1
    if ( d (p) >= tol )
        break ;
    end
end
if ( p < 1 )
    return ;
end
if ( p < 2 )
    R (1, 1) = d (1) ;
    return ;
end
z = zeros (p-1, 1) ;
large = 2 ;           % largest diagonal element
small = p ;           % smallest diagonal element
perm = [1 : p] ;      % perm (i) = location in d of i-th largest entry
invperm = [ 1 : p ] ; % maps diagonal entries to perm
sigma_bar = (prod (d (1:p)))^(1/p) ;
for k = 1 : p-1
    flag = 0 ;
    if ( d (k) >= sigma_bar )    
        i = perm (small) ;
        small = small - 1 ;
        if ( d (i) >= sigma_bar )
            flag = 1 ;
        end
    else
        i = perm (large) ;
        large = large + 1 ;
        if ( d (i) <= sigma_bar )
            flag = 1 ;
        end
    end
        
    k1 = k + 1 ;
    if ( i ~= k1 )            % Apply permutation Pi of paper
        t = d (k1) ;          % Interchange d (i) and d (k1)
        d (k1) = d (i) ;
        d (i) = t ;
        j = invperm (k1) ;    % Update perm arrays
        perm (j) = i ;
        invperm (i) = j ;
        I = [ k1 i ] ;
        J = [ i k1 ] ;
        Q (:, I) = Q (:, J) ; % interchange columns i and k+1
        P (:, I) = P (:, J) ;
    end

    delta1 = d (k) ;
    delta2 = d (k1) ;
    t = delta1 + delta2 ;
    if ( flag )
        c = 1 ;
        s = 0 ;
    else
        f = (delta1 - sigma_bar)/(delta1 - delta2) ;
        s = sqrt (f*(delta1+sigma_bar)/t) ;
        c = sqrt(1-s^2) ;
    end
    d (k1) = delta1*delta2/sigma_bar ;          % = y in paper
    z (k) = s*c*(delta2 - delta1)*t/sigma_bar ; % = x in paper
    R (k, k) = sigma_bar ;
    if ( k > 1 )
        R (1:k-1, k) = z (1:k-1)*c ; % new column of R
        z (1:k-1) = -z (1:k-1)*s ;   % new column of Z
    end

    G1 = [ c -s
           s  c ] ;
    J = [ k k1 ] ;
    P (:, J) = P (:, J)*G1 ;        % apply G1 to P

    G2 = (1/sigma_bar)*[ c*delta1 -s*delta2
                         s*delta2  c*delta1 ] ;
    Q (:, J) = Q (:, J)*G2 ;        % apply G2 to Q
end

R (p, p) = sigma_bar ;
R (1:p-1, p) = z ;
end

% -----------------------------------------------------------------------

% optimal_precoding - fully-connected at both Tx and Rx sides. 
function [P,C,H_eff]=optimal_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
%P=zeros(Nt,Nt_RF); 
comb_u=zeros(Nr,Nr_RF,Num_users);
%H_eff=zeros(Num_users*Nr,Nt,Num_users);
channel=zeros(Nr,Nt);
H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);

for user=1:1:Num_users
    channel(:,:)= H(:,:,user);
    [Uu,~,~]=svd(channel); 
    comb_u(:,:,user)=Uu(:,1:Nr_RF); % hybrid combiner matrix of user u.  
end

c_ms_cell=num2cell(comb_u,[1,2,4]);
C=blkdiag(c_ms_cell{:});

[~,~,V]=svd(H_eff);
P=V(:,1:Nt_RF); % P=precoder matrix at Tx side 

end

% --------------------------------------------------------------------------------

% ------------------------------------------------------------------------------

% SIC-based_hybrid_precoding using fully-connected structure at both Tx and Rx sides
function [P,C,H_eff]=SIC_based_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
%P=zeros(Nt,Nt_RF); 
comb_u=zeros(Nr,Nr_RF,Num_users);
%H_eff=zeros(Num_users*Nr,Nt,Num_users);
channel=zeros(Nr,Nt);
H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);
%a_ms=zeros(Nr,Nr_RF,Num_users);
%d_ms=zeros(Nr_RF,Nr_RF,Num_users);

for user=1:1:Num_users
channel(:,:)= H(:,:,user);
[Uu,~,~]=svd(channel);
a_ms=exp(1i*angle(Uu(:,1:Nr_RF)))/sqrt(Nr); % phase=a_m %since a_m amplitude=1/sqrt(M)
d_ms=norm(Uu(:,1:Nr_RF),1)/sqrt(Nr);
comb_u(:,:,user)=a_ms*d_ms; % hybrid combiner matrix of user u 
end

c_ms_cell=num2cell(comb_u,[1,2,4]);
C=blkdiag(c_ms_cell{:});

[~,~,V]=svd(H_eff); 
a_bs=exp(1i*angle(V(:,1:Nt_RF)))/sqrt(Nt); 
d_bs=norm(V(:,1:Nt_RF),1)/sqrt(Nt); 
P=a_bs*d_bs; % P=precoder matrix at Tx side

end


% -------------------------------------------------------------------------------

% optimal_precoding - sub-connected at both Tx and Rx sides. 
function [P,C,H_eff]=optimal_precoding_sub_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,SNR)
P=zeros(Nt,Nt_RF); 
comb_u=zeros(Nr,Nr_RF,Num_users);
%H_eff=zeros(Num_users*Nr,Nt,Num_users);
channel=zeros(Nr,Nt);
H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);

for user=1:1:Num_users
channel(:,:)= H(:,:,user);
for j=1:Nr_RF
    Gu=channel*inv(eye(Nt)+(SNR/Nr_RF)*channel'*comb_u(:,1:(j-1),user)*comb_u(:,1:(j-1),user)'*channel)*channel';
    cc=zeros(Nr,1); %cc=combiner matrix at user u
    Su=Gu(Mu*(j-1)+1:Mu*(j-1)+Mu,Mu*(j-1)+1:Mu*(j-1)+Mu);  
    [Uu,~,~]=svd(Su);
    u1=Uu(:,1); 
    cc(Mu*(j-1)+1:Mu*(j-1)+Mu)=u1; 
    comb_u(:,j,user)=cc; % hybrid combiner matrix of user u.  
end
end

c_ms_cell=num2cell(comb_u,[1,2,4]);
C=blkdiag(c_ms_cell{:});

for i=1:Nt_RF
    G=H_eff'*C*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C'*C)*C'*H_eff*P(:,1:(i-1))*P(:,1:(i-1))'*H_eff'*C)*C'*H_eff;
%    G=H_eff'*C*((eye(Num_users*Nr_RF)+(SNR/Nt_RF)*((C'*C)\(C'*H_eff*P(:,1:(i-1))*P(:,1:(i-1))'*H_eff'*C)))\C')*H_eff;
    f=zeros(Nt,1); %f=precoding matrix of 1 RF Rf chain at Tx side 
    S=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    f(M*(i-1)+1:M*(i-1)+M)=v1; 
    P(:,i)=f; % P=precoder matrix of 1 RF chain at Tx side  
end
end

%----------------------------------------------------------------------------

% SIC-based_hybrid_precoding using sub-connected structure at both Tx and Rx sides
function [P,C,H_eff]=SIC_based_precoding_sub_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,M,Mu,H,SNR)
P=zeros(Nt,Nt_RF); %S1=zeros(M,M,Nt_RF,Num_users);
comb_u=zeros(Nr,Nr_RF,Num_users);
%H_eff=zeros(Num_users*Nr,Nt,Num_users);
channel=zeros(Nr,Nt);
H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);
%a_ms=zeros(Nr,Nr_RF,Num_users);
%d_ms=zeros(Nr_RF,Nr_RF,Num_users);

%considering sub-connected structure at the Rx side 
for user=1:1:Num_users
channel(:,:)= H(:,:,user);
for j=1:Nr_RF
    Gu=channel*inv(eye(Nt)+(SNR/Nr_RF)*channel'*comb_u(:,1:(j-1),user)*comb_u(:,1:(j-1),user)'*channel)*channel';
    cc=zeros(Nr,1); %cc=combiner matrix of 1 RF chain at Rx side 
    Su=Gu(Mu*(j-1)+1:Mu*(j-1)+Mu,Mu*(j-1)+1:Mu*(j-1)+Mu); %temp=S 
    [Uu,~,~]=svd(Su);
    u1=Uu(:,1); 
    
    a_ms=exp(1i*angle(u1))/sqrt(Mu); % phase=a_m %since a_m amplitude=1/sqrt(M)
    d_ms=norm(u1,1)/sqrt(Mu);
    cc(Mu*(j-1)+1:Mu*(j-1)+Mu)=a_ms*d_ms;
    comb_u(:,j,user)=cc; % F=P=combiner matrix of user u
    
end
end

c_ms_cell=num2cell(comb_u,[1,2,4]);
C=blkdiag(c_ms_cell{:});

for i=1:Nt_RF
    G=H_eff'*C*inv(eye(Num_users*Nr_RF)+(SNR/(Num_users*Nr_RF))*inv(C'*C)*C'*H_eff*P(:,1:(i-1))*P(:,1:(i-1))'*H_eff'*C)*C'*H_eff;
%    G=H_eff'*C*((eye(Num_users*Nr_RF)+(SNR/Nt_RF)*((C'*C)\(C'*H_eff*P(:,1:(i-1))*P(:,1:(i-1))'*H_eff'*C)))\C')*H_eff;
    f=zeros(Nt,1); %f=precoding matrix of 1 RF chain at Tx side 
    S=G(M*(i-1)+1:M*(i-1)+M,M*(i-1)+1:M*(i-1)+M); 
    [~,~,V]=svd(S);
    v1=V(:,1); 
    a_bs=exp(1i*angle(v1))/sqrt(M); % phase=a_m %since a_bs amplitude=1/sqrt(M)
    d_bs=norm(v1,1)/sqrt(M);
    f(M*(i-1)+1:M*(i-1)+M)=a_bs*d_bs; 
    P(:,i)=f; % P=precoder matrix at Tx side 
end
end

%-------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%----------------------------------------------------------------------------------
%{
% other method 1 - not correct as it is going above the optimal precoding
% for 1 user 
% SIC-based_hybrid_precoding using SVD by using fully-connected structure at both Tx and Rx sides 
%function [a_bs,D_bs,A_ms,D_ms,H_eff]=SIC_based_precoding(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
function [P,C,H_eff]=SIC_based_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
channel=zeros(Nr, Nt);
Hu_eff=zeros(Nr_RF,Nt_RF,Num_users);
a_ms=zeros(Nr,Nr_RF,Num_users);
for user=1:1:Num_users
channel(:,:)= H(:,:,user);
[Uu,~,~]=svd(channel);
%[Uu,~,~]=svd(channel/sqrt(Nt));
a_ms(:,:,user)=exp(1i*angle(Uu(:,1:Nr_RF)))/sqrt(Nr);
end
a_ms_cell=num2cell(a_ms,[1,2,4]);
A_ms=blkdiag(a_ms_cell{:});

H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);
%{
H_comp = A_ms'*H_eff;
[~,~,V_comp]=svd(H_comp/sqrt(Nt));
a_bs=exp(1i*angle(V_comp(:,1:Nt_RF)))/sqrt(Nt);
%}

[~,~,V_eff]=svd(H_eff);
%[~,~,V_eff]=svd(H_eff/sqrt(Nt));
a_bs=exp(1i*angle(V_eff(:,1:Nt_RF)))/sqrt(Nt);

H_comp_eff=A_ms'*H_eff*a_bs;
[~,~,V_comp_eff]=svd(H_comp_eff);
d_bs=V_comp_eff(:,1:Nt_RF);

d_ms=zeros(Nr_RF,Nr_RF,Num_users);

for users=1:1:Num_users
    Hu_eff(:,:,users)=a_ms(:,:,users)'*H(:,:,users)*a_bs;
    [Uu_eff,~,~]=svd(Hu_eff(:,:,users));
    d_ms(:,:,users)=Uu_eff(:,1:Nr_RF);
end

d_ms_cell=num2cell(d_ms,[1,2,4]);
D_ms=blkdiag(d_ms_cell{:});
    
P=a_bs*d_bs;
C=A_ms*D_ms;
end
%}
% -------------------------------------------------------------------------------


%------------------------------------------------------------------------------
%----------------------------------------------------------------------------------
%{
% other method 2 (according to paper) - it is also going above the optimal
% precoding for 1 user 

% SIC-based_hybrid_precoding using SVD by using sub-connected structure at 
% the Tx side and fully-connected structure at Rx sides  (Method 4)
%function [a_bs,D_bs,A_ms,D_ms,H_eff]=SIC_based_precoding(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
function [P,C,H_eff]=SIC_based_precoding_fully_connected(Num_users,Nt_RF,Nr_RF,Nt,Nr,H)
channel=zeros(Nr, Nt);
Hu_eff=zeros(Nr_RF,Nt_RF,Num_users);
a_ms=zeros(Nr,Nr_RF,Num_users);
for user=1:1:Num_users
channel(:,:)= H(:,:,user);
[Uu,~,~]=svd(channel/sqrt(Nt));
a_ms(:,:,user)=exp(1i*angle(Uu(:,1:Nr_RF)))/sqrt(Nr);
end
a_ms_cell=num2cell(a_ms,[1,2,4]);
A_ms=blkdiag(a_ms_cell{:});
H_eff=permute(H,[1 3 2]);
H_eff=reshape(H_eff,[],size(H,2),1);

H_comp = A_ms'*H_eff;
[~,~,V_comp]=svd(H_comp/sqrt(Nt));
a_bs=exp(1i*angle(V_comp(:,1:Nt_RF)))/sqrt(Nt);

d_ms=zeros(Nr_RF,Nr_RF,Num_users);
d_bs=zeros(Nt_RF,Nr_RF,Num_users);

Vu_eff_zero=zeros(Nt_RF,Nt_RF,Num_users);

for users=1:1:Num_users
    Hu_eff(:,:,users)=a_ms(:,:,users)'*H(:,:,users)*a_bs;
    [Uu_eff(:,:,users),~,Vu_eff(:,:,users)]=svd(Hu_eff(:,:,users));
    d_ms(:,:,users)=Uu_eff(:,1:Nr_RF,users);
    
    if users==1
        d_bs(:,:,users)=Vu_eff(:,1:Nr_RF,users); %correct
%        d_bs(:,:,users)=[]; %not correct
    else
        Vu_eff_zero(:,(users-1)*Nr_RF+1:end,users)=Vu_eff(:,(users-1)*Nr_RF+1:end,users-1);
        [~,~,Vu_tilde]=svd(d_ms(:,:,users)'*Hu_eff(:,:,users)*Vu_eff_zero(:,(Num_users-1)*Nr_RF+1:end,users));
        d_bs(:,1:Nr_RF,users)=Vu_eff_zero(:,(Num_users-1)*Nr_RF+1:end,users)*Vu_tilde(:,1:Nr_RF);
     % check whether, Vu_eff_zero(:,(Num_users or users -1)*Nr_RF+1:end,users))
    end
    
end
D_bs=permute(d_bs,[1 2 3]);
D_bs=reshape(D_bs,size(d_bs,1),[],1);

d_ms_cell=num2cell(d_ms,[1,2,4]);
D_ms=blkdiag(d_ms_cell{:});
    
P=a_bs*D_bs;
C=A_ms*D_ms;
end
%}
%-----------------------------------------------------------------------------




