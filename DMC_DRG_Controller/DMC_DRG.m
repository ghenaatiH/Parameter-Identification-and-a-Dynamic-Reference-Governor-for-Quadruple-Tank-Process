%% Scenario 1 in Section 6.1
clc;
clear;
close all;

%% Global vars
global gama1 gama2 K1 K2 A1 A2 A3 A4 a1 a2 a3 a4 g V1 V2 Vmax Vmin Y1max Y2max Y3max Y4max;

D13=18;d13=1.2;
D24=15;d24=1;
gama1=0.4;gama2=0.3;
K1=13;K2=13;
A1=pi*(0.5*D13)^2;A3=pi*(0.5*D13)^2;
A2=pi*(0.5*D24)^2;A4=pi*(0.5*D24)^2;
a1=pi*(0.5*d13)^2;
a3=pi*(0.5*d13)^2;
a2=pi*(0.5*d24)^2;
a4=pi*(0.5*d24)^2;
g=981;
Vmax=12;
Vmin=0;
%% Init
load('StepResponse.mat');
Y11=Y11(2:end);Y12=Y12(2:end);
Y21=Y21(2:end);Y22=Y22(2:end);
N = 20; % number of samplint instances to construct DMC model
p = 5; % Prediction Horizon
m = p; % Control Horizon
G11=zeros(p,m);G12=zeros(p,m);G21=zeros(p,m);G22=zeros(p,m);G=zeros(2*p,2*m);
for i=1:m
    G11(: , i) = [zeros(i-1 , 1); Y11(1:p-i+1)];
    G12(: , i) = [zeros(i-1 , 1); Y12(1:p-i+1)];
    G21(: , i) = [zeros(i-1 , 1); Y21(1:p-i+1)];
    G22(: , i) = [zeros(i-1 , 1); Y22(1:p-i+1)];
end
G=[G11,G12;G21,G22];
H11=zeros(p,N);H12=zeros(p,N);H21=zeros(p,N);H22=zeros(p,N);H=zeros(2*p,2*N);
for i = 1:p
    for j = 1:N
        H11(i ,j) = Y11(i+j) - Y11(j);
        H12(i ,j) = Y12(i+j) - Y12(j);
        H21(i ,j) = Y21(i+j) - Y21(j);
        H22(i ,j) = Y22(i+j) - Y22(j);
    end
end
H=[H11,H12;H21,H22];
%% Main Loop
Nt = 300;
Y = zeros(Nt , 4);
U = zeros(Nt , 2);
Vol = zeros(Nt , 2);
ElapsedTime = zeros(Nt-1 , 1);
dU = zeros(2*m+p , 1);
dU_Rep = zeros(2*m+p , Nt);   % Future control increments Repository
Umin = zeros(2*N , 1);
h0=[5;5;2.5;2.5];
V1=0;V2=0;
Q1=20; R=1; Q2=1; Q3=20; % Coefficients of the Equation 27 in paper.
Fval = zeros(Nt , 1);
ExitFlag = zeros(Nt , 1);
CHI = zeros(Nt , p);
Sol_Init.dU=dU;Sol_Init.w1=0;Sol_Init.w2=0;Sol_Init.w3=0;Sol_Init.w4=0;Sol_Init.w5=0;Sol_Init.w6=0;Sol_Init.z1=0;Sol_Init.z2=0;Sol_Init.z3=0;Sol_Init.z4=0;Sol_Init.z5=0;Sol_Init.z6=0;Sol_Init.Z_ToTAL=0;Sol_Init.Chi=zeros(p,1);Sol_Init.DeltaChi=zeros(p,1);Sol_Init.W1=zeros(p,1);Sol_Init.W2=zeros(p,1);Sol_Init.Ref1=zeros(p,1);Sol_Init.Ref2=zeros(p,1);Sol_Init.Y1max=0;Sol_Init.Y2max=0;Sol_Init.Y3max=0;Sol_Init.Y4max=0;Sol_Init.Y1_Pred=zeros(p,1);Sol_Init.Y2_Pred=zeros(p,1);
Sol=repmat(Sol_Init,Nt , 1);
RefTrackFlg=1;
Refrences=[8.14, 5.20, 8.14, 4.60;...
    11.31, 7.20, 11.31, 4.50];
W1 = [Refrences(1,1)*ones(ceil((Nt)/4) , 1); Refrences(1,2)*ones(ceil((Nt)/4) , 1); Refrences(1,3)*ones(ceil((Nt)/4) , 1); Refrences(1,4)*ones(ceil((Nt)/4)+p , 1)];
W2 = [Refrences(2,1)*ones(ceil((Nt)/4) , 1); Refrences(2,2)*ones(ceil((Nt)/4) , 1); Refrences(2,3)*ones(ceil((Nt)/4) , 1); Refrences(2,4)*ones(ceil((Nt)/4)+p , 1)];
W = [W1,W2];
for i = 1:Nt-1
    for j = 1:N
        if i-j < 1
            Umin(j ,1) = 0;
            Umin(j+N ,1) = 0;
        else
            if i-j-1>0
                Umin(j , 1) = U(i-j,1)-U(i-j-1,1);
                Umin(j+N , 1) = U(i-j,2)-U(i-j-1,2);
            else
                Umin(j , 1) = U(i-j,1);
                Umin(j+N , 1) = U(i-j,2);
            end
        end
    end
    if i>1
        FreeResponse = H *Umin + [Y(i-1,1)*ones(p,1);Y(i-1,2)*ones(p,1)];
    else
        FreeResponse = H *Umin + [h0(1)*ones(p,1);h0(2)*ones(p,1)];
    end
    tic;
    % Dynamic matrix for converting voltage increments to absolute voltages
    S = tril(ones(m));
    % Coefficient matrix for the two pumps
    S2 = [blkdiag(S,S),zeros(2*size (S,1),size (S,2))];
    % Coefficient matrix for the two outputs
    G2 = [G,zeros(2*p,p)];
    
    % Previous voltages of the two pumps
    if i==1
        Vprev = [V1; V2];
    else
        Vprev = [U(i-1,1); U(i-1,2)];
    end
    
    % Voltage lower and upper bounds
    Vmin = 0;
    Vmax = 12;
    Y1max= 6.79;
    Y2max= 9;   
    Y3max= 6; 
    Y4max= 6; 
    % Construct linear inequality constraints:
    % Vmin <= Vprev + S2*dU <= Vmax
    % G dU + FreeResponse <= Ymax   
    A = [S2;
        -S2;
         G2]; 
    b = [Vmax*ones(2*m,1) - kron(Vprev,ones(m,1));
        -Vmin*ones(2*m,1) + kron(Vprev,ones(m,1));
         kron([Y1max;Y2max],ones(m,1))-FreeResponse];  
    
    Aeq = [];
    beq = [];
    LB = [-3*ones(2*m,1);zeros(p,1)];
    UB = [+3*ones(2*m,1);ones(p,1)];
     
    if i==1
        DU_Init = abs([Vmax*rand(m,1)-Vmax;Vmax*rand(m,1)-Vmax;rand(p,1)]);
    else
        %DU_Init = abs([Vmax*rand(m,1)-Vmax;Vmax*rand(m,1)-Vmax;rand(p,1)]);
        %DU_Init = abs(randn(size(dU)));
        DU_Init=dU;
    end
    options = optimoptions('fmincon','Display','final','Algorithm','SQP','MaxFunctionEvaluations',20000);
    [dU , Fval(i), ExitFlag(i)] = fmincon(@(dU) DMC_DRG_Cost(dU, G, FreeResponse, h0, W, RefTrackFlg, i, p, m, Q1, R, Q2, Q3), DU_Init , A , b , Aeq , beq , LB , UB ,@(dU)NLC_Scenario(dU,U,i,p,m,h0, Ts),options);
    dU_Rep(:,i)=dU;
    [~, Sol(i)]=DMC_DRG_Cost(dU, G, FreeResponse, h0, W, RefTrackFlg, i, p, m, Q1, R, Q2, Q3);   
    CHI(i,:)=dU(2*m+1:2*m+p);
    ElapsedTime(i)=toc;

    if i>1
        U(i,1) = U(i-1,1) + dU(1);
        U(i,2) = U(i-1,2) + dU(1+m);
    else
        U(i,1) = V1 + dU(1);
        U(i,2) = V2 + dU(1+m);
    end
%     V1=min(max(Vmin,U(i,1)),Vmax);
%     V2=min(max(Vmin,U(i,2)),Vmax);
    V1=U(i,1);
    V2=U(i,2);
    Vol(i,1) = V1;
    Vol(i,2) = V2;
    %----------------------------------------------------------------------
    [T_Sim,H_Sim] = ode45(@QTP,Ts*[i-1,i],h0);
    Y(i,1)=H_Sim(end,1); Y(i,2)=H_Sim(end,2); Y(i,3)=H_Sim(end,3);Y(i,4)=H_Sim(end,4);
    h0=H_Sim(end,:);
    disp([' TimeInstance = ', num2str(i)]);
end
%% plot Results
figure(1);
subplot(2 , 2 ,1);
TSIM=Ts*(1:Nt-1)';
plot(TSIM , Y(1:end-1,1) , TSIM , W(1:Nt-1,1),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output1' , 'Ref1');
grid on
subplot(2 , 2 ,2);
plot(TSIM , Y(1:end-1,2) , TSIM , W(1:Nt-1,2),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output2' , 'Ref2');
grid on
subplot(2, 2, 3);
plot(TSIM , U(1:end-1,1) ,TSIM , U(1:end-1,2), 'LineWidth' , 1.5); hold on
xlabel('Time (second)');
ylabel('U amp');
legend('U1','U2');
title('Input Signal (Effort Control)');
grid on
subplot(2, 2, 4);
plot(TSIM , Vol(1:end-1,1) ,TSIM , Vol(1:end-1,2), 'LineWidth' , 1.5); hold on
xlabel('Time (second)');
ylabel('U amp');
legend('Vol1','Vol2');
title('Input Signal (Effort Control)');
grid on

figure(2);
subplot(2 , 2 ,1);
TSIM=Ts*(1:Nt-1)';
plot(TSIM , Y(1:end-1,1) , TSIM , W(1:Nt-1,1),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output1' , 'Ref1');
grid on
subplot(2 , 2 ,2);
plot(TSIM , Y(1:end-1,2) , TSIM , W(1:Nt-1,2),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output2' , 'Ref2');
grid on
subplot(2, 2, 3);
plot(TSIM , Y(1:end-1,3),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Tank 3');
grid on
subplot(2, 2, 4);
plot(TSIM , Y(1:end-1,4),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Tank 4');
grid on

ET=sum(ElapsedTime)
FV=sum(Fval)
save(sprintf('Results_DMC_DRG_Scenario1_%d.mat', randi([10, 99])));

%% Scenario 1 in Section 6.2
clc;
clear;
close all;

%% Global vars
global gama1 gama2 K1 K2 A1 A2 A3 A4 a1 a2 a3 a4 g V1 V2 Vmax Vmin Y1max Y2max Y3max Y4max;

D13=18;d13=1.2;
D24=15;d24=1;
gama1=0.4;gama2=0.3;
K1=13;K2=13;
A1=pi*(0.5*D13)^2;A3=pi*(0.5*D13)^2;
A2=pi*(0.5*D24)^2;A4=pi*(0.5*D24)^2;
a1=pi*(0.5*d13)^2;
a3=pi*(0.5*d13)^2;
a2=pi*(0.5*d24)^2;
a4=pi*(0.5*d24)^2;
g=981;
Vmax=12;
Vmin=0;
%% Init
load('StepResponse.mat');
Y11=Y11(2:end);Y12=Y12(2:end);
Y21=Y21(2:end);Y22=Y22(2:end);
N = 20; % number of samplint instances to construct DMC model
p = 5; % Prediction Horizon
m = p; % Control Horizon
G11=zeros(p,m);G12=zeros(p,m);G21=zeros(p,m);G22=zeros(p,m);G=zeros(2*p,2*m);
for i=1:m
    G11(: , i) = [zeros(i-1 , 1); Y11(1:p-i+1)];
    G12(: , i) = [zeros(i-1 , 1); Y12(1:p-i+1)];
    G21(: , i) = [zeros(i-1 , 1); Y21(1:p-i+1)];
    G22(: , i) = [zeros(i-1 , 1); Y22(1:p-i+1)];
end
G=[G11,G12;G21,G22];
H11=zeros(p,N);H12=zeros(p,N);H21=zeros(p,N);H22=zeros(p,N);H=zeros(2*p,2*N);
for i = 1:p
    for j = 1:N
        H11(i ,j) = Y11(i+j) - Y11(j);
        H12(i ,j) = Y12(i+j) - Y12(j);
        H21(i ,j) = Y21(i+j) - Y21(j);
        H22(i ,j) = Y22(i+j) - Y22(j);
    end
end
H=[H11,H12;H21,H22];
%% Main Loop
Nt = 300;
Y = zeros(Nt , 4);
U = zeros(Nt , 2);
Vol = zeros(Nt , 2);
ElapsedTime = zeros(Nt-1 , 1);
dU = zeros(2*m+p , 1);
dU_Rep = zeros(2*m+p , Nt);   % Future control increments Repository
Umin = zeros(2*N , 1);
h0=[5;5;2.5;2.5];
V1=0;V2=0;
Q1=20; R=1; Q2=1; Q3=20; % Coefficients of the Equation 27 in paper.
Fval = zeros(Nt , 1);
ExitFlag = zeros(Nt , 1);
CHI = zeros(Nt , p);
Sol_Init.dU=dU;Sol_Init.w1=0;Sol_Init.w2=0;Sol_Init.w3=0;Sol_Init.w4=0;Sol_Init.w5=0;Sol_Init.w6=0;Sol_Init.z1=0;Sol_Init.z2=0;Sol_Init.z3=0;Sol_Init.z4=0;Sol_Init.z5=0;Sol_Init.z6=0;Sol_Init.Z_ToTAL=0;Sol_Init.Chi=zeros(p,1);Sol_Init.DeltaChi=zeros(p,1);Sol_Init.W1=zeros(p,1);Sol_Init.W2=zeros(p,1);Sol_Init.Ref1=zeros(p,1);Sol_Init.Ref2=zeros(p,1);Sol_Init.Y1max=0;Sol_Init.Y2max=0;Sol_Init.Y3max=0;Sol_Init.Y4max=0;Sol_Init.Y1_Pred=zeros(p,1);Sol_Init.Y2_Pred=zeros(p,1);
Sol=repmat(Sol_Init,Nt , 1);
RefTrackFlg=1;
Refrences=[8.14, 5.20, 8.14, 4.60;...
    11.31, 7.20, 11.31, 4.50];
W1 = [Refrences(1,1)*ones(ceil((Nt)/4) , 1); Refrences(1,2)*ones(ceil((Nt)/4) , 1); Refrences(1,3)*ones(ceil((Nt)/4) , 1); Refrences(1,4)*ones(ceil((Nt)/4)+p , 1)];
W2 = [Refrences(2,1)*ones(ceil((Nt)/4) , 1); Refrences(2,2)*ones(ceil((Nt)/4) , 1); Refrences(2,3)*ones(ceil((Nt)/4) , 1); Refrences(2,4)*ones(ceil((Nt)/4)+p , 1)];
W = [W1,W2];
for i = 1:Nt-1
    for j = 1:N
        if i-j < 1
            Umin(j ,1) = 0;
            Umin(j+N ,1) = 0;
        else
            if i-j-1>0
                Umin(j , 1) = U(i-j,1)-U(i-j-1,1);
                Umin(j+N , 1) = U(i-j,2)-U(i-j-1,2);
            else
                Umin(j , 1) = U(i-j,1);
                Umin(j+N , 1) = U(i-j,2);
            end
        end
    end
    if i>1
        FreeResponse = H *Umin + [Y(i-1,1)*ones(p,1);Y(i-1,2)*ones(p,1)];
    else
        FreeResponse = H *Umin + [h0(1)*ones(p,1);h0(2)*ones(p,1)];
    end
    tic;
    % Dynamic matrix for converting voltage increments to absolute voltages
    S = tril(ones(m));
    % Coefficient matrix for the two pumps
    S2 = [blkdiag(S,S),zeros(2*size (S,1),size (S,2))];
    % Coefficient matrix for the two outputs
    G2 = [G,zeros(2*p,p)];
    
    % Previous voltages of the two pumps
    if i==1
        Vprev = [V1; V2];
    else
        Vprev = [U(i-1,1); U(i-1,2)];
    end
    
    % Voltage lower and upper bounds
    Vmin = 0;
    Vmax = 12;
    Y1max= 6.79;
    Y2max= 9;   
    Y3max= 3; 
    Y4max= 3; 
    % Construct linear inequality constraints:
    % Vmin <= Vprev + S2*dU <= Vmax
    % G dU + FreeResponse <= Ymax   
    A = [S2;
        -S2;
         G2]; 
    b = [Vmax*ones(2*m,1) - kron(Vprev,ones(m,1));
        -Vmin*ones(2*m,1) + kron(Vprev,ones(m,1));
         kron([Y1max;Y2max],ones(m,1))-FreeResponse];  
    
    Aeq = [];
    beq = [];
    LB = [-3*ones(2*m,1);zeros(p,1)];
    UB = [+3*ones(2*m,1);ones(p,1)];
     
    if i==1
        DU_Init = abs([Vmax*rand(m,1)-Vmax;Vmax*rand(m,1)-Vmax;rand(p,1)]);
    else
        %DU_Init = abs([Vmax*rand(m,1)-Vmax;Vmax*rand(m,1)-Vmax;rand(p,1)]);
        %DU_Init = abs(randn(size(dU)));
        DU_Init=dU;
    end
    options = optimoptions('fmincon','Display','final','Algorithm','SQP','MaxFunctionEvaluations',20000);
    [dU , Fval(i), ExitFlag(i)] = fmincon(@(dU) DMC_DRG_Cost(dU, G, FreeResponse, h0, W, RefTrackFlg, i, p, m, Q1, R, Q2, Q3), DU_Init , A , b , Aeq , beq , LB , UB ,@(dU)NLC_Scenario(dU,U,i,p,m,h0, Ts),options);
    dU_Rep(:,i)=dU;
    [~, Sol(i)]=DMC_DRG_Cost(dU, G, FreeResponse, h0, W, RefTrackFlg, i, p, m, Q1, R, Q2, Q3);   
    CHI(i,:)=dU(2*m+1:2*m+p);
    ElapsedTime(i)=toc;
    
    if i>1
        U(i,1) = U(i-1,1) + dU(1);
        U(i,2) = U(i-1,2) + dU(1+m);
    else
        U(i,1) = V1 + dU(1);
        U(i,2) = V2 + dU(1+m);
    end
%     V1=min(max(Vmin,U(i,1)),Vmax);
%     V2=min(max(Vmin,U(i,2)),Vmax);
    V1=U(i,1);
    V2=U(i,2);
    Vol(i,1) = V1;
    Vol(i,2) = V2;
    %----------------------------------------------------------------------
    [T_Sim,H_Sim] = ode45(@QTP,Ts*[i-1,i],h0);
    Y(i,1)=H_Sim(end,1); Y(i,2)=H_Sim(end,2); Y(i,3)=H_Sim(end,3);Y(i,4)=H_Sim(end,4);
    h0=H_Sim(end,:);
    disp([' TimeInstance = ', num2str(i)]);
end
%% plot Results
figure(1);
subplot(2 , 2 ,1);
TSIM=Ts*(1:Nt-1)';
plot(TSIM , Y(1:end-1,1) , TSIM , W(1:Nt-1,1),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output1' , 'Ref1');
grid on
subplot(2 , 2 ,2);
plot(TSIM , Y(1:end-1,2) , TSIM , W(1:Nt-1,2),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output2' , 'Ref2');
grid on
subplot(2, 2, 3);
plot(TSIM , U(1:end-1,1) ,TSIM , U(1:end-1,2), 'LineWidth' , 1.5); hold on
xlabel('Time (second)');
ylabel('U amp');
legend('U1','U2');
title('Input Signal (Effort Control)');
grid on
subplot(2, 2, 4);
plot(TSIM , Vol(1:end-1,1) ,TSIM , Vol(1:end-1,2), 'LineWidth' , 1.5); hold on
xlabel('Time (second)');
ylabel('U amp');
legend('Vol1','Vol2');
title('Input Signal (Effort Control)');
grid on

figure(2);
subplot(2 , 2 ,1);
TSIM=Ts*(1:Nt-1)';
plot(TSIM , Y(1:end-1,1) , TSIM , W(1:Nt-1,1),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output1' , 'Ref1');
grid on
subplot(2 , 2 ,2);
plot(TSIM , Y(1:end-1,2) , TSIM , W(1:Nt-1,2),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Output2' , 'Ref2');
grid on
subplot(2, 2, 3);
plot(TSIM , Y(1:end-1,3),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Tank 3');
grid on
subplot(2, 2, 4);
plot(TSIM , Y(1:end-1,4),  'LineWidth' , 2); hold on
xlabel('Time (second)');
ylabel('Y Amp');
title('DMC Algorithm');
legend('Tank 4');
grid on

ET=sum(ElapsedTime)
FV=sum(Fval)
save(sprintf('Results_DMC_DRG_Scenario2_%d.mat', randi([10, 99])));
