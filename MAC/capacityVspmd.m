Tg = 55; %us
M = 2; % 2 bands
Bc = 540; % chirp bandwidth: 540MHz
d = 1000; %m
c = 300; % m/us
fH = 10; % MHz
h= 12; %MHz/us
L = 200;
L0 = 15; % for UAFH
N = 450;
X = 10; % number of targets
p0 = 0.033; % p0 for UAPC
p0_2 = 0.03; % p0 for UAFHPC
%% for UARFDM
r0_ua = 170;  
pn = 0.221;
% 3.12MHz -> n0 =140 -> r0 =60, pn = 0.221
% 5.12MHz -> n0 =230 -> r0 =170, pn = 0.221
% 8.12MHz -> n0 =365 -> r0 =435, pn = 0.221
%% for UAFH
r0fh = 100;  
% 3.12MHz -> n0 =140 -> r0 =40
% 5.12MHz -> n0 =230 -> r0 =100
% 8.12MHz -> n0 =365 -> r0 =265
%% for UARFDMPC
r0_pc = 130;
% 3.12MHz -> n0 =140 -> r0 =40
% 5.12MHz -> n0 =230 -> r0 =130
% 8.12MHz -> n0 =365 -> r0 =265
%% for UAFHPC
r0_fhpc = 85;
% 3.12MHz -> n0 =140 -> r0 =25
% 5.12MHz -> n0 =230 -> r0 =85
% 8.12MHz -> n0 =365 -> r0 =200
%% closed-form mis-detection probability
lambdaPmdUARFDM = [];
lambdaPmdUAFH = [];
lambdaPmdUARFDMPC = [];
lambdaPmdUAFHPC = [];
maxMisDetProb = 0.4/10^3;
for prob = 0:0.4*10^(-6):maxMisDetProb
    lambdaPmdUA_temp = -log(1-prob)*N*Tg*M*h/fH/(1-(1-1/N * pn)^X)/r0_ua/2;
    lambdaPmdUARFDM = [lambdaPmdUARFDM, lambdaPmdUA_temp];
end
for prob = 0:4*10^(-7):maxMisDetProb
    lambdaPmdUAFH_temp = -log(1-prob)*N*Tg*h/fH/(1-(1-1/N)^X)/r0fh/2;
    lambdaPmdUAFH = [lambdaPmdUAFH, lambdaPmdUAFH_temp];
end
for prob = 0:4*10^(-7):maxMisDetProb
    lambdaPmdUAPC_temp = -log(1-prob)*N*Tg*M*h/fH/(1-(1-1/N)^X)/r0_pc/2;
    lambdaPmdUARFDMPC = [lambdaPmdUARFDMPC, lambdaPmdUAPC_temp];
end
for prob = 0:4*10^(-7):maxMisDetProb
    lambdaPmdUAFHPC_temp = -log(1-prob)*N*Tg*h/fH/(1-(1-1/N)^X)/r0_fhpc/2;
    lambdaPmdUAFHPC = [lambdaPmdUAFHPC, lambdaPmdUAFHPC_temp];
end
%% Capacity
CUARFDM = lambdaPmdUARFDM.* exp(-2*d*lambdaPmdUARFDM*fH/h/Tg/M);
CUAFH = lambdaPmdUAFH.* exp(-2*d*lambdaPmdUAFH*fH/h/Tg*(1-(1-1/N)^X));
CUARFDMPC = lambdaPmdUARFDMPC.* exp(-2*d*lambdaPmdUARFDMPC*fH/h/Tg/M*p0);
CUAFHPC = lambdaPmdUAFHPC.* exp(-2*d*lambdaPmdUAFHPC*fH/h/Tg*p0_2*(1-(1-1/N)^X));
%% Plot
plot([0:0.4*10^(-6):maxMisDetProb],CUARFDM,'.',[0:4*10^(-7):maxMisDetProb],CUAFH,'+',[0:4*10^(-7):maxMisDetProb],CUARFDMPC,'o',[0:4*10^(-7):maxMisDetProb],CUAFHPC,'--')