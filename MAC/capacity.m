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
p0 = 0.033; % p0 for UAPC
p0_2 = 0.03; % p0 for UAFHPC
lambdaUARFDM = 1/d/2*M*Tg/(fH/h);
CUARFDMmax = lambdaUARFDM/exp(1);
lambdaUAFH5 = 1/d/2*Tg/(fH/h)/(1-(1-1/N)^5);
CUAFH5max = lambdaUAFH5/exp(1);
lambdaUAFH10 = 1/d/2*Tg/(fH/h)/(1-(1-1/N)^10);
CUAFH10max = lambdaUAFH10/exp(1);
lambdaUARFMDPC = 1/d/2*M*Tg/(fH/h)/p0;
CUARFDMPCmax = lambdaUARFMDPC/exp(1);
lambdaUAFHPC5 = 1/d/2*Tg/(fH/h)/p0_2/(1-(1-1/N)^5);
CUAFHPC5max = lambdaUAFHPC5/exp(1);
lambdaUAFHPC10 = 1/d/2*Tg/(fH/h)/p0_2/(1-(1-1/N)^10);
CUAFHPC10max = lambdaUAFHPC10/exp(1);
CUARFDM = [];
CUAFH5 = []; % number of targets = 5
CUAFH10 = []; % number of targets = 10
CUARFDMPC = [];
CUAFHPC5 = []; % number of targets = 5
CUAFHPC10 = []; % number of targets = 10
maxDensity = 10;
for lambda = 0:0.001:maxDensity
    CUARFDM = [CUARFDM, lambda * exp(-2*d*lambda*fH/h/Tg/M)];
end
for lambda = 0:0.5:maxDensity
    CUAFH5 = [CUAFH5, lambda * exp(-2*d*lambda*fH/h/Tg*(1-(1-1/N)^5))]; 
    CUAFH10 = [CUAFH10, lambda * exp(-2*d*lambda*fH/h/Tg*(1-(1-1/N)^10))]; 
    CUARFDMPC = [CUARFDMPC, lambda * exp(-2*d*lambda*fH/h/Tg/M*p0)];
    CUAFHPC5 = [CUAFHPC5, lambda * exp(-2*d*lambda*fH/h/Tg*p0_2*(1-(1-1/N)^5))];
    CUAFHPC10 = [CUAFHPC10, lambda * exp(-2*d*lambda*fH/h/Tg*p0_2*(1-(1-1/N)^10))];
end
plot([0:0.001:maxDensity], CUARFDM, 'b', [0:0.5:maxDensity], CUAFH5, 'y',[0:0.5:maxDensity], CUAFH10, 'o',[0:0.5:maxDensity], CUARFDMPC, 'r', [0:0.5:maxDensity], CUAFHPC5, 'g', [0:0.5:maxDensity], CUAFHPC10, '-*')
hold on
CchatMax = Bc/(fH + d/c*h)*M/(2*d);
Cchat = 0:0.001:CchatMax;
Cchat = [Cchat,CchatMax*ones(1,10000)];
plot([0:0.001:maxDensity],Cchat(1:10001))

%plot([0:0.5:10], CUA, 'b', [0:0.5:10], CUAFH, 'y', [0:0.5:10], CUAPC, '-o', [0:0.5:10], CUAFHPC, '-*');