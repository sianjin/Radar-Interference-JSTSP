Tg = 55; %us
M = 2; % 2 bands
d = 1000; %m
c = 300; % m/us
fH = 10; % MHz
h= 12; %MHz/us
L = 200;
N = 450;
X = 10; % number of targets
%% for UARFDM
r0_ua = 435;  
pn = 0.221;
% 3.12MHz -> n0 =140 -> r0 =60, pn = 0.221
% 5.12MHz -> n0 =230 -> r0 =170, pn = 0.221
% 8.12MHz -> n0 =365 -> r0 =435, pn = 0.221
%% for UAFH
r0fh = 265;  
% 3.12MHz -> n0 =140 -> r0 =40
% 5.12MHz -> n0 =230 -> r0 =100
% 8.12MHz -> n0 =365 -> r0 =265
%% for UARFDMPC
r0_pc = 265;
% 3.12MHz -> n0 =140 -> r0 =40
% 5.12MHz -> n0 =230 -> r0 =130
% 8.12MHz -> n0 =365 -> r0 =265
%% for UAFHPC
r0_fhpc = 200;
% 3.12MHz -> n0 =140 -> r0 =25
% 5.12MHz -> n0 =230 -> r0 =85
% 8.12MHz -> n0 =365 -> r0 =200
%% closed-form
PmdUARFDM = [];
PmdUAFH = [];
PmdUARFDMPC = [];
PmdUAFHPC = [];
maxDensity = 10;
for lambda = 0:0.5:maxDensity
    PmdUA_temp = 1- exp(-2*r0_ua*lambda*fH/h/Tg/N/M*(1-(1-1/N * pn)^X));
    PmdUARFDM = [PmdUARFDM, PmdUA_temp];
end
for lambda = 0:0.5:maxDensity
    PmdUAFH_temp = 1- exp(-2*r0fh*lambda*fH/h/Tg/N*(1-(1-1/N)^X));
    PmdUAFH = [PmdUAFH, PmdUAFH_temp];
end
for lambda = 0:0.5:maxDensity
    PmdUAPC_temp = 1 -exp(-2*r0_pc*lambda*fH/h/Tg/N/M*(1-(1-1/N)^X));
    PmdUARFDMPC = [PmdUARFDMPC, PmdUAPC_temp];
end
for lambda = 0:0.5:maxDensity
    PmdUAFHPC_temp = 1 -exp(-2*r0_fhpc*lambda*fH/h/Tg/N*(1-(1-1/N)^X));
    PmdUAFHPC = [PmdUAFHPC, PmdUAFHPC_temp];
end
plot([0:0.5:maxDensity], PmdUARFDM, '.',[0:0.5:maxDensity], PmdUAFH, [0:0.5:maxDensity], PmdUARFDMPC, '--',[0:0.5:maxDensity],PmdUAFHPC,'.-')
hold on
densityChatMax = Bc/(fH + d/c*h)*M/(2*d);
nonIntDensity = 0:0.0001:densityChatMax;
density = [densityChatMax+0.0001:0.0001:maxDensity];
midDetecChat = [nonIntDensity*0,1-densityChatMax./density];
plot([0:0.0001:maxDensity],midDetecChat)