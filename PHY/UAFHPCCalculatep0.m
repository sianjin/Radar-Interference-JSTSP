%% range setup
c = 300; % m/us
h = 12; % slope: 12 MHz/us
Bc = 540; % chirp bandwidth: 540MHz
Tc = Bc/h; % 45 us, chirp duration
fH = 10; % LPF cutoff: 10MHz 
Fs_r = fH; % complex sample frequency on bin bins
N = Tc/ (1/fH); % number of range bins = Tc/ (1/fH)
t_r = (0: N-1)* (1/fH); % fast time vector, us
% >4 strong interference in the same range bin happens rarely. 
% calculate such prob at density 3
% average # int: 3m^{-1}*2000m=6000 
% average # int fall into range <=200m: 6000*0.2=1200
% such prob is 1- binocdf(4,1200,1/N) = 0.1320 
% 0.1320 is a small probability -> 
% So, the case under consideration is <=3 interference case for each range bin 
% When beam size is considered or the interference range is smaller, 
% this condition holds more
fr1 = 8.11; % IF frequency: 8.11MHz, 5.11MHz, 3.11MHz, ghost
fr2 = 8.12; % IF frequency: 8.12MHz, 5.12MHz, 3.12MHz, target   or  7.05
fr3 = 8.09;
fr4 = 8.13;
fr5 = 8.10;
nr1 = fr1/fH*N; % range bin of ghost 
nr2 = fr2/fH*N; % range bin of target
%% Doppler setup
fc = 77* 1000; % (us)^{-1}, 77k MHz 
Tg = 55;           % 55us, inter-chirp duration
Fs_D = 1/Tg;         % complex sample frequency on Doppler bins
L = 200;           % number of chirps per frame
t_l = (0:L-1)*Tg;    % slow time vector, us
% Doppler frequency fD = 2vf_c/c; 
% maximum Doppler: -pi <= 2 pi * fD * Tg <= pi   
% -1/2 * 1/Tg <= fD <= 1/2 * 1/Tg
% as Tg is at the order of tens of us, fD is at the order of kHz
fDmax = 1/2 * 1/Tg*1000; % max positive Doppler: kHz
fD1 = 4.21; % Doppler: 4.21kHz -> ghost 16.40m/s  
fD2 = 6.21; % Doppler: 6.21kHz -> target 12.10m/s
fD3 = 1.2;
fD4 = 0.2;
fD5 = -2.2;
lD1 = (fD1+fDmax)/fDmax*L/2; % Doppler bin of ghost
lD2 = (fD2+fDmax)/fDmax*L/2; % Doppler bin of target
%% power setup
Pt_dB = 12; % dBm
G_t = 12; % dBi
G_r = 12; % dBi
RCS = 20; %m^2
dmax = 1000; % m
d_g = 140; %m
d_t = fr2 * c /2/h; % target distance, m
SignalPower0_dBm = Pt_dB + G_t + G_r;
SignalPower0 = 10^(SignalPower0_dBm/10); 
ghostPowerLinear = SignalPower0 * (c*10^6)^2/(4*pi)^2/(fc*10^6)^2/ (d_g)^2; % ghost power, mw
ghostAmpLinear = sqrt(ghostPowerLinear);
targetPowerLinear = SignalPower0 * RCS * (c*10^6)^2/(4*pi)^3/(fc*10^6)^2 / (d_t)^4; % target power, mW 
targetAmpLinear = sqrt(targetPowerLinear);
%% unslotted ALOHA with frequency hopping and phase coding
%% for loop to calculate p0
p0Count = 0;
MonteNum = 1000;
for simuTimes = 1: MonteNum
    % example: 3 ghosts and 1 target
    PC1 = exp(i*pi*randi([0,1],1,L));
    PC3 = exp(i*pi*randi([0,1],1,L));
    PC4 = exp(i*pi*randi([0,1],1,L));
    PC5 = exp(i*pi*randi([0,1],1,L));
    M = 2;
    FH1 = binornd(1,1.0/M,1,L);
    FH3 = binornd(1,1.0/M,1,L);
    FH4 = binornd(1,1.0/M,1,L);
    FH5 = binornd(1,1.0/M,1,L);
    X_2D = [];
    for l = 1: L
        signal_l = targetAmpLinear*exp(i*2*pi*fr2*t_r)*exp(i*2*pi*fD2*t_l(l)/1000);
        signal_l = signal_l + ghostAmpLinear*exp(i*2*pi*fr1*t_r)*exp(i*2*pi*fD1*t_l(l)/1000)*FH1(l)*PC1(l);
        signal_l = signal_l + ghostAmpLinear*exp(i*2*pi*fr3*t_r)*exp(i*2*pi*fD3*t_l(l)/1000)*FH3(l)*PC3(l);
        signal_l = signal_l + ghostAmpLinear*exp(i*2*pi*fr4*t_r)*exp(i*2*pi*fD4*t_l(l)/1000)*FH4(l)*PC4(l);
        signal_l = signal_l + ghostAmpLinear*exp(i*2*pi*fr5*t_r)*exp(i*2*pi*fD5*t_l(l)/1000)*FH5(l)*PC5(l);
        X_2D = [X_2D; signal_l];
        % X_2D = [X_2D;  targetAmpLinear*exp(i*2*pi*fr2*t_r)*exp(i*2*pi*fD2*t_l(l)/1000)];
        %X_2D = [X_2D; ghostAmpLinear*exp(i*2*pi*fr1*t_r)*exp(i*2*pi*fD1*t_l(l)/1000)*FH1(l)];
    end
    %% Range FFT
    Y_r_2D  = [];
    for l = 1: L
        X_2D_l = X_2D(l,:);
        Y_r_2D = [Y_r_2D; fft(X_2D_l,N)];
    end
    Prange = abs(Y_r_2D(1,:)/(N)).^2;
    PrangedB = 10*log10(Prange);
    f_r = Fs_r*(0:N-1)/N; % MHz
    range = f_r * c/2/h;
    Y_r_2D = Y_r_2D/N;
    Power_Y_r_2D_linear = abs(Y_r_2D).^2;
    % ------------------------------------------------------------------------
    % threshold detector
    Power_Y_r_2D_dB = 10*log10(Power_Y_r_2D_linear); % calculate power for threshold detector
    countY_r_2D_dB = Power_Y_r_2D_dB > -128; % threshold detector
    % ------------------------------------------------------------------------
    % cfar detector
    % countY_r_2D_dB = zeros(L,N);
    % for l = 1: L
    %     Amp_Y_r_2D_linear_l = sqrt(Power_Y_r_2D_linear(l,:));
    %     Detect_l = cfar_ca1D(Amp_Y_r_2D_linear_l,10,2,5,0); % cfar detector
    %     countY_r_2D_dB(l, Detect_l) = 1;
    % end
    % ------------------------------------------------------------------------
    classifyY_r_2D = sum(countY_r_2D_dB); % classify, entry n means number of signals over L chirp cycles on range bin n
    %% Doppler FFT
    Y_rD_2D = [];
    for n = 1: N
        if classifyY_r_2D(n) == L
            Y_r_2D_n = Y_r_2D(:,n);
            Y_rD_2D = [Y_rD_2D, fftshift(fft(Y_r_2D_n,L))];
        else
            Y_rD_2D = [Y_rD_2D, zeros(L,1)];
        end
    end
    Y_rD_2D = Y_rD_2D/L;
    %% detector
    % ------------------------------------------------------------------------
    Amp_Y_rD_2D_linear = abs(Y_rD_2D);
    Pow_Y_rD_2D_linear = Amp_Y_rD_2D_linear.^2;
    Pow_Y_rD_2D_dB = 10*log10(Pow_Y_rD_2D_linear);
    Pow_Y_rD_2D_dB_nr1 = Pow_Y_rD_2D_dB(:,ceil(nr1));
    countY_rD_2D_nr1 = zeros(L,1);
    % ------------------------------------------------------------------------
    % my cfar detector
    % 1D-CFAR
    for n = ceil(nr1)
        Pow_Y_rD_2D_linear_n = Pow_Y_rD_2D_linear(:,n);
        [Detect_n, noisePower_n]= cfar_ca1D(Pow_Y_rD_2D_linear_n,20,2,10,0); % cfar detector
        if size(Detect_n,2)>0
            Detect_n_index = Detect_n(1,:);
            countY_rD_2D_nr1(Detect_n_index,1) = 1;
        else
            continue
        end
    end
    countY_rD_2D_nr1 = countY_rD_2D_nr1.* (Pow_Y_rD_2D_dB_nr1>-128);
    countY_rD_2D_nr1([ceil(lD2)-2: ceil(lD2)+2]) = 0; % exclude the peak caused by target
    if sum(countY_rD_2D_nr1) >0
        p0Count = p0Count + 1;
    end
end
p0 = p0Count/MonteNum;
