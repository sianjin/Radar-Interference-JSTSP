%% range setup
c = 300; % m/us
h = 12; % slope: 12 MHz/us
Bc = 540; % chirp bandwidth: 540MHz
Tc = Bc/h; % 45 us, chirp duration
fH = 10; % LPF cutoff: 10MHz    
Fs_r = fH; % complex sample frequency on bin bins
N = Tc/ (1/fH); % number of range bins = Tc/ (1/fH)
t_r = (0: N-1)* (1/fH); % fast time vector, us
fr1 = 8.11; % IF frequency: 8.11MHz, 5.11MHz, 3.11MHz, ghost
fr2 = 8.12; % IF frequency: 8.12MHz, 5.12MHz, 3.12MHz, target   or  7.05
nr1 = fr1/fH*N; % range bin of ghost 
nr2 = fr2/fH*N; % range bin of target
dres = c/2/Bc; % range resolution
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
lD1 = (fD1+fDmax)/fDmax*L/2; % Doppler bin of ghost
lD2 = (fD2+fDmax)/fDmax*L/2; % Doppler bin of target
%% power setup
Pt_dB = 12; % dBm
G_t = 12; % dBi
G_r = 12; % dBi
RCS = 20; %m^2
dmax = 1000; % m
d_g = 265; %m 30, 10, 1000
d_t = fr2 * c /2/h; % target distance, m
SignalPower0_dBm = Pt_dB + G_t + G_r;
SignalPower0 = 10^(SignalPower0_dBm/10); 
ghostPowerLinear = SignalPower0 * (c*10^6)^2/(4*pi)^2/(fc*10^6)^2/ (d_g)^2; % ghost power, mw
ghostAmpLinear = sqrt(ghostPowerLinear);
targetPowerLinear_ref = SignalPower0 * RCS * (c*10^6)^2/(4*pi)^3/(fc*10^6)^2;
targetPowerLinear = targetPowerLinear_ref/ (d_t)^4; % target power, mW 
targetAmpLinear = sqrt(targetPowerLinear);
%% unslotted ALOHA with frequency hopping 
% example: 1 ghost and 1 target
M = 2; % number of subbands
FH1 = binornd(1,1.0/M,1,L);
sumFH1 = sum(FH1);
X_2D = [];
for l = 1: L
    X_2D = [X_2D; ghostAmpLinear*exp(i*2*pi*fr1*t_r)*exp(i*2*pi*fD1*t_l(l)/1000)*FH1(l) + targetAmpLinear*exp(i*2*pi*fr2*t_r)*exp(i*2*pi*fD2*t_l(l)/1000)];
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
subplot(3,1,1)
plot(range,PrangedB) 
xlabel('Range (m)')
ylabel('Power (dBm)')
Y_r_2D = Y_r_2D/N;
Power_Y_r_2D_linear = abs(Y_r_2D).^2; 
% ------------------------------------------------------------------------
% threshold detector for range detection
Power_Y_r_2D_dB = 10*log10(Power_Y_r_2D_linear); % calculate power for threshold detector
countY_r_2D_dB = Power_Y_r_2D_dB > -128; % threshold detector
% ------------------------------------------------------------------------
% cfar detector for range detection
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
%% detector for range-velocity detection
% ------------------------------------------------------------------------
Amp_Y_rD_2D_linear = abs(Y_rD_2D);
Pow_Y_rD_2D_linear = Amp_Y_rD_2D_linear.^2;
Pow_Y_rD_2D_dB = 10*log10(Pow_Y_rD_2D_linear);
countY_rD_2D = zeros(L,N);
noisePower = []; 
% targetThrePowerdB = getRefTargetThrePowerdB(targetPowerLinear_ref, dres, N, L);
% ------------------------------------------------------------------------
% MATLAB 1D-CFAR 
% detector = phased.CFARDetector('NumTrainingCells',10,'NumGuardCells',2, ...
%     'ProbabilityFalseAlarm',0.005,'OutputFormat','Detection index');
% cutidx = 1:L;
% for n = 1: N
%     dets_n = detector(Pow_Y_rD_2D_linear(:,n),cutidx);
%     if length(dets_n)>0
%        countY_rD_2D(dets_n, n) = 1;
%     else 
%         continue
%     end
% end
% countY_rD_2D = countY_rD_2D.* (Pow_Y_rD_2D_dB>-128);
% ------------------------------------------------------------------------
% my cfar detector & detect-and-classify
% 1D-CFAR
for n = 1: N
    % targetThrePowerdB_n = targetThrePowerdB(:,n);
    Pow_Y_rD_2D_linear_n = Pow_Y_rD_2D_linear(:,n);
    [Detect_n, noisePower_n]= cfar_ca1D(Pow_Y_rD_2D_linear_n,20,2,10,0); % cfar detector
    noisePower =[noisePower, noisePower_n'];
    if size(Detect_n,2)>0
        Detect_n_index = Detect_n(1,:);
        countY_rD_2D(Detect_n_index, n) = 1;
        %countY_rD_2D(:,n) = countY_rD_2D(:,n).* (Pow_Y_rD_2D_dB_n<targetThrePowerdB_n);
    else 
        continue
    end
end
countY_rD_2D = countY_rD_2D.* (Pow_Y_rD_2D_dB>-128);
noisePower = noisePower.* (Pow_Y_rD_2D_dB>-128);
% ------------------------------------------------------------------------
% MATLAB 2D-CFAR -> poor performance
% detector = phased.CFARDetector2D('TrainingBandSize',[5,5], ...
%     'ThresholdFactor','Auto','GuardBandSize',[5,5], ...
%     'ProbabilityFalseAlarm',5e-4,'Method','SOCA','ThresholdOutputPort',true);
% Ngc = detector.GuardBandSize(2);
% Ngr = detector.GuardBandSize(1);
% Ntc = detector.TrainingBandSize(2);
% Ntr = detector.TrainingBandSize(1);
% cutidx = [];
% colstart = Ntc + Ngc + 1;
% colend = N - ( Ntc + Ngc);
% rowstart = Ntr + Ngr + 1;
% rowend = L - ( Ntr + Ngr);
% for m = colstart:colend
%     for n = rowstart:rowend
%         cutidx = [cutidx,[n;m]];
%     end
% end
% dets = detector(Pow_Y_rD_2D_linear,cutidx);
% for k = 1 : length(dets)
%     countY_rD_2D(cutidx(1,k), cutidx(2,k)) = dets(k);
% end
% countY_rD_2D = countY_rD_2D.* (Pow_Y_rD_2D_dB>-128);
% ------------------------------------------------------------------------
%% plot
f_D = Fs_D*(-(L/2):(L/2)-1)/L *10^3; % kHz
velocity = f_D * 1000 * c/fc/2;
subplot(3,1,2)
[Range, Velocity] = meshgrid(range, velocity);
mesh(Range, Velocity, 10*log10(abs(Y_rD_2D).^2)) % 2D-FFT result
colorbar
colormap(hot)
cmap = colormap;
cmapNew = flip(cmap);
colormap(cmapNew)
xlim([0, 130])
ylim([-20,20])
zlim([-128, -40]);
xlabel('Range (m)')
ylabel('Velocity (m/s)')
zlabel('Power (dBm)')
subplot(3,1,3)
mesh(Range, Velocity, countY_rD_2D) % cfar detection result
colorbar
xlim([0, 130])
ylim([-20,20])
zlim([0, 1]);
xlabel('Range (m)')
ylabel('Velocity (m/s)')
% conclusion: 
% 1. when target power is small (target is far) or ghost power is large (ghost is near), if target and ghost fall in the same range bin (except same range-Doppler bin), fa & misdetection
% 2. when target power is large (target is close) or ghost power is small (ghost is far), if target and ghost fall
% in the same range bin (except same range-Doppler bin), fa only
% 3. ghosts' side lobe can also influence but not destroy target's spectrum 