%% Seismic inversion Driver %%
% In this script we apply the Bayesian linearized AVO inversion method
% (Buland and Omre, 2003) to predict the elastic properties (P- and S-wave
% velocity and density) from seismic data.
clear
close all
%% Available data and parameters
% Load data (seismic data and time)
addpath(genpath('../SeReM/'))
load Data/data3.mat

%% Initial parameters
% number of samples (elastic properties)
nm = size(Snear,1)+1;
% number of variables
nv = 3;
% reflection angles 
theta = [15, 30, 45];
ntheta = length(theta);
% time sampling
dt = TimeSeis(2)-TimeSeis(1);
% error variance
varerr = 10^-4;
sigmaerr = varerr*eye(ntheta*(nm-1));

%% Wavelet
% wavelet 
freq = 45;
ntw = 64;
[wavelet, tw] = RickerWavelet(freq, dt, ntw);

%% Prior model (filtered well logs)
Vp = repmat(Vp,5,1);
Vs = repmat(Vs,5,1);
Rho = repmat(Rho,5,1);

nfilt = 3;
cutofffr = 0.04;
[b, a] = butter(nfilt, cutofffr);
Vpprior = filtfilt(b, a, Vp);
Vsprior = filtfilt(b, a, Vs);
Rhoprior = filtfilt(b, a, Rho);
logVpprior = log(Vpprior);
logVsprior = log(Vsprior);
logRhoprior = log(Rhoprior);
mprior = [logVpprior; logVsprior; logRhoprior];

nm = size(Vp,1);

%% Spatial correlation matrix
corrlength = 5*dt;
trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = abs(trow-tcol);
%sigmatime = exp(-(tdis./corrlength).^2);
sigmatime = exp(-abs(tdis./corrlength));
sigma0 = cov([log(Vp),log(Vs),log(Rho)]);
sigmaprior = kron(sigma0, sigmatime);

%% Operator in freq - test
% Time:
theta = [15, 30, 45];
ntheta = length(theta);
A = AkiRichardsCoefficientsMatrix(4*ones(size(Vpprior)), 2.5*ones(size(Vpprior)), theta, nv);
a = AkiRichardsCoefficientsMatrix([4 4], [2.5 2.5], theta, nv)
D = DifferentialMatrix(nm,nv);
W = WaveletMatrix(wavelet, nm, ntheta);

% forward operator
m = log([Vp;Vs;Rho]);
m_ = [fft(log(Vp(1:end-1)')); fft(log(Vs(1:end-1)')); fft(log(Rho(1:end-1)'))];

G = W*A*D;
seismic_time = G*m;

%freq:
Fs = 1;
N = nm-1;
w_shifted = 2 * pi * (fftshift((0:N-1) - floor(N/2)) * (Fs / N));
w_shifted = w_shifted';
s_ = [fft(wavelet',N); fft(wavelet',N); fft(wavelet',N)] ;

[~,A_m,~] = svd(sigmatime);
[~,A_d,~] = svd(sigmaerr);
A_m = diag(A_m);
A_d = diag(A_d);
A_d = A_d(1);
mprior_= [fft((logVpprior(1:end-1)')); fft((logVsprior(1:end-1)')); fft((logRhoprior(1:end-1)'))];
Seis = [Snear; Smid; Sfar];
for k=1:length(w_shifted)
    G_ = [];
    for ang = 1:ntheta
       G_ = [G_; i.*w_shifted(k)*s_(ang,k)*a(ang,:)] ;
    end    
    
    if k==22
        stop=1
    end
    
    d(:,k) =  G_*m_(:,k);        
    
    if s_(1,k)<0.01
        mu_post(:,k) = mprior_(:,k);
    else
        C_m_(:,:,k) = sigma0*A_m(k);
        C_d_(:,:,k) = G_*C_m_(:,:,k)*G_' + A_d;
        mu_d_(:,k) = G_*mprior_(:,k);
        
        OPERATOR_ = (G_*C_m_(:,:,k))'*inv(C_d_(:,:,k));        
        mu_post(:,k) = mprior_(:,k) + OPERATOR_*( d(:,k) - mu_d_(:,k) );
    end
    
end
mu_post(:,1) = mprior_(:,1);

Vp_ref = real(ifft(m_(1,:)));
Vs_ref = real(ifft(m_(2,:)));
Rho_ref = real(ifft(m_(3,:)));
Vp_prior = real(ifft(mprior_(1,:)));
Vs_prior = real(ifft(mprior_(2,:)));
Rho_prior = real(ifft(mprior_(3,:)));
Vp_inv = real(ifft(mu_post(1,:)));
Vs_inv = real(ifft(mu_post(2,:)));
Rho_inv = real(ifft(mu_post(3,:)));



for ang = 1:ntheta
    d(ang,:) = real(ifft(d(ang,:)));
end
d = circshift(d,-length(wavelet)/2,2);
d = d';
d = d(:);

figure
subplot(311)
plot(seismic_time)
hold all
plot(d)

figure
subplot(311)
%plot(log(Vp))
plot(Vp_ref)
hold all
plot(Vp_prior)
plot(Vp_inv)
subplot(312)
%plot(log(Vs))
plot(Vs_ref)
hold all
plot(Vs_prior)
plot(Vs_inv)
subplot(313)
%plot(log(Rho))
plot(Rho_ref)
hold all
plot(Rho_prior)
plot(Rho_inv)


%% Seismic inversion
% Seis = [Snear; Smid; Sfar];
% [mmap] = SeismicInversion(Seis, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, sigmaerr, wavelet, theta, nv);
% 
% Vpmap = mmap(1:nm);
% Vsmap = mmap(nm+1:2*nm);
% Rhomap = mmap(2*nm+1:end);
