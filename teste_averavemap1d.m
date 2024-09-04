
clear
%close all

%% average map test
nm = 20;
std_prior = 0.1;
mean_prior = 0.2;

average_map = 0.1;
average_map_std = 0.01;

mprior = mean_prior*ones(nm,1);

dt = 1;
corrlength = 3*dt;
trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = abs(trow-tcol);
sigmatime = exp(-abs((tdis./corrlength)).^2);
%sigmatime = exp(-abs((tdis./corrlength)));
sigmaprior = sigmatime * std_prior^2;

sigmaerr = average_map_std^2;

G = ones(1,nm);
G = G/sum(G);

mdobs = G*mprior;
% covariance matrix
sigmadobs = G*sigmaprior*G'+sigmaerr;

% posterior mean
mpost = mprior+(G*sigmaprior)'*(sigmadobs\(average_map-mdobs));
% posterior covariance matrix
sigmapost = sigmaprior-(G*sigmaprior)'*(sigmadobs\(G*sigmaprior));


samples = mvnrnd(mpost,sigmapost,1000)';
samples_average_map = G*samples;

x = linspace(-0.5,1,200);
likelihood_pdf = normpdf(x,average_map,average_map_std);
prior_pdf = normpdf(x,mean_prior,std_prior);

figure(1)
clf
subplot(121)
imagesc(samples)
subplot(122)
plot(samples)

figure(2)
clf
subplot(131)
histogram(samples_average_map,'Normalization','pdf')
hold all
plot(x,likelihood_pdf)
plot(x,prior_pdf)
xlim([0 0.5])

subplot(132)
histogram(samples(:),'Normalization','pdf')
hold all
plot(x,prior_pdf)
xlim([-0.2 0.7])

subplot(133)
%samples_ = samples - mean(samples(:));
samples_ = samples - mpost;
corr_func = real(ifft( fft(samples_).*conj(fft(samples_)) ));
corr_func = mean(corr_func');
plot(corr_func(1:round(nm/2))/corr_func(1))
hold all
plot(sigmapost(:,1)/sigmapost(1,1))
plot(sigmaprior(:,1)/sigmaprior(1,1))





