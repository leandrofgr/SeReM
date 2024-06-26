function sgsim = GaussianSimulation(xcoord, dcoords, dvalues, xmean, xvar, l, type, krig, angles)

% GAUSSIAN SIMULATION  generates a realization of the random variable 
% conditioned on the available measurements
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       xmean = prior mean
%       xvar = prior variance
%       h = distance
%       l = correlation length
%       type = function ype ('exp', 'gau', 'sph')
%       krig = kriging type (0=simple, 1=ordinary)
%       angles = angles for anisotropic variogram , 3x1  
% OUTPUT sgsim = realization

% Written by Dario Grana (August, 2020)

if krig == 0
    [krigmean, krigvar] = SimpleKriging(xcoord, dcoords, dvalues, xmean, xvar, l, type, angles);
else
    [krigmean, krigvar] = OrdinaryKriging(xcoord, dcoords, dvalues, xvar, l, type, angles);
end

% realization
sgsim = krigmean+sqrt(krigvar)*randn(1);
