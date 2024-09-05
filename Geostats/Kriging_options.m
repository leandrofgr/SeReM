function [x_est, x_var] = Kriging_options(xcoord, dcoords, dvalues, xvar, l, type, krig, angles)
%%  This function works only by limiting the number of max points for conditioning sorted by distande.
% TO DO: % Kriging_options is being implemented to account for additional options
% such as max number of conditioning points, searching neighborhood and type
% of kriging, including the for loop over all the points 
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       xvar = prior variance
%       l = correlation length, 1x1 for isotropic or 3x1 for anisotropic
%       type = function ype ('exp', 'gau', 'sph')
%       angles = angles for anisotropic variogram , 3x1  
% OUTPUT x_est = kriging estimate
%        x_var = kriging variance

x_est = zeros(size(xcoord,1),1);
x_var = zeros(size(xcoord,1),1);

nmax = 30;
xmean = mean(dvalues);
if numel(dvalues) > nmax
    
    parfor i=1:size(xcoord,1)        
        dist = sqrt(sum((xcoord(i,:) - dcoords).^2,2));
        [~, sorted_indices] = sort(dist);
        lowest_indices = sorted_indices(1:nmax);
        if krig == 0
            [x_est(i,:), x_var(i,:) ] = SimpleKriging(xcoord(i,:), dcoords(lowest_indices,:), dvalues(lowest_indices,:), xmean, xvar, l, type, angles);
        else
            [x_est(i,:), x_var(i,:) ] = OrdinaryKriging(xcoord(i,:), dcoords(lowest_indices,:), dvalues(lowest_indices,:), xvar, l, type, angles);
        end
    end
    
else
    
    parfor i=1:size(xcoord,1)
        if krig == 0
            [x_est(i,:), x_var(i,:) ] = SimpleKriging(xcoord(i,:), dcoords, dvalues, xmean, xvar, l, type, angles);
        else
            [x_est(i,:), x_var(i,:) ] = OrdinaryKriging(xcoord(i,:), dcoords, dvalues, xvar, l, type, angles);
        end
    end
    
end


end