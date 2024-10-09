function [logs_simulated_all] = DMS(I,J, range, type, angles, cell_size, ref_variables, cond_pos, cond_value, num_of_sims, cond_variables)
% DIRECT MULTIVARIATE SIMULATION:
% INPUT
%   I,J - Simulation grid size
%   range = correlation length, 1x1 for isotropic or 3x1 for anisotropic
%   type = function ype ('exp', 'gau', 'sph')
%   angles = angles for anisotropic variogram , 3x1  
%   cell_size = Cell size of joint distribution. The default value is 0.05 (5% of the difference of min and max values of each variable)
%   ref_variables = Reference variables that is used to define the non-parametric joing distribution
%   dcoords = coordinates of the measurements, hard data, well data points (ns, ndim)
%   dvalues = values of the measurements, hard data, well data points  (ns, 1)
%   num_of_sims = Number of simulations.
%   cond_variables = Conditioned variables for petrophysical inversion applications (Example: Vp, Vs, Rho to estimate Phi, Sw, V_clay)
% Obs: Aiming a better perfomance, this implementation uses FFTMA and kriging/PFS. The methodology also allows to use sgs instead

%The simplest application is the unconditioned multivariate non parametric simulation:
% [simulations] = DMS(50,50, 3, 'exp', [0 0 0], 0.05, ref_variables, [], [], 10, [])

% Different applications for petrophysical inversion of elastic properties cond_variables
 
% When running with n_sim=1, cond_pos = [] and cond_value = [], it considers PHI(g) = 0.5 to estimate the median model:
% [simulations] = DMS(50,50, 3, 'exp', [0 0 0], 0.05, ref_variables, [], [], 1, cond_variables)

% When running with n_sim~=1, cond_pos = [] and cond_value = [], it considers g as a geostatistica simulation unconditioned to the transformed hard data cond_value
% n_sim = 2;
% [simulations] = DMS(50,50, 3, 'exp', [0 0 0], 0.05, ref_variables, [], [], n_sim, cond_variables)

% When running with n_sim~=1, cond_pos ~= [] and cond_value ~= [], it considers g as a geostatistica simulation conditioned to the transformed hard data cond_value
% [simulations] = DMS(50,50, 3, 'exp', [0 0 0], 0.05, ref_variables, [], [], n_sim, cond_variables)

n_vars = size(ref_variables,2);
n_cond_vars = size(cond_variables,2);


krig_mean = zeros(n_vars, I, J);
krig_std = ones(n_vars, I, J);
%% if we have conditioning points, use kriging to condition FFTMA simulation
if size(cond_value, 1) > 0
    
    % Gaussian transforming the conditional points
    cond_value_uniform = nonParametric_to_uniform( cond_value, ref_variables , cell_size);
    cond_value_gaussian_2krig = norminv(cond_value_uniform);
    
    [X,Y] = meshgrid(1:J,1:I);
    xcoords = [ Y(:) X(:)];
    for i = n_cond_vars + 1:n_vars
        krig = 1;
        [mean_krig, var_krig] = Kriging_options(xcoords, cond_pos, cond_value_gaussian_2krig(:,i), 1, range, type, krig, angles);
        krig_mean(i,:) = mean_krig(:);
        krig_std(i,:) = sqrt(var_krig(:));
    end
end

% correlation function to use FFTMA in DMS. It is also possible to use SGS instead
[correlation_function] = construct_correlation_function(zeros(I,J), range, type, angles, 1);

%% DMS
logs_simulated_all = cell(num_of_sims,1);
for n = 1:num_of_sims
    
    % simulating Gaussian realizations using FFTMA
    simulations2D = zeros(n_vars, I, J);
    if num_of_sims ~= 1
        for i = n_cond_vars + 1:n_vars
            white_noise = randn(I,J);
            simulations_gaussian = FFT_MA_3D(correlation_function,white_noise);
            simulations_gaussian = make_it_gaussian(simulations_gaussian(:)); % It makes the simulation to be "Perfectly" Gaussian distributed, therefore, better non-parametric simulations
            simulations_gaussian = reshape(simulations_gaussian,I,J);
            simulations2D(i,:,:) = simulations_gaussian;
        end
    end
    
    simulations_gaussian = simulations2D(:,:)';
    simulations_conditioned = zeros(size(simulations_gaussian));
    
    for i = n_cond_vars + 1:n_vars
        mean_krig = reshape(krig_mean(i,:), I, J);
        var_krig = reshape(krig_std(i,:), I, J);
        simulation_gaussian = reshape(simulations_gaussian(:,i), I, J);
        simulation_gaussian = simulation_gaussian.*var_krig +  mean_krig;
        simulations_conditioned(:,i) = reshape(simulation_gaussian, I * J, 1);
    end
    
    % Transforming the conticioned Gaussian realizations to the non parametris pdf
    simulations_conditioned = normcdf(simulations_conditioned);
    simulations_conditioned_nonParametric = uniform_to_nonParametric( simulations_conditioned, ref_variables, cell_size, cond_variables);    
    
    result = zeros(size(simulations_conditioned_nonParametric,2), I, J);
    for i = 1:n_vars
        result(i,:,:) = reshape(simulations_conditioned_nonParametric(:,i), I, J);
    end
    logs_simulated_all{n} = result;
end


