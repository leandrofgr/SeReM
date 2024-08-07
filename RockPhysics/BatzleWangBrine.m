function [K_fl, rho_fl] = BatzleWangBrine(temperature, pressure, salinity)
% Batzle and Wang model to compute Brine/water properties
% Mavko, G., Mukerji, T. and Dvorkin, J., 2009. The rock physics handbook. 2nd Edition. Cambridge university press. Page 343.           
% INPUT temperature = Temperature in C (ex. 90)
%       pressure = Pressure in MPa (ex. 60)
%       salinity = Salinity in PPM (ex. 0.08)
% OUTUPT K_fl = Oil bulk modulus in GPa
%        rho_oil = Oil density in g/cm^3

    rho_fl = density(temperature, pressure, salinity);
    K_fl = calculate_bulk_modulus(temperature, pressure, salinity);
    
end

function bulk_modulus = calculate_bulk_modulus(temp, pres, salin)
    density_value = density(temp, pres, salin);
    acoustic_velocity_value = acoustic_velocity(temp, pres, salin);
    bulk_modulus = density_value * acoustic_velocity_value^2 * 1e-6; % in GPa
end

function density_value = density(temp, pres, salin)
    pure_water_density_value = pure_water_density(temp, pres);
    density_value = pure_water_density_value + salin * ( ...
        0.668 + 0.44 * salin + (10^-6) * (300 * pres - 2400 * pres * salin + ...
        temp * (80 + 3 * temp - 3300 * salin - 13 * pres + 47 * pres * salin)) );
end

function velocity = acoustic_velocity(temp, pres, salin)
    velocity = pure_water_acoustic_velocity(temp, pres) + salin * ( ...
        1170 - 9.6 * temp + 0.055 * (temp^2) - 8.5 * (10^-5) * (temp^3) + ...
        2.6 * pres - 0.0029 * temp * pres - 0.0476 * (pres^2) ) + ...
        (salin^(3/2)) * (780 - 10 * pres + 0.16 * (pres^2)) - 1820 * (salin^2);
end

function pure_water_density_value = pure_water_density(temp, pres)
    pure_water_density_value = 1 + (10^-6) * (-80 * temp - 3.3 * (temp^2) + 0.00175 * (temp^3) + ...
        489 * pres - 2 * temp * pres + 0.016 * (temp^2) * pres - 1.3 * (10^-5) * (temp^3) * pres - ...
        0.333 * (pres^2) - 0.002 * temp * (pres^2));
end

function velocity = pure_water_acoustic_velocity(temp, pres)
    OMEGA_COEFFICIENTS = [ ...
        1402.85, 1.524, 3.437e-3, -1.197e-5; ...
        4.871, -0.0111, 1.739e-4, -1.628e-6; ...
        -0.04783, 2.747e-4, -2.135e-6, 1.237e-8; ...
        1.487e-4, -6.503e-7, -1.455e-8, 1.327e-10; ...
        -2.197e-7, 7.987e-10, 5.230e-11, -4.614e-13];
    
    velocity = 0;
    for i = 1:5
        for j = 1:4
            velocity = velocity + OMEGA_COEFFICIENTS(i, j) * (temp^(i-1)) * (pres^(j-1));
        end
    end
end
