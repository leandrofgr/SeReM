function [K_oil, rho_oil] = BatzleWangOil(temperature, pressure, gas_oil_ratio, api, gas_gravity)
% Batzle and Wang model to compute oil properties
% Mavko, G., Mukerji, T. and Dvorkin, J., 2009. The rock physics handbook. 2nd Edition. Cambridge university press. Page 343.           
% INPUT temperature = Temperature in C (ex. 90)
%       pressure = Pressure in MPa (ex. 60)
%       gas_oil_ratio = Gas Oil Ratio (ex. 70)
%       api = Oil API (ex. 17)
%       gas_gravity = Gas Gravity (ex. 0.6)
% OUTUPT K_oil = Oil bulk modulus in GPa
%        rho_oil = Oil density in g/cm^3

    use_max_gor = isnan(gas_oil_ratio);
    gas_oil_ratio = max_gor(temperature, pressure, api, gas_gravity) * use_max_gor + gas_oil_ratio * ~use_max_gor;

    rho_oil = density(temperature, pressure, gas_oil_ratio, api, gas_gravity);
    K_oil = bulk_modulus(temperature, pressure, gas_oil_ratio, api, gas_gravity);
    
end

function bulk_modulus = bulk_modulus(temp, pres, gas_oil_ratio, api, gas_gravity)
    density_value = density(temp, pres, gas_oil_ratio, api, gas_gravity);
    acoustic_velocity_value = acoustic_velocity(temp, pres, gas_oil_ratio, api, gas_gravity);
    bulk_modulus = real(density_value * acoustic_velocity_value^2 * 1e-6); % in GPa
end

function density_value = density(temp, pres, gas_oil_ratio, api, gas_gravity)
    if is_live_oil(gas_oil_ratio)
        rho = true_density(temp, pres, gas_oil_ratio, api, gas_gravity);
    else
        rho = reference_density(api);
    end

    if (is_live_oil(gas_oil_ratio) && gas_gravity ~= 0 && pres < pb(gas_oil_ratio, api, gas_gravity, temp)) || ~is_live_oil(gas_oil_ratio)
        rhop = density_pressure_dependent(rho, pres);  % correct for pressure
        rho = real(rhop / (0.972 + 0.000381 * (temp + 17.8)^1.175));  % correct for temp
    end
    density_value = rho; % in g/cm^3
end

function velocity = acoustic_velocity(temp, pres, gas_oil_ratio, api, gas_gravity)
    if is_live_oil(gas_oil_ratio)
        d = pseudo_density(temp, pres, gas_oil_ratio, api, gas_gravity);
    else
        d = reference_density(api);
    end

    velocity = 2096 * (d / (2.6 - d))^(0.5) ...
        - 3.7 * temp + 4.64 * pres + 0.0115 * (4.12 * (1.08 / d - 1)^0.5 - 1) * temp * pres;
end

function max_gor_value = max_gor(temp, pres, api, gas_gravity)
    max_gor_value = 2.03 * gas_gravity * ((pres * exp((0.02878 * api) - 0.00377 * temp))^1.205);
end

function live = is_live_oil(gas_oil_ratio)
    live = gas_oil_ratio > 1e-10;
end

function density_value = density_pressure_dependent(density, pres)
    density_value = density ...
        + (0.00277 * pres - 1.71 * 10^-7 * pres^3) * (density - 1.15)^2 ...
        + 3.49 * 10^-4 * pres;
end

function b0_value = b0(temp, gas_oil_ratio, gas_gravity, api)
    b0_value = 0.972 + 0.000381 * (2.495 * gas_oil_ratio * (gas_gravity / reference_density(api))^0.5 ...
        + temp + 17.8)^1.175;
end

function pseudo_density_value = pseudo_density(temp, pres, gas_oil_ratio, api, gas_gravity)
    pseudo_density_value = (reference_density(api) / b0(temp, gas_oil_ratio, gas_gravity, api)) ...
        * (1 + 0.001 * gas_oil_ratio)^-1;
end

function reference_density_value = reference_density(api)
    reference_density_value = 141.5 / (api + 131.5);
end

function true_density_value = true_density(temp, pres, gas_oil_ratio, api, gas_gravity)
    true_density_value = (reference_density(api) + 0.0012 * gas_gravity * gas_oil_ratio) ...
        / b0(temp, gas_oil_ratio, gas_gravity, api);
end

function pb_value = pb(gas_oil_ratio, api, gas_gravity, temp)
    pb_value = ((gas_oil_ratio / (2.03 * gas_gravity))^0.8299) ...
        * exp(-0.02878 * api + 0.00377 * temp) - 0.176;
end
