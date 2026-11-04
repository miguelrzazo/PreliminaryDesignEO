function days = referenceLifetime(alt)
    % Aproximación basada en la tabla y la tendencia logarítmica
    % 400 km: 1 año (365 días), 500 km: 10 años (3650 días), 
    % 700 km: 100 años (36500 días), 900 km: 1000 años (365000 días)
    % Ajuste logarítmico entre los puntos conocidos
    
    if alt < 400
        days = 0;
    elseif alt <= 500
        % Interpolación logarítmica entre 400 y 500 km
        days = exp(interp1([log(400), log(500)], [log(365), log(3650)], log(alt)));
    elseif alt <= 700
        days = exp(interp1([log(500), log(700)], [log(3650), log(36500)], log(alt)));
    elseif alt <= 900
        days = exp(interp1([log(700), log(900)], [log(36500), log(365000)], log(alt)));
    else
        % Por encima de 900 km, extrapolamos logarítmicamente hasta 2000 km
        days = exp(interp1([log(900), log(2000)], [log(365000), log(3650000)], log(alt)));
    end
end

