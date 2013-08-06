function [] = enhancedlandingplots()
clear all;

% Generate charts and data for each planet.
% Uses landingcalc3

planets = [
    mkplanet('Eve', 8.1717302e12, 7e5, 96708.574, 5, 7000, 80500)
    mkplanet('Kerbin', 3.5316e12, 6e5, 69077.553, 1, 5000, 21600)
    mkplanet('Duna', 3.0136321e11, 3.2e5, 41446.532, 0.2, 3000, 65517.859)
    mkplanet('Laythe', 1.962e12, 5e5, 55262.042, 0.8, 4000, 52980.879)
    ];
for i = 1:numel(planets)
    genplot(planets(i));
end
end

function [ planet ] = mkplanet( name, mu, Rmin, Ratm, P0, H0, Trot )
    planet = struct('name',name,'mu',mu,'Rmin',Rmin,'Ratm',Ratm,'P0',P0,'H0',H0,'Trot',Trot);
end

function [] = genplot( planet )
    figure();
    boundfunc = @(orbit_alt, landing_alt) landingcalc3( planet.mu, planet.Rmin, planet.Ratm, planet.P0, planet.H0, planet.Trot, orbit_alt, landing_alt);
    alts = linspace(1.01*planet.Ratm, 1e6, 50);
    landing_alts = arrayfun(@(o_alt) fzero(@(l_alt) boundfunc(o_alt, l_alt), planet.Ratm/2), alts);
    errorPlus = arrayfun(@(o_alt, l_alt) boundfunc(o_alt, l_alt+1), alts, landing_alts);
    errorMinus = arrayfun(@(o_alt, l_alt) boundfunc(o_alt, l_alt-1), alts, landing_alts);
    errors = (planet.Rmin*pi/180).*(abs(errorPlus) + abs(errorMinus))./2;
    [AX H1 H2] = plotyy(alts, landing_alts, alts, errors);
    xlabel('Initial Orbit Altitude (m)');
    set(get(AX(1), 'YLabel'), 'String', 'Landing Periapsis (m) (Place Above Target)');
    set(get(AX(2), 'YLabel'), 'String', 'Landing Site Error per m (m)');
    titlestr = sprintf('Enhanced Landing Chart (%s)', planet.name);
    title(titlestr);
    grid minor;
    tableTitle = sprintf('Landing PE Values for %s\n', planet.name);
    fprintf(tableTitle);
    fprintf('===================================================\n');
    fprintf('Orbit Altitude (m)\tLanding PE (m)\tError per m (m)\n');
    fprintf('===================================================\n');
    for i=1:numel(alts)
        fprintf('%.0f\t\t\t\t%.0f\t\t\t%.0f\n', alts(i), landing_alts(i), errors(i));
    end
end