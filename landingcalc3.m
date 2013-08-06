function [phi_out] = landingcalc3( mu, Rmin, Ratm, P0, H0, Trot, orbit_alt, landing_pe )
% MODIFIED FOR USE BY ENHANCED LANDING CHART GENERATOR
% CHANGED EDGE-CASE RETURN VALUES
% Returns the landing phase angle, where the ship starts in a circular
% orbit and drops its PE down to landing_pe.

% Set to true to log a bunch of details.
details = false;

% Kerbin planetary data from Wiki
dt = 1; % Timestep. Tweak to change accuracy (and speed!)
m = 1;
A = 1;
d = 0.2;
omega_rot = 2*pi/Trot;    % Rotational frequency (rad/s)

% Measure from the planet's centre
r0 = orbit_alt + Rmin;
rpe = landing_pe + Rmin;
Ratm = Ratm + Rmin;

if (rpe > Ratm)
    phi_out = -1;
    fprintf('rpe > Ratm, no atmosphere encounter.\n');
    return;
end

% circ. orbit
v0 = sqrt(mu/r0);

% Step 1: Find the delta-v that when applied retrograde drops you into the
% orbit with the target landing PE.
% I worked this out on paper, too long for this margin to contain, etcetera
% . . .
v1_mag = sqrt(2*mu*(rpe/r0)*(r0-rpe)/(r0^2-rpe^2));
delta_v = abs(v1_mag-v0);

% We now have our orbital parameters for the actual landing orbit.
% (v,r,rpe). I convert these into more useful parameters here.

r = [r0 0];
v = [0 v1_mag];

[ ep ec a hmag rpe rap ] = get_orbit_params( r, v, mu );

% Calculate vcontact and rcontact where the orbit intersects the
% atmosphere.

cos_theta_contact = (1/ec)*(a*(1-ec^2)/Ratm - 1);
theta_contact = acos(cos_theta_contact);

vcontact_mag = sqrt(2*(ep + mu/Ratm));
theta_1 = asin(hmag/(Ratm*vcontact_mag));

vcontact = vcontact_mag * [-cos(theta_1+theta_contact) -sin(theta_1+theta_contact)];
rcontact = Ratm * [cos(theta_contact) sin(theta_contact)];

% Find eccentric anomaly
E = atan2(sqrt(1-ec^2)*sin(theta_contact), ec+cos(theta_contact));

% Find transit time until atmospheric contact.
% Man, "eccentric anomaly" is about archaic as it gets.
% Note that M0 = pi, and the expression is effectively negated to give a
% positive time.
Ttransit = sqrt(a^3/mu)*(pi-E+ec*sin(E));

% Simulate craft in the atmosphere.
F = in_atmo_force( Rmin, P0, H0, d, mu, m, A, omega_rot );
[rfinal vfinal Tatm] = integrate_path( F, m, rcontact, vcontact, dt, Rmin, Ratm );

% Angle through which the planet rotates and angle through which the craft
% travels through the atmosphere.
theta_atm = 180/pi*acos((rcontact*rfinal')/(norm(rcontact)*norm(rfinal)));
theta_rot = 180/pi * omega_rot * (Tatm + Ttransit);

% NOTE: All values measured from planet centre
% Log a bunch of useful orbital data.
if (details)
    fprintf('Landing Maneuver Data\n------------------------------\n');
    fprintf('Landing maneuver delta-v: %f m/s\n', delta_v);
    fprintf('\n');
    fprintf('Initial Orbit Data\n------------------------------\n');
    fprintf('=> Sp. orbital energy: %f\n', ep);
    fprintf('=> Eccentricity: %f\n', ec);
    fprintf('=> Sp. angular momentum: %f\n', hmag);
    fprintf('=> Semi-major axis: %f\n', a);
    fprintf('=> Eccentric anomaly at atmosphere contact: %f\n', E);
    fprintf('=> PE: %f m\n', rpe);
    fprintf('=> AP: %f m\n', (1+ec)*a);
    fprintf('=> Transit time: %f\n', Ttransit);
    fprintf('=> Atmos. contact angle: %f degrees\n', theta_contact*180/pi);
    fprintf('=> Atmospheric entry velocity: %f m/s\n', norm(vcontact));
    fprintf('\n');
    fprintf('Atmospheric Maneuver Data\n------------------------------\n');
    fprintf('=> Time in atmosphere: %f s\n', Tatm);
    fprintf('=> Phase angle of manoeuver: %f degrees\n', ...
        theta_atm);
    if (norm(rfinal) > Rmin)
        fprintf('=> Escaped atmosphere.\n');
        fprintf('=> Exit velocity: %f m/s\n', norm(vfinal));
        fprintf('=> Final Orbit Data\n------------------------------\n');
        fprintf('=> => Sp. orbital energy: %f\n', ep);
        fprintf('=> => Eccentricity: %f\n', ec);
        fprintf('=> => Sp. angular momentum: %f\n', hmag);
        fprintf('=> => Semi-major axis: %f\n', a);
        fprintf('=> => PE: %f\n', rpe);
        fprintf('=> => AP: %f\n', (1+ec)*a);
    else
        fprintf('=> Surface impact!\n');
        fprintf('=> Impact velocity: %f m/s\n', norm(vfinal));
    end
    fprintf('\n');
    fprintf('Planet rotation angle: %f degrees\n', theta_rot);
end

if (norm(rfinal) > Rmin)    % Escaped atmosphere
    phi_out = -1;    % No landing!
    return;
end

% Phase angle to put PE ahead of target for landing.
phi_out = theta_rot+(theta_contact*180/pi-theta_atm);

end

function [ F ] = in_atmo_force( R0, P0, H0, d, mu, m, A, omega_rot )
% Returns a function which gives the total vector force in atmosphere
% as a function of position and velocity.

Kp = 1.2230948554874*0.008;

% Surface velocity for drag calculation.
vs = @(r,v) v - sign(r(1)*v(2)-r(2)*v(1))*[-omega_rot*r(2) omega_rot*r(1)];

% Total force = drag + gravity
F = @(r,v) -0.5*Kp*P0*exp((R0-norm(r))/H0)*norm(vs(r,v))*d*m*A*vs(r,v) - m*(mu/norm(r).^3)*r;

end

function [r v t] = integrate_path( F, m, r0, v0, dt, Rmin, Ratm )
% Integrate the equation of motion to find the path in atmosphere.
% Terminate on collision with surface (sea level) or atmosphere escape.

t = 0;
r = r0;
v = v0;

a = @(rin, vin) F(rin, vin)./m;

firstrun = true;

while firstrun || (norm(r) <= Ratm && norm(r) >= Rmin)
    r_old = r;
    v_old = v;
    a_t = a(r_old, v_old);
    r = r_old + v_old*dt + 0.5*a_t*dt.^2;
    v_est = v_old + 0.5*dt*(a_t + a(r,v_old+dt*a_t));
    v = v_old + 0.5*(a_t + a(r,v_est))*dt;
    t = t + dt;
    firstrun = false;
end
end

function [ ep ec a hmag rpe rap ] = get_orbit_params( r, v, mu )
% Get useful orbit parameters from position vector, velocity vector, and
% gravitational parameter.

% sp. orbital energy
ep = dot(v, v)/2 - mu/norm(r);

%sp angular momentum
r3d = [r(1) r(2) 0];
v3d = [v(1) v(2) 0];
hmag = cross(r3d, v3d);
hmag = hmag(3);

%eccentricity
ec = sqrt(1+2*ep*hmag.^2./(mu.^2));

%semi-major axis
a = -mu/(2*ep);

%periapse distance
rpe = -a*(ec-1);

% Apoapse distance (valid iff ec < 0)
rap = (1+ec)*a;

end