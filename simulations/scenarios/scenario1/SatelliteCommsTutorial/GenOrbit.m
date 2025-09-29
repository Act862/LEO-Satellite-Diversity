function Orbit = GenOrbit(Param, Elem)
    % Extract parameters
    mu = Param.mu;
    t = Param.t;

    % Extract elements
    a = Elem.a;
    e = Elem.e;
    i = deg2rad(Elem.Inc);
    omega = deg2rad(Elem.omega);
    RAAN = deg2rad(Elem.RAAN);
    TA0 = deg2rad(Elem.TA);

    % Precompute constants
    n = sqrt(mu / a^3); % Mean motion

    % Prepare output
    Orbit.r = zeros(length(t), 3);

    for idx = 1:length(t)
        % Mean anomaly
        M = n * t(idx);

        % Solve Kepler's Equation for Eccentric Anomaly
        E = KeplerSolver(M, e);

        % True Anomaly
        TA = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

        % Radius
        r_mag = a * (1 - e * cos(E));

        % Position in Perifocal frame
        r_perifocal = [r_mag * cos(TA); r_mag * sin(TA); 0];

        % Rotation matrices
        R3_W = [cos(RAAN), -sin(RAAN), 0; sin(RAAN), cos(RAAN), 0; 0, 0, 1];
        R1_i = [1, 0, 0; 0, cos(i), -sin(i); 0, sin(i), cos(i)];
        R3_w = [cos(omega), -sin(omega), 0; sin(omega), cos(omega), 0; 0, 0, 1];

        Q_pX = R3_W * R1_i * R3_w;

        % Position in ECI frame
        r_ECI = Q_pX * r_perifocal;

        % Save to orbit structure
        Orbit.r(idx, :) = r_ECI';
    end
end
