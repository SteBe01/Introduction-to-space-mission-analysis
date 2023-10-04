function plotApsidalAxis(a, e, i, OM, om, mu)

% Function to plot the apsidal axis of an orbit
%
% plotApsidalAxis(a, e, i, OM, om, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [Km]
% e             [1x1]   eccentricity                    [-]
% i             [1x1]   inclination                     [rad]
% OM            [1x1]   RAAN                            [rad]
% om            [1x1]   pericenter anomaly              [rad]
% mu            [1x1]   gravitational parameters        [Km^3/s^2]

    r_peri = parorb2rv(a, e, i, OM, om, 0, mu);
    r_apo = parorb2rv(a, e, i, OM, om, pi, mu);

    line([r_peri(1) r_apo(1)], [r_peri(2) r_apo(2)], [r_peri(3) r_apo(3)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 1.25);
end