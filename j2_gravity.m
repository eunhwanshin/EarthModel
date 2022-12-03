function [g,N] = j2_gravity(r)
% function [g,N] = j2_gravity(r)
% Computes the normal gravity vector in the ECEF
% Ref: Hsu, D.Y. (1996). Comparison of Four Gravity Models. IEEE PLANS.
% @param[in]   r   ECEF position vector (m)
% @param[out]  g   Normal gravity vector in the ECEF (m/s^2)
% @param[out]  N   partial g/partial r

WGS84_A2 = 4.0680631590769e+13; % Semi-major axis squared (m^2)
WGS84_J2 = 1.0826298213133048e-3; % -sqrt(5)*C20
WGS84_GM = 3986004.418e+8; % (gravitational constant)*(mass of earth) (m^3/s^2)
WGS84_WE2 = 5.3174941173224997e-9; % earth rotation rate squared (rad/s)^2

z2 = r(3)^2;
r2 = r(1)^2 + r(2)^2 + z2;
inv_r2 = 1.0/r2;
a2_r2 = WGS84_A2 * inv_r2;
z2_r2 = z2 * inv_r2;
r3 = r2^1.5;
t1 = 1.5 * WGS84_J2 * a2_r2;
t2 = 1.0 - 5.0* t1 * z2_r2;
t3 = WGS84_GM / r3;
t4 = WGS84_WE2 - t3 * (t1+t2);

g = [r(1)*t4; r(2)*t4 ; -r(3)*t3*(3*t1+t2)];

N = zeros(3,3);

N(1,1) = t3*(3.0*r(1)^2 * inv_r2 - 1) +  WGS84_WE2;
N(1,2) = t3 * 3.0 * r(1) * r(2) * inv_r2;
N(1,3) = t3 * 3.0 * r(1) * r(3) * inv_r2;

N(2,1) = N(1,2);
N(2,2) = t3 * (3.0*r(2)^2 * inv_r2 - 1.0) +  WGS84_WE2;
N(2,3) = t3 * 3.0 * r(2) * r(3) * inv_r2;

N(3,1) = N(1,3);
N(3,2) = N(2,3);
N(3,3) = t3 * (3.0 * z2_r2 - 1.0);

end