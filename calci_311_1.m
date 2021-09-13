M = 2;                        %change Mach Number
l = 22;                       %change Length of Fuselage (metre)
gamma = 0;                    %change Flight Path Angle (degrees)
theta = 0;                    %change Ray-Path Angle (degrees)
hv = 15000;                   %change height above ground (m)
w = 45000;                    %change weight of Aircraft (kg)
pv = 12110;                   %atmospheric pressure at aircraft (vehicle) altitude (Pa)
av = 295;                     %speed of sound at aircraft (vehicle) altitude (m/sec)
hg = 760;                     %change altitude of ground above sea level (km)
pg = 92300;                   %atmospheric pressure at ground level (Pa) 
kr = 2.0;                     %reflection factor, assumed to be 2.0
h = hv-hg;                    %altitude of aircraft above ground, hv - hg (km)

temp = -56.5 + 273;
av = sqrt(1.4*287*temp);

kl = (sqrt(M^2 - 1)*w*cosd(gamma)*cosd(theta))/(1.4*pv*M^2*l^2);
ks = 0.08;                    %change from graph with Kl
                              %Or use this for approximate value ks = 0.74*(sqrt(kl + 0.027))

A = 1/(cosd(gamma) * sqrt(M^2 -1)) ;
B = 1/(cosd(theta) * sqrt(M^2 -1)) ;
C = tand(theta)/sqrt(M^2 -1);   
D = (tand(gamma))^2 +1 ;
phi = atand((tand(theta)*cosd(gamma)*D)/(tand(gamma)+A));

Me = sqrt( 1 + ((A*(1-B*tand(gamma)))^2 / ((A*(tand(gamma)+B))^2 + (C*D)^2)));

kd = 1.1;                     %change from graph for M
d = kd*(h/sqrt(Me^2 - 1));
dx = d*cosd(phi);
dy = d*cosd(phi);
he = sqrt( dy^2 + (h*cosd(gamma) + dx*sind(gamma))^2);


kp = 1.1;                     %change from graph for M
kt = 0.85;                    %change from graph for M

delta_Pmax = kp*kr* sqrt(pv*pg) * ((M^2 - 1)^(1/8)) * (he^(-3/4)) * (l)^(3/4) * ks ;

delta_t = kt *(3.42 / av) * (M/(M^2 - 1)^(3/8))* he^(1/4) * l^(3/4) * ks;


