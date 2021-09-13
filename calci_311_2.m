M = 2.5;   %change 
l = 19.43;  %change  ( %metre)
gamma = 0;  %change (%degrees)
theta = 0;  %change (%degrees)
hv = 15240;    %change (m)
pv = 11600;             %pa 
av = 295;    
hg = 760;    %change (km)
pg = 92300;      %pa  
kr = 2.0;
h = hv-hg;
w = 13166*10;    %change %newton

kl = (sqrt(M^2 - 1)*w*cosd(gamma)*cosd(theta))/(1.4*pv*M^2*l^2);
ks = 0.08;          %change from graph after KL


A = 1/(cosd(gamma) * sqrt(M^2 -1)) ;
B = 1/(cosd(theta) * sqrt(M^2 -1)) ;
C = tand(theta)/sqrt(M^2 -1);   
D = (tand(gamma))^2 +1 ;
phi = atand((tand(theta)*cosd(gamma)*D)/(tand(gamma)+A));

Me = sqrt( 1 + ((A*(1-B*tand(gamma)))^2 / ((A*(tand(gamma)+B))^2 + (C*D)^2)));


d = kd*(h/sqrt(Me^2 - 1));
dx = d*cosd(phi);
dy = d*cosd(phi);
he = sqrt( dy^2 + (h*cosd(gamma) + dx*sind(gamma))^2);

kd = 1.05;      %change from graph 
kp = 1.15;      %change from graph 
kt = 0.84;      %change from graph 

delta_Pmax = kp*kr* sqrt(pv*pg) * ((M^2 - 1)^(1/8)) * (he^(-3/4)) * (l)^(3/4) * ks ;

delta_t = kt *(3.42 / av) * (M/(M^2 - 1)^(3/8))* he^(1/4) * l^(3/4) * ks;


