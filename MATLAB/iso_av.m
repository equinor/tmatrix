function [t_bar] = iso_av(t)
% Returns a (6x6) matrix t_bar averaged over all the orientations (isotropy).

% t : 6x6 matrix which has a HTI symmetry

% 09.03.2012
% Remy Agersborg 

lambda = (t(1,1)+t(3,3)+5*t(1,2)+8*t(1,3)-2*t(4,4))/15;
mu = (7*t(1,1)+2*t(3,3)-5*t(1,2)-4*t(1,3)+6*t(4,4))/30;
c11 = lambda + 2*mu;
c12 = lambda;
c44 = mu;
t_bar = [   
   c11 c12 c12 0     0     0;
   c12 c11 c12 0     0     0;
   c12 c12 c11 0     0     0;
   0   0   0   2*c44 0     0;
   0   0   0   0     2*c44 0;
   0   0   0   0     0     2*c44;
];