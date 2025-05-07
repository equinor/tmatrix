function fun=betaf(alpha,u,mu,delta)
%Returns the value of the generalizd beta probability 
%function a given alpha (aspect ratio), u (end member of the distribution),
%mu (mean aspect ratio of the distribution), delta (weight parameter of the distribution). 

%the generalizd beta probability functionwhich can be found in Agersborg (2007), phd thesis:
%https://bora.uib.no/handle/1956/2422

%
% Warning: The generalized beta probability function is very unstable. 
%          Consider to use other functions to describe the cracks in the
%          system.

% 09.03.2012
% Remy Agersborg 

p=(u-mu-(mu*delta^2))/(u*delta^2);
q=((u-mu)/mu)*p;
B=(gamma(p)*gamma(q))/gamma(p+q);
en=q-1;
to=p-1;
argen=1-alpha./u;
argto=alpha./u;
fun=(1/(u*B))*(argen).^(en).*(argto).^(to);
