function [ b1 , b2, tau_out, gamma_out ] = check_for_nan_enteries_connected( in_vector1, in_vector2, tau, gamma )

%A check if there is a NaN entry in vector 1 or corresponding vector 2. If
%there is, the entery is removed in both vectors in addition to gamma and
%tau. If in_vector1 is the concentrations, v, the in_vector2 must be the '
%aspect ratios, alpha, or opposite.  

% 09.03.2012
% Remy Agersborg 

%checking vector 1
a1 = in_vector1(~isnan(in_vector1)); %vector 1 is now only numbers
a2 = in_vector2(~isnan(in_vector1)); %there may still be some NaN enteries in a2 (vector 2)
%sorting the corresponding tau and gamma vector
tau_1   = tau(~isnan(in_vector1)); %there may still be some NaN enteries in a2 (vector 2)
gamma_1 =  gamma(~isnan(in_vector1)); %there may still be some NaN enteries in a2 (vector 2)

%checking vector 2
b1 = a1(~isnan(a2));
b2 = a2(~isnan(a2));  
tau_out   =  tau_1(~isnan(a2)); 
gamma_out =  gamma_1(~isnan(a2)); 

%both vectors should be only numbers and b1(i) corresponds to b2(i). 

end