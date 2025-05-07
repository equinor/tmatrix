function [ alpha,v ] = check_alpha_concentration( alpha, v )
% Checks if v(i) < or = alpha(i). 
% if v(i) > alpha(i) display a warning and set v(i) = alpha(i).

% 09.03.2012
% Remy Agersborg 

for j = 1: length(alpha)
    if alpha(j) < v(j)
        name1 = ['v(' , num2str(j) , ') is  larger than alpha(', num2str(j) ,')' ];
        warning(name1);
        name2 = ['v(' , num2str(j) , ') is  set to alpha(', num2str(j) , ')/2'];
        warning(name2);
        v(j) = alpha(j)/2;
    end
end
end