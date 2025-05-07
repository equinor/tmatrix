function [new_alpha, v_c] = betadistr(delta, alpha_mean, u, crack_density)
%Returns  a crack distribution made from a generalizd beta probability
%function which can be found in Agersborg (2007), phd thesis:
%https://bora.uib.no/handle/1956/2422

%delta	        : Weight parameter of the distribution 
%alpha_mean	    : Mean aspect ratio of the distribution
%u	            : End member of the crack distribution
%crack_density	: Total crack density


% Plots both the crack distribution as function of crack density and
% concentration
%
%The length of new_alpha and v_c is set to 1001
%

%
% Warning: The generalized beta probability function is very unstable. 
%          Consider to use other functions to describe the cracks in the
%          system.

% 09.03.2012
% Remy Agersborg 

total=10001;                       % total number of alphas
alpha=linspace(0,u,total);         % alpha vector

%The beta function normalized 
fun = betaf(alpha,u,alpha_mean, delta);
normalized_fun = fun/sum(fun);

total_crack_density  = sum(normalized_fun)*crack_density;           % check: total crack density

%The crack density function
epsilon_original=normalized_fun.*total_crack_density;

% Resample the crack density function to reduce the numbers of cracks.
res = 10; %Every 'res' the alpha is picked. 
new_alpha = alpha(1:10:end);
% Must normalize the "new" beta function
fun_resampled = normalized_fun(1:10:end);
normalized_fun_resampled = fun_resampled / sum(fun_resampled);
%the resampled crack density function
fun_resampled_result = normalized_fun_resampled.*total_crack_density;


v_c = (4/3)*pi.*new_alpha.*fun_resampled_result;


%plotting the results
figure(50);

xlabel('alpha');
ylabel('\epsilon (Crack density)')
hold on;
plot(new_alpha,fun_resampled_result,'-b');
Titlename = ['Total crack density original: ' , num2str(total_crack_density), ' New: ' num2str(sum(fun_resampled_result))];
title(Titlename)


figure(51)
plot(new_alpha,v_c,'-b');
Titlename = ['New distribution ' ];
title(Titlename)
xlabel('\alpha');
ylabel('v (concentration)')
end




