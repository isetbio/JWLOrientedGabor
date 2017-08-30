function y = ogWeibull(vals,x)
%y = Weibull(p,x)
%
%Parameters:  p.b slope
%             p.t threshold yeilding ~80% correct
%             x   intensity values.

%% Examples

exampleOn = 0;

if exampleOn
    vals = [0.5, 0.1];
    x    = [0 : 0.05 : 1];
end

%% pre-defined variables

p.b = vals(1);
p.t = vals(2);
g   = 0.5;        %chance performance
e   = (.5)^(1/3); %threshold performance ( ~80%)

%% compute

k = (-log( (1-e)/(1-g)))^(1/p.b);
y = 1- (1-g)*exp(- (k*x/p.t).^p.b);

%% plot results

if exampleOn
    figure, plot(x, y)
    title('Example function'), xlabel('stimulus levels'), ylabel('% Correct')
end