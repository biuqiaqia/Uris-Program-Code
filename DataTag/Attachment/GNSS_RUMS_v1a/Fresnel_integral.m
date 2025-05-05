function [F] = Fresnel_integral(x)

% Fresnel_integral
% Numerical method by:
% [1] McNamara, D. A., C. W. I. Pistorius, and J. A. G. Malherbe. 
%     "Introduction to The Uniform Geometrical theory of diffraction." 
%     Artech House, London (1990).

if x >= 0 && x < 0.3
    F = (sqrt(pi*x)-2*x*exp(1i*pi/4)-2/3*x^2*exp(-1i*pi/4))*exp(1i*pi/4)*exp(1i*x);
% elseif x >= 0.3 && x < 0.5
%     F0 = (sqrt(pi*0.3)-2*0.3*exp(1i*pi/4)-2/3*0.3^2*exp(-1i*pi/4))*exp(1i*pi/4)*exp(1i*0.3);
%     F = (0.5729+0.2677i)*(x-0.3);
elseif x >= 0.3 && x < 0.5
    F = (0.6768+0.2682i) + (0.5195+0.0025i)*(x-0.5);
elseif x >= 0.5 && x < 0.7
    F = (0.7439+0.2549i) + (0.3355-0.0665i)*(x-0.7);
elseif x >= 0.7 && x < 1.0
    F = (0.8095+0.2322i) + (0.2187-0.0757i)*(x-1.0);
elseif x >= 1.0 && x < 1.5
    F = (0.8730+0.1982i) + (0.1270-0.068i)*(x-1.5);
elseif x >= 1.5 && x < 2.3
    F = (0.9240+0.1577i) + (0.0638-0.0506i)*(x-2.3);
elseif x >= 2.3 && x < 4.0
    F = (0.9658+0.1073i) + (0.0246-0.0296i)*(x-4.0);
elseif x >= 4.0 && x <= 5.5
    F = (0.9797+0.0828i) + (0.0093-0.0163i)*(x-5.5);
elseif x > 5.5
    F = 1+1i/2/x-3/4/x/x-1i*15/(8*x^3)+75/(16*x^4);    
end


