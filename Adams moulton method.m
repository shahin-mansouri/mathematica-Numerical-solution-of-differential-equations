% Shahin Mansouri
% adams moulton method
% ===========================================================

clear; clc; format compact;

% -------------------------Code ------------------------------

f=@(t,y) 3*t+y/t;
alpha=5;
a=1;
b=2;
n=3;
[t, w, h] = abs2(f, a, b, alpha, n)

% ------------------------function----------------------------------

function [t, w, h] = abs2(f, a, b, alpha, n)
    %AB2 Two-step Adams Bashforth method
    % [t, w, h] = ab2(f, a, b, alpha, n) performs the two-step Adams Bashforth 
    % method for solving the IVP y' = f(t,y) with initial condition y(a) = alpha
    % taking n steps from t = a to t = b. The first step from t = a to t = a + h
    % is performed using the modified Euler method.
    h=(b-a)/n;
    t=a:h:b;
    w=zeros(1,length(t));
    w(1)=alpha;
    %modified euler method
    for i=2
        k1=h*f(t(i),w(i));
        k2=h*f(t(i)+h,w(i)+k1);
        w(i+1)=w(i)+1/2*(k1+k2);
    end
    for i=3:length(t)
        w(i+1)=w(i)+(3/2)*h*f(t(i),w(i))-.5*h*f(t(i-1),w(i-1));
    end
end