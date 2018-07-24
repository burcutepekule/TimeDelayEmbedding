function [X,T] = mckyGlss(a,b,tau,x0,deltat,sample_n,var,solver)
% RELATED ARTICLE : http://www.scholarpedia.org/article/Mackey-Glass_equation

time  = 0;
index = 1;
history_length = floor(tau/deltat);
x_history = zeros(history_length, 1); % here we assume x(t)=0 for -tau <= t < 0
x_t = x0;

X = zeros(sample_n+1, 1); % vector of all generated x samples
T = zeros(sample_n+1, 1); % vector of time samples

for i = 1:sample_n+1
    X(i) = x_t;
    if tau == 0
        x_t_minus_tau = 0.0;
    else
        x_t_minus_tau = x_history(index);
    end
    
    if (strcmp(solver,'rk4')==1)
        x_t_plus_deltat = mackeyglass_rk4(x_t, x_t_minus_tau, deltat, a, b);
    else
        x_t_plus_deltat = mackeyglass_euler(x_t, x_t_minus_tau, deltat, a, b);
    end
    if (tau ~= 0)
        x_history(index) = x_t_plus_deltat;
        index = mod(index, history_length)+1;
    end
    time = time + deltat;
    T(i) = time;
    x_t = x_t_plus_deltat + sqrt(var)*randn;
end
X = X';
end

