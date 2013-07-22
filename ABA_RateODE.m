function x = ABA_RateODE(k,tf,initial_integrations,kays)
%Kinetic model of primer extension reaction
%Models each step as a pseudo-first order reaction
%Takes array kays that contains fitted rate constants for primer -> n
%Takes k, rate constant for n -> n+1
%Works for any n

%Invoke ODE Solver

[t, x] = ode15s(@derivatives, [0 tf], initial_integrations);

%Defining ODE system for ode solver
function xdot = derivatives(t, x)
    if kays(1) == 0 %fitting primer -> n+1
        xdot = zeros(2,1);
        xdot(1) = -k*x(1); %primer consumption
        xdot(2) = k*x(1); %production of n+1 and above
    elseif length(kays) == 1 %fitting n+1 -> n+2
        xdot = zeros(length(kays)+2,1);
        xdot(1) = -kays(1)*x(1); %primer consumption
        xdot(2) = -k*x(2)+kays(1)*x(1); %n+1 -> n+2
        xdot(3) = k*x(2); %n+2 -> n+3 and above
    else %fitting n+m-1 -> n+m
        xdot = zeros(length(kays)+2,1);
        xdot(1) = -kays(1)*x(1); %primer consumption
        for i=2:length(kays)
            xdot(i) = -kays(i)*x(i) + kays(i-1)*x(i-1); %n+i-1 -> n+i
        end
        xdot(length(kays)+1) = -k*x(length(kays)+1)+kays(length(kays))*x(length(kays)); %n+m-2 -> n+m-1
        xdot(length(kays)+2) = k*x(length(kays)+1); %n+m-1 -> n+m
    end
end

end
