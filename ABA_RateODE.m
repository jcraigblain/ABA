function x = ABA_RateODE(kays,tf,initial_integrations)
%Kinetic model of primer extension reaction
%Models each step as a pseudo-first order reaction
%Works for any n

%Invoke ODE Solver

[t, x] = ode15s(@derivatives, [0 tf], initial_integrations);

%Defining ODE system for ode solver
function xdot = derivatives(t,x)
    num_bands = length(kays)+1;
    xdot = zeros(num_bands,1);
    xdot(1) = -kays(1)*x(1); %consumption of primer
    for i=2:num_bands-1
       xdot(i) = kays(i-1)*x(i-1)-kays(i)*x(i); %production and consumption of intermediate bands
    end
    xdot(num_bands) = kays(end)*x(num_bands-1); %production of last product
end

end
