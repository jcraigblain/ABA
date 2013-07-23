
function kmin = ABA_kScan(k0, kstep, kf, time_points, integrations,kays)
%scans a range of k values for a particular step in the primer extension
%reaction and selections the one that gives the least squared error
%from the measured integrations

k_scan = k0:kstep:kf;
r = zeros(length(k_scan));

initial_integrations = integrations(:,1);
for i=1:length(k_scan)
    test_k=k_scan(i);
    for j=2:length(time_points)
        model = ABA_RateODE([kays test_k],time_points(j),initial_integrations);
        r(i) = r(i) + sum((integrations(:,j)-model(end,:)').^2);
    end
end

[rminn, index] = min(r);
kmin = k_scan(index(1));

end