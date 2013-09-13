function kays = ABA(title, k0, kstep, kf, time_points, integrations)
%Fit k values to multiple steps in a primer extension reaction
%Takes a matrix of gel integrations (increasing time left to right,
%increasing extension length bottom to top)
%Fits a k value to each extension step one at a time

%normalize gel lanes
column_sums=sum(integrations,1);
for i=1:size(integrations,2)
   integrations(:,i)=integrations(:,i)/column_sums(i); 
end

times = time_points - time_points(1); %t0 now 0

flipped_integrations = flipud(integrations); %primer now at top

%for each product length n, add all integrations for >=n+1 into n+1 band
%then fit a kay value for n->n+1
kays=[];
[num_bands,num_timepoints] = size(flipped_integrations);
for b=2:num_bands
    reduced_integrations = [flipped_integrations(1:b-1,:); sum(flipped_integrations(b:end,:),1)];
    kays(b-1) = ABA_kScan(k0, kstep, kf, times, reduced_integrations, kays);
end

%Model fitted k values
initial_integrations = reduced_integrations(:,1);
model_output = initial_integrations;
for t=2:length(times)
    model_result = ABA_RateODE(kays,times(t),initial_integrations);
    model_output(:,t) = model_result(end,:)';
end

%Plot model output
model_output = flipud(model_output);
model_heatmap_data = flipud(model_output);
y_labels = 0:length(kays);
hm1 = HeatMap(model_heatmap_data,'RowLabels',y_labels,'ColumnLabels',time_points,'Symmetric',false,'DisplayRange',1,'Colormap',colormap(flipud(gray)),'Standardize','none');
addXLabel(hm1,'Time (h)');
addYLabel(hm1,'Extension');
addTitle(hm1,strcat(title,' Model'));
plot(hm1,figure(1))

%Plot input data
input_heatmap_data = flipud(integrations);
y_labels = 0:length(kays);
hm2 = HeatMap(input_heatmap_data,'RowLabels',y_labels,'ColumnLabels',time_points,'Symmetric',false,'DisplayRange',1,'Colormap',colormap(flipud(gray)),'Standardize','none');
addXLabel(hm2,'Time (h)');
addYLabel(hm2,'Extension');
addTitle(hm2,strcat(title,' Data'));
plot(hm2,figure(2))
end

