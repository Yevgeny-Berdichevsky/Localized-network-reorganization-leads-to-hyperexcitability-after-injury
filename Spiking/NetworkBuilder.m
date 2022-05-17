close all;
clc; clear all;

%% Inputs 
% By the default values this code will create a network similar to figure 3Bi, intact sphere and 95% neuronal loss 
% This code will yield to 6 figures representing total and average pre- and post-synaptic
% strengths before cut, immediately after cut, and post sprouting after cut
% To choose what scenario the model to run with, user can change variable
% 'scenario'. 
% .scenario = 1 is the 'most confined' condition: short sprouting distance
% (Fig. 3B(i))
% .scenario = 2 is the 'intermediate' condition: intermediate sprouting
% distance (Fig. 3B(ii))
% .scenario = 3 is the 'same as before' condition: long (developmental)
% sprouting distance (Fig. 3B(iii))
% To test for different number of nodes to be removed, change No_rmnodes (This will create all the plots represented in figure 3A and B)
N_exc = 1800; % Number of excitatory neurons
N_inh = 200; % Number of inhibitory neurons
N = N_exc + N_inh; % Total number of neurons
No_rmnodes = 1900; % How many nodes to remove from the north pole
No_exc_conn = 300; % The setpoint value of connections on excitatory neurons
No_inh_conn = 140; % The setpoint value of connections on inhibitory neurons
radius = 950; % radius of the sphere in um
Loop_No = 1000; % How many loops to run for the development
Loop_No2 = 50000; % How many loops to run for the development after the node removal
N_post = 1:N; % Number of postsynaptic neurons
N_pre = 1:N; % Number of presynaptic neurons
sig_mean = 230; % Mean value for sigmoid function that constrains the nodes ability
%to accumulate connections too rapidly
scenario = 1; % 1 for maximum constraint after cut, 2 for intermediate, and 3 for the same as before
st_node = 1; % starting node to remove the nodes (1 to start from the north pole)

%% Finding the coordinates of each nodes on the surface of the sphere
% [V] = ParticleSampleSphere('N', N); %This function is developed by ...
... Anton Semechko (2022). Suite of functions to perform uniform sampling of a sphere (https://github.com/AntonSemechko/S2-Sampling-Toolbox), GitHub. Retrieved April 19, 2022.
% % 3D Coordinates
% X = V(:,1); 
% Y = V(:,2); 
% Z = V(:,3);
% save('Sphere_coordinates.mat','V');
% Demo coordinates used in the model
V = load('Sphere_coordinates.mat');
X = V.V(:,1)*radius;
Y = V.V(:,2)*radius;
Z = V.V(:,3)*radius;

%% Place inhibitory neurons at equidistant locations to create a uniform distribution
% To find coordinates of inhibitory neurons, a second sphere with the same
% radius but fewer number of nodes is created
[V_inh] = ParticleSampleSphere('N', N_inh);
% Coordinates of inhibitory neurons
X_inh = V_inh(:,1)*radius; 
Y_inh = V_inh(:,2)*radius;
Z_inh = V_inh(:,3)*radius;

% Closest nodes from the original sphere to the second set of nodes are
% found and then designated as inhibitory neurons
for i=1:N
    for j=1:N_inh
        X_diff = X_inh(j) - X(i);
        Y_diff = Y_inh(j) - Y(i);
        Z_diff = Z_inh(j) - Z(i);
        eudist(i,j) = sqrt((X_diff^2)+(Y_diff^2)+(Z_diff^2));
    end
end
min_eudist = min(eudist);

for i=1:N_inh
    inh_indx(i) = find(min_eudist(i)==eudist(:,i)); % Indices of inhibitory neurons
end

%% Set points
max_deg_out = No_exc_conn * ones(N,1); % Matrix that indicates the number of connections each neuron sends
max_deg_in = No_exc_conn * ones(N,1); % Matrix that indicates the number of connections each neuron recieves 
max_deg_out(inh_indx) = No_inh_conn; % Number of connections each inhibitory neuron sends
max_deg_in(inh_indx) = No_inh_conn; % Number of connections each inhibitory neuron recieves 
exc_neuron_indx = find(max_deg_out>No_inh_conn+1); % index of excitatory neurons in NbN matrix
inh_neuron_indx = find(max_deg_out<No_exc_conn); % index of inhibitory neurons in NbN matrix

figure;colormap 'white'; %Plot the inh neurons in red and exc neurons in blue to check for the distributions
[xs,ys,zs] = sphere; surf(radius*xs,radius*ys,radius*zs);hold on;
for i = 1:N
    if max_deg_out(i) < No_exc_conn
        plot3(X(i),Y(i),Z(i),'.', 'MarkerSize', 15,'color', 'r');
    else
        plot3(X(i),Y(i),Z(i),'.', 'MarkerSize', 15,'color', 'b');
    end
end

%% Find arc distances between any two nodes
[Conn_Probability, Dist, Conn_Probability_aftercut] = ...
    ArcDistCalc(N,X,Y,Z,exc_neuron_indx,inh_neuron_indx,radius,scenario);

%% Development
% Target nodes are being picked based on their distance (weighted probability)
tic;
% figure('color','white'); plot(connections,Pick_Prob, 'LineWidth', 3); 
[Conn_dev, syn_strength_dev] = ...
    ConnDev(N,Loop_No,max_deg_out,max_deg_in,Conn_Probability,sig_mean);
toc;

%% Remove nodes and edges
[N_rmv,Affected_prenodes,Affected_postnodes,max_deg_out_rmv,max_deg_in_rmv,...
    Conn_Probability_aftercut,Conn_rmv, syn_strength_rmv,...
    inh_neuron_indx_rmv,exc_neuron_indx_rmv,max_inhin_aftercut,...
    max_excin_aftercut,Dist_rmv,X_rmv,Y_rmv,Z_rmv] = ...
    rmvnodedge(N, No_rmnodes,syn_strength_dev,Conn_dev,X,Y,Z,st_node,...
    Conn_Probability_aftercut,Dist,radius,max_deg_out,max_deg_in);

%% Development after cut
[Conn_dev_aftercut, syn_strength_dev_aftercut] = ...
    ConnDevAftercut(N_rmv,Affected_prenodes,Affected_postnodes,Loop_No2,...
    max_deg_out_rmv,max_deg_in_rmv,Conn_Probability_aftercut, Conn_rmv, syn_strength_rmv,sig_mean,...
    radius,X_rmv,Y_rmv,Z_rmv,max_inhin_aftercut,max_excin_aftercut,inh_neuron_indx_rmv,exc_neuron_indx_rmv);
toc;

%% Save .mat files to hold both the geometry and the connectivity map of intact and lesioned networks 
%These files can be retrieved in the LIF.m to run integrate-and-fire over the network
save('Conn_dev_aftercut_model1.mat','Conn_dev_aftercut'); %0-1 connectivity matrix/map post-lesion sprout, indicating which neurons are connected 
save('syn_strength_dev_aftercut_model1.mat','syn_strength_dev_aftercut'); %matrix containg the synaptic strength post-lesion sprout
save('Conn_dev_model1.mat','Conn_dev'); %0-1 connectivity matrix/map of intact sphere, indicating which neurons are connected 
save('syn_strength_dev_model1.mat','syn_strength_dev'); %matrix containg the synaptic strength of intact network
save('exc_indx_model1.mat','exc_neuron_indx'); %index of excitatory neurons in the intact network
save('inh_indx_model1.mat','inh_neuron_indx'); %index of inhibitory neurons in the intact network
save('Dist_model1.mat','Dist'); %arc distances between enurons

%% Figures
% Plot the developed network 
figure('color','white'); subplot(1,2,1);
scatter3(X,Y,Z,200,sum(abs(syn_strength_dev),2),'filled');
set(gca, 'fontsize', 14); title('Total Post-Synaptic Weight Pre-injury'); 
box off; colormap jet;colorbar;caxis([1 300]);
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
average_post = sum(abs(syn_strength_dev),2)./sum(Conn_dev,2);
subplot(1,2,2); scatter3(X,Y,Z,200,average_post,'filled');
title('Average Post-Synaptic Weight Pre-injury');
box off; set(gca, 'fontsize', 14); set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
colormap jet;colorbar;caxis([1 2]);

figure('color','white');subplot(1,2,1);
scatter3(X,Y,Z,100,sum(abs(syn_strength_dev)),'filled');
set(gca, 'fontsize', 14); title('Total Pre-Synaptic Weight Pre-injury');
box off; colormap jet;colorbar;caxis([1 300]);
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
average_pre = sum(abs(syn_strength_dev))./sum(Conn_dev);
subplot(1,2,2); scatter3(X,Y,Z,100,average_pre,'filled');
set(gca, 'fontsize', 14); title('Average Pre-Synaptic Weight Pre-injury');
box off; set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
colormap jet;colorbar;caxis([1 2]);

% Network immediately after node removal
figure('color','white'); subplot(1,2,1);
scatter3(X_rmv,Y_rmv,Z_rmv,200,sum(abs(syn_strength_rmv),2),'filled');
set(gca, 'fontsize', 14); title('Total Post-Synaptic Weight Pre-sprouting');
box off; colormap jet;colorbar;caxis([1 300]);
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);caxis([1 300]);
average_post = sum(abs(syn_strength_rmv),2)./sum(Conn_rmv,2);
subplot(1,2,2); scatter3(X_rmv,Y_rmv,Z_rmv,200,average_post,'filled');
title('Average Post-Synaptic Weight Pre-sprouting');
box off; set(gca, 'fontsize', 14); 
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
colormap jet;colorbar;caxis([1 2]);

figure('color','white'); subplot(1,2,1);
scatter3(X_rmv,Y_rmv,Z_rmv,200,sum(abs(syn_strength_rmv)),'filled');
set(gca, 'fontsize', 14); title('Total Pre-Synaptic Weight Pre-sprouting');
box off; colormap jet;colorbar;caxis([1 300]);
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);caxis([1 300]);
average_pre = sum(abs(syn_strength_rmv))./sum(Conn_rmv);
subplot(1,2,2); scatter3(X_rmv,Y_rmv,Z_rmv,100,average_pre,'filled');
set(gca, 'fontsize', 14); title('Average Pre-Synaptic Weight Pre-sprouting');
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
box off; 
colormap jet;colorbar;caxis([1 2]);

% Network after sprouting
figure('color','white');subplot(1,2,1);
scatter3(X_rmv,Y_rmv,Z_rmv,200,sum(abs(syn_strength_dev_aftercut),2),'filled');
set(gca, 'fontsize', 14); title('Total Post-Synaptic Weight Post-sprouting');
box off; colormap jet;colorbar;
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);caxis([1 300]);
average_post_aftercut = sum(abs(syn_strength_dev_aftercut),2)./sum(Conn_dev_aftercut,2);
subplot(1,2,2);scatter3(X_rmv,Y_rmv,Z_rmv,200,average_post_aftercut,'filled');
title('Average Post-Synaptic Weight Post-sprouting');
box off; set(gca, 'fontsize', 14); 
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
colormap jet;colorbar;caxis([1 2]);

figure('color','white');subplot(1,2,1);
scatter3(X_rmv,Y_rmv,Z_rmv,100,sum(abs(syn_strength_dev_aftercut)),'filled');
set(gca, 'fontsize', 14); title('Total Pre-Synaptic Weight Post-sprouting');
box off; colormap jet;colorbar;caxis([min(sum(syn_strength_dev_aftercut)) max(sum(syn_strength_dev_aftercut))]);
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);caxis([1 300]);
average_pre_aftercut = sum(abs(syn_strength_dev_aftercut))./sum(Conn_dev_aftercut);
subplot(1,2,2);
scatter3(X_rmv,Y_rmv,Z_rmv,100,average_pre_aftercut,'filled');
set(gca, 'fontsize', 14); title('Average Pre-Synaptic Weight Post-sprouting');
set(gca,'TickDir','out','fontsize', 14, 'linewidth', 1.2);
box off; 
colormap jet;colorbar;caxis([1 2]);

%% Probability of connection vs distance for excitatory neurons
% Conn_dev_exc = Conn_dev;
% Conn_dev_exc(inh_neuron_indx,:)=0;
% Conn_distance = Conn_dev_exc .* Dist; Conn_distance(Conn_distance == 0) = [];
% all_exc_dist = unique(Conn_distance); all_exc_dist(all_exc_dist==0)=[];
% bin_edges = all_exc_dist; 
% Dist2 = Dist(exc_neuron_indx,:);
% Dist2(Dist2 == 0) = []; hist_Dist = histcounts(Dist2,bin_edges); hist_Conn = histcounts(Conn_distance,bin_edges);
% Conn_Percent = (hist_Conn ./ hist_Dist);
% Conn_Percent2 = Conn_Percent(~isnan(Conn_Percent));
% NaNs = isnan(Conn_Percent) == 1;  bin_edges(1,NaNs) = 0; bin_edges(bin_edges == 0) = [];
% bins = bin_edges(1:end-1);
% % Conn_Fit = fit(bins', Conn_Percent2', 'gauss1'); % Find the best fit of Gaussian distribution
% figure('color','white');
% plot(bins, Conn_Percent2, 'k.', 'MarkerSize', 12); hold on; 
% % plot(Conn_Fit, 'b'); hold on;
% set(gca, 'fontsize', 16, 'LineWidth', 1); 
% set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);


%% Probability of connection vs distance for inhibitory neurons
% Conn_dev_inh = Conn_dev;
% Conn_dev_inh(exc_neuron_indx,:)=0;
% Conn_distance = Conn_dev_inh .* Dist; Conn_distance(Conn_distance == 0) = [];
% all_inh_dist = unique(Conn_distance); all_inh_dist(all_inh_dist==0)=[];
% bin_edges = all_inh_dist; 
% Dist2 = Dist(inh_neuron_indx,:);
% Dist2(Dist2 == 0) = []; hist_Dist = histcounts(Dist2,bin_edges); hist_Conn = histcounts(Conn_distance,bin_edges);
% Conn_Percent = (hist_Conn ./ hist_Dist);
% Conn_Percent2 = Conn_Percent(~isnan(Conn_Percent));
% NaNs = find(isnan(Conn_Percent) == 1);  bin_edges(1,NaNs) = 0; bin_edges(bin_edges == 0) = [];
% bins = bin_edges(1:end-1);
% % Conn_Fit = fit(bins', Conn_Percent2', 'gauss1'); % Find the best fit of Gaussian distribution
% % figure('color','white');
% plot(bins, Conn_Percent2, 'r.', 'MarkerSize', 12); hold on; 
% % plot(Conn_Fit, 'b'); hold on;
% set(gca, 'fontsize', 16, 'LineWidth', 1); 
% set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% box off; xlabel('Distance(um)'); ylabel('Probability of connection')
% legend('Excitatory','Inhibitory'); legend boxoff;

%% one neuron tracker

% neuron1_post = Conn_dev_aftercut(1,:) .* Dist_rmv(1,:);
% neuron1_post_indx = find(Conn_dev_aftercut(1,:)==1);
% figure;colormap 'white';
% [xs,ys,zs] = sphere; surf(radius*xs,radius*ys,radius*zs);
% hold on;plot3(X(1),Y(1),Z(1),'.', 'MarkerSize', 60,'color', 'k'); 
% hold on;plot3(X(neuron1_post_indx),Y(neuron1_post_indx),Z(neuron1_post_indx),'.', 'MarkerSize', 30,'color', 'r'); 
% 
% 
% neuron1_pre = Conn_dev_aftercut(:,1) .* Dist_rmv(:,1);
% neuron1_pre_indx = find(Conn_dev_aftercut(:,1)==1);
% figure;colormap 'white';
% [xs,ys,zs] = sphere; surf(radius*xs,radius*ys,radius*zs);
% hold on;plot3(X(1),Y(1),Z(1),'.', 'MarkerSize', 60,'color', 'k'); 
% hold on;plot3(X(neuron1_pre_indx),Y(neuron1_pre_indx),Z(neuron1_pre_indx),'.', 'MarkerSize', 30,'color', 'b');

%% Histogram of connections 
% N = 2000;
% inh_exc_mat = ones(N);
% inh_mat = zeros(N);
% exc_mat = zeros(N);
% inh_exc_mat(inh_neuron_indx,:) = -1;
% inh_mat(inh_neuron_indx,:) = -1;
% exc_mat(exc_neuron_indx,:) = 1;
% syn_strength_inh = syn_strength_dev .* inh_mat;
% syn_strength_exc = syn_strength_dev .* exc_mat;
% total_post_inh_syn = sum(syn_strength_inh,2);
% total_post_inh_syn(total_post_inh_syn==0) = [];
% total_pre_inh_syn = sum(syn_strength_inh);
% total_post_exc_syn = sum(syn_strength_exc,2);
% total_post_exc_syn(total_post_exc_syn==0) = [];
% total_pre_exc_syn = sum(syn_strength_exc);
% 
% 
% syn_strength_dev = syn_strength_dev .* inh_exc_mat;
% inh_in = sum(syn_strength_dev.*(syn_strength_dev<0));
% exc_in = sum(syn_strength_dev.*(syn_strength_dev>0));
% inh_inh_in = inh_in(inh_neuron_indx);
% exc_inh_in = exc_in(inh_neuron_indx);
% ratio_inh = [-1*inh_inh_in;exc_inh_in]';
% ratio_inh_vect = abs(inh_inh_in) ./ exc_inh_in;
% mean_ratio_inh = mean(ratio_inh_vect);
% stdev_ratio_inh = std(ratio_inh_vect);
% 
% inh_exc_in = inh_in(exc_neuron_indx);
% exc_exc_in = exc_in(exc_neuron_indx);
% ratio_exc = [-1*inh_exc_in;exc_exc_in]';
% ratio_exc_vect = abs(inh_exc_in) ./ exc_exc_in;
% mean_ratio_exc = mean(ratio_exc_vect);
% stdev_ratio_exc = std(ratio_exc_vect);
% 
% figure('color','white');
% bar(ratio_inh,'stacked');
% set(gca, 'fontsize', 16); 
% ylabel('Synaptic weight onto inh');
% xlabel('Neuron No');
% box off; set(gca,'TickDir','out','fontsize', 24, 'linewidth', 1.2);
% legend('inh','exc'); legend boxoff;
% 
% figure('color','white');
% bar(ratio_exc,'stacked');
% set(gca, 'fontsize', 16); 
% ylabel('Synaptic weight onto exc');
% xlabel('Neuron No');
% box off; set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% legend('inh','exc'); legend boxoff;
% 
% % [h,p] = kstest2(ratio_inh_vect,ratio_exc_vect);
% 
% inh_exc_mat_rmv = ones(N_new);
% inh_mat_rmv = zeros(N_new);
% exc_mat_rmv = zeros(N_new);
% inh_exc_mat_rmv(inh_neuron_indx_rmv,:) = -1;
% inh_mat_rmv(inh_neuron_indx_rmv,:) = -1;
% exc_mat_rmv(exc_neuron_indx_rmv,:) = 1;
% syn_strength_inh_rmv = syn_strength_dev_aftercut .* inh_mat_rmv;
% syn_strength_exc_rmv = syn_strength_dev_aftercut .* exc_mat_rmv;
% 
% syn_strength_dev_aftercut = syn_strength_dev_aftercut .* inh_exc_mat_rmv;
% inh_in = sum(syn_strength_dev_aftercut.*(syn_strength_dev_aftercut<0));
% exc_in = sum(syn_strength_dev_aftercut.*(syn_strength_dev_aftercut>0));
% inh_inh_in_rmv = inh_in(inh_neuron_indx_rmv);
% exc_inh_in_rmv = exc_in(inh_neuron_indx_rmv);
% ratio_inh_rmv = [-1*inh_inh_in_rmv;exc_inh_in_rmv]';
% ratio_inh_vect_rmv = abs(inh_inh_in_rmv) ./ exc_inh_in_rmv;
% inh_exc_in_rmv = inh_in(exc_neuron_indx_rmv);
% exc_exc_in_rmv = exc_in(exc_neuron_indx_rmv);
% ratio_exc_rmv = [-1*inh_exc_in_rmv;exc_exc_in_rmv]';
% ratio_exc_vect_rmv = abs(inh_exc_in_rmv) ./ exc_exc_in_rmv;
% 
% figure('color','white');
% bar(ratio_inh_rmv,'stacked');
% set(gca, 'fontsize', 24); 
% ylabel('Synaptic weight onto inh');
% xlabel('Neuron No');
% box off; set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% legend('inh','exc'); legend boxoff;
% figure('color','white');
% bar(ratio_exc_rmv,'stacked');
% set(gca, 'fontsize', 24); 
% ylabel('Synaptic weight onto exc');
% xlabel('Neuron No');
% box off; set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% legend('inh','exc'); legend boxoff;
% % 
% % ratio_beforepost = [ratio_exc_vect;ratio_exc_vect_rmv];
% N_exc = 1800;
% N_inh = 200;
% ratio_exc_vect_rmv_all = zeros(N_exc,1);
% ratio_exc_vect_rmv_all((length(ratio_exc_vect)-length(ratio_exc_vect_rmv)+1):end)=ratio_exc_vect_rmv;
% ratio_exc_intact_lesioned = [ratio_exc_vect;ratio_exc_vect_rmv_all'];
% figure('color','white');
% bar(ratio_exc_intact_lesioned',1);
% set(gca, 'fontsize', 16); ylabel('Ratio inh/exc'); xlabel('neurons')
% box off; set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% legend('intact','lesioned'); legend boxoff;
% 
% % 
% % [h1,p1] = kstest2(ratio_exc_vect(100:400),ratio_exc_vect_rmv(1:300));
% % [h2,p2] = kstest2(ratio_inh_vect,ratio_inh_vect_rmv);
% 
% % s1 = mwwtest(ratio_exc_vect,ratio_exc_vect_rmv);
% % s2 = mwwtest(ratio_inh_vect,ratio_inh_vect_rmv);
% mean_rmv = mean(ratio_exc_vect_rmv);
% stdev_rmv = std(ratio_inh_vect_rmv);
% 
% variance_exc = var(ratio_exc_vect);
% variance_inh = var(ratio_inh_vect);
% 
% variance_exc_rmv = var(ratio_exc_vect_rmv);
% variance_inh_rmv = var(ratio_inh_vect_rmv);
% 
% set(gcf,'renderer','Painters');
% 
% exc_syn_strength = sum(syn_strength_dev(exc_neuron_indx,:));
% inh_syn_strength = sum(syn_strength_dev(inh_neuron_indx,:));
% ratio = abs(inh_syn_strength)./exc_syn_strength;
% ratio_mean = mean(ratio);
% ratio_std = std(ratio);
% exc_syn_strength_aftercut = sum(syn_strength_dev_aftercut(exc_neuron_indx_rmv,:));
% inh_syn_strength_aftercut = sum(syn_strength_dev_aftercut(inh_neuron_indx_rmv,:));
% ratio_aftercut = abs(inh_syn_strength_aftercut)./exc_syn_strength_aftercut;
% ratio_mean_aftercut = mean(ratio_aftercut);
% ratio_aftercut_std = std(ratio_aftercut);
% 
% result = [ratio_mean ratio_mean_aftercut; ratio_std ratio_aftercut_std];