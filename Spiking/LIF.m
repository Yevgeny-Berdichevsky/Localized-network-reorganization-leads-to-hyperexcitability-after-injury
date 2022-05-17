close all;
clc; clear all;

%% DEFINE PARAMETERS
%The default parameters will yield in a raster plot similar to what is shown in figure 4F
%Since the code is generated based on random selections, each run is slightly different
%This code will yield two figures: one indicates the membrane potential of a neuron
%indexed 1 in the map and the other represent the raster plot of neurons sorted from first neuron closest to the edge. 
%To assess the network susceptibility to generate ictal-like events under different number of spontaneous action potentials,
%the variable 'spont_frequency_exc' can be adjusted. 
dt = 1; %time step [ms]
t_end = 10000; %total time of run [ms]
E_L = -65; %resting membrane potential [mV]

V_th = -50; %spike threshold [mV]
V_reset = -75; %value to reset voltage to after a spike [mV]
V_spike = 20; %value to draw a spike to, when cell spikes [mV]
E_syn_exc = 0;  %EPSP reversal potential [mV]
E_syn_inh = -70; %IPSP reversal potential [mV]

Rm_exc = 150; %membrane resistance of excitatory neuron[MOhm]
Rm_inh = 400; %membrane resistance of inhibitory neuron[MOhm]
tau_exc = 50; %membrane time constant [ms]
tau_inh = 5; %membrane time constant [ms]

g_syn_max_exc = 0.0015; %max conductivity of excitatory synapse [uS]
g_syn_max_inh = 0.04; % max conductivity of inhibitory synapse [uS]
tau_syn_exc = 2; %time constant of synapse [ms]
tau_syn_inh = 2; %time constant of synapse [ms]

AP_delay = 5; %time delay for action potential propagation across synapse[ms]
axon_delay = 0.05; %time delay for action potential propagation down axon [ms/um]
spont_frequency_exc = 0.3; %frequency of spontaneous APs, [Hz]
no_spontaneous_APs_exc = floor((t_end/1000)*spont_frequency_exc); %figure out the number of spontaneous APs
spont_frequency_inh = 0.3; %frequency of spontaneous APs, [Hz]
no_spontaneous_APs_inh = floor((t_end/1000)*spont_frequency_inh); %figure out the number of spontaneous APs

syn_decay = 0.975; %factor by which synaptic gmax decays per use
tau_decay_inh = 1000;  %time constant for synaptic gmax recovery [ms]
tau_decay_exc = 1000;  %time constant for synaptic gmax recovery [ms]

%% Retrieve .mat files containing geometry and connectivity maps
Conn_dev_aftercut = load('Demo_Aftercut.mat'); %0-1 connectivity matrix/map post-lesion sprout, indicating which neurons are connected 
Conn_dev_aftercut = Conn_dev_aftercut.Conn_dev_aftercut;
Conn_dev = load('Demo.mat'); %0-1 connectivity matrix/map of intact sphere, indicating which neurons are connected 
Conn_dev = Conn_dev.Conn_dev;
syn_strength_dev_aftercut = load('Demo_Aftercut_SynStrength.mat'); %matrix containg the synaptic strength post-lesion sprout
syn_strength_dev_aftercut = syn_strength_dev_aftercut.syn_strength_dev_aftercut;
syn_strength_dev = load('Demo_SynStrength.mat'); %matrix containg the synaptic strength of intact network
syn_strength_dev = syn_strength_dev.syn_strength_dev;
exc_neuron_indx = load('Demo_ExcIndx.mat'); %index of excitatory neurons in the intact network
exc_neuron_indx = exc_neuron_indx.exc_neuron_indx;
inh_neuron_indx = load('Demo_InhIndx.mat'); %index of inhibitory neurons in the intact network
inh_neuron_indx = inh_neuron_indx.inh_neuron_indx;
Dist = load('Demo_Dist.mat'); %arc distances between enurons
Dist = Dist.Dist;
V = load('Sphere_coordinates.mat');
radius = 950; % radius of the sphere in um
X = V.V(:,1)*radius;
Y = V.V(:,2)*radius;
Z = V.V(:,3)*radius;
st_node = 1; % Where the lesion was made? 1 indicates lesion from the northpole
No_rmnodes = size(Conn_dev,2)-size(Conn_dev_aftercut,2); % How many nodes were removed?
% To find the inhibitory and excitatory neurons indices in the network
% (NbyN matrix)
N = 2000; % Number of neurons
No_exc_conn = 300; % The setpoint value of connections on excitatory neurons
No_inh_conn = 140; % The setpoint value of connections on inhibitory neurons
dummy = No_exc_conn * ones(2000,1); dummy(inh_neuron_indx) = No_inh_conn; max_deg_out_rmv = dummy;
max_deg_out_rmv(st_node:No_rmnodes+st_node-1) = []; 
inh_neuron_indx_rmv = find(max_deg_out_rmv < No_exc_conn);
exc_neuron_indx_rmv = find(max_deg_out_rmv > No_inh_conn+1);
% To adjust for the node removal in the matrix of distances and coordinates
Dist_rmv = Dist; Dist_rmv(st_node:No_rmnodes+st_node-1,:) = []; Dist_rmv(:,st_node:No_rmnodes+st_node-1) = [];

X_rmv = X; Y_rmv = Y; Z_rmv = Z;
X_rmv(st_node:No_rmnodes+st_node-1) = []; % coordinates adjustment
Y_rmv(st_node:No_rmnodes+st_node-1) = [];
Z_rmv(st_node:No_rmnodes+st_node-1) = [];

%% DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect = 0:dt:t_end; %will hold vector of times
% Ask user what network to run
prompt = input('1.intact or 2.lesioned to run? Type the number');
if prompt == 1 % Intact
    neuron_postsynaptic_distance = Dist; % Arc distances between neurons
    Conn_dev_run = Conn_dev; % Connectivity matrix
    neuron_presynaptic_no = sum(Conn_dev_run,2); % Number of presynaptic connections to each neuron
    neuron_postsynaptic_no = sum(Conn_dev_run,1); % Number of postsynaptic connections for each neuron
    weighted_connectivity_map = syn_strength_dev; % Synaptic strength between any two nodes in a NbyN matrix
    inh_indx = inh_neuron_indx;
    exc_indx = exc_neuron_indx;
elseif prompt == 2 % Lesioned
    neuron_postsynaptic_distance = Dist_rmv;
    Conn_dev_run = Conn_dev_aftercut;
    neuron_presynaptic_no = sum(Conn_dev_run,2);
    neuron_postsynaptic_no = sum(Conn_dev_run,1);
    weighted_connectivity_map = syn_strength_dev_aftercut;
    N = N - No_rmnodes;
    inh_indx = inh_neuron_indx_rmv;
    exc_indx = exc_neuron_indx_rmv;
end

N_exc = length(exc_indx);
N_inh = length(inh_indx);
% max_syn_no = max(neuron_postsynaptic_no(:));
neuron_presynaptic_gmax = cell(N, 1);
neuron_v = zeros(N, length(t_vect));
neuron_g_exc = zeros(N, N); %initialize excitatory synaptic conductance to zero
neuron_g_inh = zeros(N, N); %initialize inhibitory synaptic conductance to zero 
neuron_v_plot = zeros(N, length(t_vect));
neuron_spont_APs_no_exc = no_spontaneous_APs_exc*ones(N, 1);
neuron_spont_APs_no_inh = no_spontaneous_APs_inh*ones(N, 1);
max_postsyn = max(neuron_postsynaptic_no(:));
neuron_presynaptic_APs = zeros(N, N, 200);
neuron_presynaptic_APs_ind = zeros(N, N);
excinh_check = zeros(N,1);
excinh_check(inh_indx)=1;
I = zeros(length(t_vect),N);

%% INTEGRATE THE EQUATION tau*dv/dt = -V + E_L + I_e*R_m
tic;
neuron_v(:, 1) = E_L; %First element of v at t=0
neuron_spont_APs_exc = sort(rand(N_exc, no_spontaneous_APs_exc) * t_end, 2, 'descend');  %Figure out the times of spontaneous APs and sort them based on their time of occurance   
neuron_spont_APs_inh = sort(rand(N_inh, no_spontaneous_APs_inh) * t_end, 2, 'descend');  %Figure out the times of spontaneous APs and sort them based on their time of occurance   
for i = 1 : N
    Presyn_conn = [];
    Presyn_conn = find(Conn_dev_run(:,i)>0); 
    for syn_no = 1 : neuron_postsynaptic_no(i)
        j = Presyn_conn(syn_no);
        if excinh_check(j) == 1
            neuron_presynaptic_gmax{i}(syn_no) = g_syn_max_inh; %Initialize all synapses to gmax 
        else
            neuron_presynaptic_gmax{i}(syn_no) = g_syn_max_exc; %Initialize all synapses to gmax 
        end
    end
end
for k = 1 : length(t_vect)
    if mod(k, 500) == 0
        k
    end
    t = k * dt;   
    neuron_g_next = zeros(N, N);
    for i = 1 : N       
        I_syn = 0; %Initialize total synaptic current at this t
        Presyn_conn = [];
        Presyn_conn = find(Conn_dev_run(:,i)>0); % Find all the presynaptic connections to neuron i
        for syn_no = 1 : neuron_postsynaptic_no(i) % Loop for all the presynaptic connections
            j = Presyn_conn(syn_no); % j is the presynaptic partner
            if excinh_check(j) == 1 % The presynaptic neuron is inhibitory
                I_syn = I_syn + (neuron_g_inh(j, i) * (neuron_v(i, k) - E_syn_inh)) * weighted_connectivity_map(j,i); % Calculate the synaptic current at time k, for synapse from neuron jth to neuron ith 
                if (neuron_presynaptic_APs_ind(j, i) > 0) %if an action potential is traveling on presynaptic axon
                    if (abs(t-neuron_presynaptic_APs(j, i, neuron_presynaptic_APs_ind(j, i)))<dt) %if AP in queue has a timestamp = current time
                        neuron_presynaptic_APs_ind(j, i) = neuron_presynaptic_APs_ind(j, i) - 1;
                        neuron_g_next(j, i) = neuron_presynaptic_gmax{i}(syn_no); %Fully activate synapse
                        neuron_presynaptic_gmax{i}(syn_no) = neuron_presynaptic_gmax{i}(syn_no) * syn_decay;  %Activation-dependent decay of gmax
                    else
                        neuron_g_next(j, i) = neuron_g_inh(j, i)*exp(-dt/tau_syn_inh); %Decay synapse
                    end
                end
                if neuron_presynaptic_gmax{i}(syn_no) < g_syn_max_inh %If synapse is still depressed
                    neuron_presynaptic_gmax{i}(syn_no) = neuron_presynaptic_gmax{i}(syn_no)*exp(dt/tau_decay_inh); %Recover synaptic gmax
                end  
                neuron_g_inh(j,i) = neuron_g_next(j, i);
            else  % The presynaptic neuron is excitatory 
                I_syn = I_syn + (neuron_g_exc(j, i) * (neuron_v(i, k) - E_syn_exc)) * weighted_connectivity_map(j,i); % Calculate the synaptic current at time k, for synapse from neuron jth to neuron ith 
                if (neuron_presynaptic_APs_ind(j, i) > 0) %if an action potential is traveling on presynaptic axon
                    if (abs(t-neuron_presynaptic_APs(j, i, neuron_presynaptic_APs_ind(j, i)))<dt) %if AP in queue has a timestamp = current time
                        neuron_presynaptic_APs_ind(j, i) = neuron_presynaptic_APs_ind(j, i) - 1;
                        neuron_g_next(j, i) = neuron_presynaptic_gmax{i}(syn_no); %Fully activate synapse
                        neuron_presynaptic_gmax{i}(syn_no) = neuron_presynaptic_gmax{i}(syn_no) * syn_decay;  %Activation-dependent decay of gmax
                    else
                        neuron_g_next(j, i) = neuron_g_exc(j, i)*exp(-dt/tau_syn_inh); %Decay synapse
                    end
                end
                if neuron_presynaptic_gmax{i}(syn_no) < g_syn_max_exc %If synapse is still depressed
                    neuron_presynaptic_gmax{i}(syn_no) = neuron_presynaptic_gmax{i}(syn_no)*exp(dt/tau_decay_exc); %Recover synaptic gmax
                end 
                neuron_g_exc(j,i) = neuron_g_next(j, i);
            end
            I(k,i) = I_syn;
        end 
%         I_k = neuron_g_k(i, k)*(neuron_v(i, k)-V_reset); %calculate K+ current to simulate recovery phase 
        if excinh_check(i) == 1 % if neuron i is inhibitory
            V_inf = E_L - I_syn*Rm_inh;
            neuron_v(i, k + 1) = V_inf + (neuron_v(i, k) - V_inf)*exp(-dt/tau_exc); %Calculate neuron potential
            spont_AP_flag = 0;
            ii = find(inh_indx==i);
           while (neuron_spont_APs_no_inh(ii) > 0) && (abs(t - neuron_spont_APs_inh(ii, neuron_spont_APs_no_inh(ii))) < dt)
                spont_AP_flag = 1;
                neuron_spont_APs_no_inh(ii) = neuron_spont_APs_no_inh(ii) - 1;
           end
           if (neuron_v(i, k) > V_th) || spont_AP_flag %Cell spiked
                neuron_v(i, k + 1) = V_reset; %Set voltage back to V_reset
                neuron_v_plot(i, k) = V_spike; %Plot spike
                %neuron_g_k(t+1) = g_k_max;  %activate K+ conductance
                Postsyn_conn = [];
                Postsyn_conn = find(Conn_dev_run(i,:)>0);
                for n_post = 1 : neuron_presynaptic_no(i)
                    j = Postsyn_conn(n_post);   
                    delay = AP_delay + axon_delay * neuron_postsynaptic_distance(i, j); %AP delay is based on axon propagation and synaptic delay
                    neuron_presynaptic_APs(i, j, (2 : neuron_presynaptic_APs_ind(i, j) + 1)) = neuron_presynaptic_APs(i, j, (1 : neuron_presynaptic_APs_ind(i, j))); %Shift first row to the next, to empty the first row 
                    neuron_presynaptic_APs_ind(i, j) = neuron_presynaptic_APs_ind(i, j) + 1;
                    neuron_presynaptic_APs(i, j, 1)= t + delay; %Pass AP time to postsynaptic neuron, (in the first row)
                end
            else
            neuron_v_plot(i, k) = neuron_v(i, k); %Plot the actual voltage if no spike
    %         neuron_g_k(i, k) = neuron_g_k(i, k)*exp(-dt/tau_k);
           end 
        else % if i is excitatory
            V_inf = E_L - I_syn*Rm_exc;
            neuron_v(i, k + 1) = V_inf + (neuron_v(i, k) - V_inf)*exp(-dt/tau_exc); %Calculate neuron potential
            spont_AP_flag = 0;
            ii = find(exc_indx==i);
           while (neuron_spont_APs_no_exc(ii) > 0) && (abs(t - neuron_spont_APs_exc(ii, neuron_spont_APs_no_exc(ii))) < dt)
                spont_AP_flag = 1;
                neuron_spont_APs_no_exc(ii) = neuron_spont_APs_no_exc(ii) - 1;
           end
           if (neuron_v(i, k) > V_th) || spont_AP_flag %Cell spiked
                neuron_v(i, k + 1) = V_reset; %Set voltage back to V_reset
                neuron_v_plot(i, k) = V_spike; %Plot spike
                %neuron_g_k(t+1) = g_k_max;  %activate K+ conductance
                Postsyn_conn = [];
                Postsyn_conn = find(Conn_dev_run(i,:)>0);
                for n_post = 1 : neuron_presynaptic_no(i)
                    j = Postsyn_conn(n_post);   
                    delay = AP_delay + axon_delay * neuron_postsynaptic_distance(i, j); %AP delay is based on axon propagation and synaptic delay
                    neuron_presynaptic_APs(i, j, (2 : neuron_presynaptic_APs_ind(i, j) + 1)) = neuron_presynaptic_APs(i, j, (1 : neuron_presynaptic_APs_ind(i, j))); %Shift first row to the next, to empty the first row 
                    neuron_presynaptic_APs_ind(i, j) = neuron_presynaptic_APs_ind(i, j) + 1;
                    neuron_presynaptic_APs(i, j, 1)= t + delay; %Pass AP time to postsynaptic neuron, (in the first row)
                end
            else
            neuron_v_plot(i, k) = neuron_v(i, k); %Plot the actual voltage if no spike
    %         neuron_g_k(i, k) = neuron_g_k(i, k)*exp(-dt/tau_k);
           end 
        end
    end
end
toc

%% Make plots:
figure('color','white');
plot(t_vect, squeeze(neuron_v_plot(1, :)), 'Color', [0 0 0]);
title('Neuron1');
xlabel('Time [ms]'); 
% % xlim([340 460]);
ylabel('Voltage [mV]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
box off; 
% subplot(4, 2, 3); plot(t_vect, squeeze(I(:,1)), 'b');
% xlabel('Time [ms]'); 
% % xlim([340 460]);
% ylabel('Current [nA]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off; 
% % figure('color','white');
% subplot(4, 2, 2);plot(t_vect, squeeze(neuron_v_plot(500, :)), 'Color', [0 0 0]);
% title('Neuron500');
% % xlim([340 460]);
% xlabel('Time [ms]');
% ylabel('Voltage [mV]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off;  
% subplot(4, 2, 4); plot(t_vect, squeeze(I(:,500)), 'b');
% xlabel('Time [ms]');
% % xlim([340 460]);
% ylabel('Current [nA]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off; 
% % figure('color','white');
% subplot(4, 2, 5);plot(t_vect, squeeze(neuron_v_plot(1000, :)), 'Color', [0 0 0]);
% title('neuron1000');
% xlabel('Time [ms]');
% % xlim([340 460]);
% ylabel('Voltage [mV]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off;
% subplot(4, 2, 7); plot(t_vect, squeeze(I(:,1000)), 'b');
% xlabel('Time [ms]');
% % xlim([340 460]);
% ylabel('Current [nA]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off; 

% figure('color','white');
% subplot(4, 2, 6);plot(t_vect, squeeze(neuron_v_plot(1500, :)), 'Color', [0 0 0]);
% title('neuron1500');
% xlabel('Time [ms]');
% % xlim([340 460]);
% ylabel('Voltage [mV]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off;
% subplot(4, 2, 8); plot(t_vect, squeeze(I(:,1500)), 'b');
% xlabel('Time [ms]');
% % xlim([340 460]);
% ylabel('Current [nA]');set(gca, 'fontsize', 16); set(gca,'TickDir', 'out');
% box off; 
if prompt == 1
    x = X;
    y = Y;
    z = Z;
elseif prompt == 2
    x = X_rmv;
    y = Y_rmv;
    z = Z_rmv;
end
max_Z = max(z(:)); %Find Z of the highest neuron
max_Z_indx = find(z == max_Z);
max_Z_indx(max_Z_indx == 0) = [];
edge_point = [x(max_Z_indx),y(max_Z_indx),z(max_Z_indx)]; %Coordinates of the neuron located at the highest Z, closest to the lesion area
Z_edge_distance = edge_point(3) - z;
[sorted_z,sorted_z_indx] = sort(Z_edge_distance);
% figure('color','white'); 
% imagesc(neuron_v_plot(sorted_z_indx,:),[-90 20]);
% colormap jet;colorbar; 

V_raster = neuron_v_plot;
% V_raster(1:length(inh_indx),:) = neuron_v_plot(inh_indx,:);
% V_raster(length(inh_indx)+1:end,:) = neuron_v_plot(exc_indx,:);
V_raster_plot = V_raster(:,1:end) > -50;
figure('color','white'); 
[yy,xx] = find(V_raster_plot==1);
% [y,x] = find(V_raster_plot(1:length(inh_indx),:)==1);
scatter(xx,yy,'k.','linewidth',1); hold on;
title('raster plot');
% scatter(x,y,'r.','linewidth',1); hold on;
xlabel('time(ms)');
ylabel('neurons'); 
% set(gca,'ytick',[])
set(gca, 'fontsize', 24); 
box off; set(gca,'TickDir','out','fontsize', 24, 'linewidth', 1.2);
% legend('inh','exc'); legend boxoff;
% set(gcf,'renderer','Painters');


