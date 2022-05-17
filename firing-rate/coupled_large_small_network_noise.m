close all;
clear all;
%Code used to generate Figure 2 (C, G, H, J, K).  
%Different results may be produced every time the code is run due to
%use of Poisson noise and stochastic network generation.   Probability in Fig. 2(C) was
%determined with n = 3 simulations.  Probabilities in Fig. 2(H) and
%Fig. 2(J) were determined using n = 10 simulations.   
%Simulations for isolated networks can be obtained by setting
%connection_weight_to_replace_array = 0
%Simulations for coupled networks can be obtained by setting
%connection_weight_to_replace_array = 200
%Noise Amplitude can be varied by setting stimulation_amplitude_array =
%desired noise amplitude

%%
%Parameter definition
tau_m = 10e-3; %membrane time constant [sec] 
t_end = 5; %end time for the stimulation [sec]
dt = 2e-5; %time step [sec]
threshold = 0.1; %input threshold for piece-wise linear neuron to fire
exc_slope = 0.01; %increase in firing rate with increasing input
inh_slope = exc_slope*4.3;
max_rate = 20; %maximum firing rate [Hz]
stochastic_input_scale = 0.01;
num_exc_neurons_large = 1000; %number of excitatory neurons in large network
num_inh_neurons_large = 100;  %number of inhibitory neurons in large network
num_exc_neurons_small = 100; %number of excitatory neurons in small network
num_inh_neurons_small = 10;  %number of inhibitory neurons in small network

number_of_neurons_large = num_exc_neurons_large+num_inh_neurons_large; %maximum number of neurons to not run out of memory is 40,000
max_exc_input = 150;
max_inh_input = 10;

K_exc = 1;
K_inh = 1;

stimulation_amplitude_array = [max_rate/400]; %amplitude of the Poisson noise events in individual neurons
stim_freq = 0.5; %average frequency of Poisson events [Hz]
connection_weight_to_replace_array = 200; %number of connections replaced to create links between networks

%%
%%Using different large-small connectivity strengths, obtain stimulation
%%threshold for generating max rate firing.
distance_between_networks = (num_exc_neurons_small+num_inh_neurons_small)/(2*pi)*0.5; %equal to a proportion of the radius of the small network
num_simulations_per_conn_weight = 1; %number of simulations to run per "number of connections replaced" to account for different network conns, stimulated neurons
smallest_stim_max_rate_large = zeros(length(connection_weight_to_replace_array), num_simulations_per_conn_weight); %stim threshold results for large network
smallest_stim_max_rate_small = zeros(length(connection_weight_to_replace_array), num_simulations_per_conn_weight); %stim threshold results for small network

tic
for conn_index = 1:length(connection_weight_to_replace_array)
    connection_weight_to_replace = connection_weight_to_replace_array(conn_index);
    for simulation_index = 1:num_simulations_per_conn_weight
        
    

%% Create Large Network
radius_large = number_of_neurons_large/(2*pi);
dtheta = 2*pi/number_of_neurons_large;
neuron_theta_large = 0:dtheta:(2*pi-dtheta);
neuron_type_large = ones(1,length(neuron_theta_large));
for i=1:round(number_of_neurons_large/num_inh_neurons_large):length(neuron_theta_large)
    neuron_type_large(i) = -1;  %inhibitory neurons marked with '-1', excitatory with '1'
end

distance_between_neighbors = radius_large*dtheta;
lambda_exc = 2000*distance_between_neighbors;
lambda_inh = 16*distance_between_neighbors;

W_large = zeros(number_of_neurons_large, number_of_neurons_large);

%stochastic connectivity

for iter = 1:1000
    disp(['Network builder iteration = ', num2str(iter)]);
    for postsyn_neuron = 1:number_of_neurons_large
        set_of_presyn_neurons = datasample (1:number_of_neurons_large,max_exc_input/10);
        for index = 1:length(set_of_presyn_neurons)
            presyn_neuron = set_of_presyn_neurons(index);
            
            distance1 = abs((neuron_theta_large(postsyn_neuron) - neuron_theta_large(presyn_neuron)))*radius_large;
            distance2 = abs(2*pi - (neuron_theta_large(postsyn_neuron) - neuron_theta_large(presyn_neuron)))*radius_large;
            if distance1 > distance2 %pick the shortest distance around the circumference
                distance = distance2;
            else
                distance = distance1;
            end
            if neuron_type_large(presyn_neuron) == 1 && sum(W_large(postsyn_neuron,neuron_type_large == 1)) < max_exc_input
                if rand(1)<(K_exc*exp(-distance/lambda_exc)) && presyn_neuron ~= postsyn_neuron %do not make synapses to self
                    W_large(postsyn_neuron,presyn_neuron) = W_large(postsyn_neuron,presyn_neuron) + 1;
                end
            elseif neuron_type_large(presyn_neuron) == -1 && abs(sum(W_large(postsyn_neuron,neuron_type_large == -1))) < max_inh_input
                if  rand(1)<(K_inh*exp(-distance/lambda_inh)) && neuron_type_large(postsyn_neuron) ~= -1   %do not make inhibitory-inhbitory synapses
                    W_large(postsyn_neuron,presyn_neuron) = W_large(postsyn_neuron,presyn_neuron) - 1;
                end
            end
        end
    end
end

x_coords_large = radius_large*sin(neuron_theta_large);
y_coords_large = radius_large*cos(neuron_theta_large);


W_large = W_large - diag(diag(W_large)); %set neuron self-connectivity to zero

%% Create Small Network
number_of_neurons_small = num_exc_neurons_small+num_inh_neurons_small; %maximum number of neurons to not run out of memory is 40,000

radius_small = number_of_neurons_small/(2*pi);
dtheta = 2*pi/number_of_neurons_small;
neuron_theta_small = 0:dtheta:(2*pi-dtheta);
neuron_type_small = ones(1,length(neuron_theta_small));
for i=1:round(number_of_neurons_small/num_inh_neurons_small):length(neuron_theta_small)
    neuron_type_small(i) = -1;  %inhibitory neurons marked with '-1', excitatory with '1'
end

W_small = zeros(number_of_neurons_small, number_of_neurons_small);

%stochastic connectivity

for iter = 1:1000
    disp(['Network builder iteration = ', num2str(iter)]);
    for postsyn_neuron = 1:number_of_neurons_small
        set_of_presyn_neurons = datasample (1:number_of_neurons_small,max_exc_input/10);
        for index = 1:length(set_of_presyn_neurons)
            presyn_neuron = set_of_presyn_neurons(index);
            distance1 = abs(neuron_theta_small(postsyn_neuron) - neuron_theta_small(presyn_neuron))*radius_small;
            distance2 = abs(2*pi - (neuron_theta_small(postsyn_neuron) - neuron_theta_small(presyn_neuron)))*radius_small;
            if distance1 > distance2
                distance = distance2;
            else
                distance = distance1;
            end
            if neuron_type_small(presyn_neuron) == 1 && sum(W_small(postsyn_neuron,neuron_type_small == 1)) < max_exc_input
                if rand(1)<(K_exc*exp(-distance/lambda_exc)) && presyn_neuron ~= postsyn_neuron %do not make synapses to self
                    W_small(postsyn_neuron,presyn_neuron) = W_small(postsyn_neuron,presyn_neuron) + 1;
                end
            elseif neuron_type_small(presyn_neuron) == -1 && abs(sum(W_small(postsyn_neuron,neuron_type_small == -1))) < max_inh_input
                if  rand(1)<(K_inh*exp(-distance/lambda_inh)) && neuron_type_small(postsyn_neuron) ~= -1   %do not make inhibitory-inhbitory synapses
                    W_small(postsyn_neuron,presyn_neuron) = W_small(postsyn_neuron,presyn_neuron) - 1;
                end
            end
        end
    end
end

x_coords_small = radius_small*sin(neuron_theta_small);
y_coords_small = radius_small*cos(neuron_theta_small);


time_vec = 0:dt:t_end;

W_small = W_small - diag(diag(W_small)); %set neuron self-connectivity to zero

%% Create combined network
number_of_neurons = number_of_neurons_large + number_of_neurons_small;
neuron_type = [neuron_type_large neuron_type_small];
W = zeros (number_of_neurons, number_of_neurons);
W(1:number_of_neurons_large, 1:number_of_neurons_large) = W_large;
W(number_of_neurons_large+1:end, number_of_neurons_large+1:end) = W_small;

W_unconnected = W;

%replace old connections in large network with new connections to small
%network



% replace connections in the large network
connection_weight_replaced_large = 0;
iteration_large = 1;
while connection_weight_replaced_large < connection_weight_to_replace && iteration_large < connection_weight_to_replace*75
    neuron_no = ceil(number_of_neurons_large*rand(1));  %pick a random neuron in the large network
    synapse_no = ceil(number_of_neurons_large*rand(1));  %pick a random synapse in the large network
    weight = W(neuron_no, synapse_no);
    if weight > 0  %if a synapse exists and is not inhibitory
        pre_syn_neuron_no = ceil((number_of_neurons_small)*rand(1))+number_of_neurons_large;
        %networks are position such that neuron_theta_large = 0 and
        %neuron_theta_small = pi are closest points between two networks
        %distance is calculated as path from neuron_no along the
        %circumference of large network to neuron_theta_large = 0, then
        %across the gap between networks, and then as a path along the
        %circumference of small network starting at neuron_theta_small = pi
        %and then to the pre_syn_neuron_no
        distance_large1 = neuron_theta_large(neuron_no)*radius_large;
        distance_large2 = abs((2*pi - neuron_theta_large(neuron_no)))*radius_large;
        if distance_large1 > distance_large2 %pick smallest path around circumference of large network
            distance_large = distance_large2;
        else
            distance_large = distance_large1;
        end
        distance_small = abs(neuron_theta_small(pre_syn_neuron_no-number_of_neurons_large) - pi)*radius_small;
        distance = distance_large + distance_between_networks + distance_small;
        if neuron_type(pre_syn_neuron_no) ~= -1 && rand(1)<(K_exc*exp(-distance/lambda_exc))
            W(neuron_no, pre_syn_neuron_no) = W(neuron_no,pre_syn_neuron_no) + weight; %add weight to synapse if presynaptic neuron is excitatory
            W(neuron_no, synapse_no) = 0; %eliminate the synapse in the large network
            connection_weight_replaced_large = connection_weight_replaced_large+weight;
        end
    end
    iteration_large = iteration_large+1;
end

% replace connections in the small network
connection_weight_replaced_small = 0;
iteration = 1;

while connection_weight_replaced_small < connection_weight_to_replace && iteration < connection_weight_to_replace*10
    neuron_no = ceil(number_of_neurons_small*rand(1))+number_of_neurons_large;  %pick a random neuron in the small network
    synapse_no = ceil(number_of_neurons_small*rand(1))+number_of_neurons_large;  %pick a random synapse in the small network
    weight = W(neuron_no, synapse_no);
    if weight > 0  %if a synapse exists and is not inhibitory
        pre_syn_neuron_no = ceil((number_of_neurons_large)*rand(1));
        %networks are position such that neuron_theta_large = 0 and
        %neuron_theta_small = pi are closest points between two networks
        %distance is calculated as path from neuron_no along the
        %circumference of the small network to neuron_theta_small = pi, then
        %across the gap between networks, and then as a path along the
        %circumference of large network starting at neuron_theta_large = 0
        %and then to the pre_syn_neuron_no
        distance_small = (neuron_theta_small(neuron_no-number_of_neurons_large)-pi)*radius_small;
        distance_large1 = neuron_theta_large(pre_syn_neuron_no)*radius_large;
        distance_large2 = abs(2*pi - neuron_theta_large(pre_syn_neuron_no))*radius_large;
        if distance_large1 > distance_large2 %pick smallest path around circumference of large network
            distance_large = distance_large2;
        else
            distance_large = distance_large1;
        end
        distance = distance_small + distance_between_networks + distance_large;
        if neuron_type(pre_syn_neuron_no) ~= -1 
            W(neuron_no, pre_syn_neuron_no) = W(neuron_no, pre_syn_neuron_no) + weight; %replace the synapse if presynaptic neuron is excitatory
            W(neuron_no, synapse_no) = 0; %eliminate the synapse in the small network
            connection_weight_replaced_small = connection_weight_replaced_small+weight;
        end
    end
    iteration = iteration+1;
end
        
%% Plot neuron positions and connections of the combined network
y_coords_small = y_coords_small+(max(y_coords_large)-min(y_coords_small))+distance_between_networks; %move the small network to the top of large one, with margin

x_coords = [x_coords_large x_coords_small];
y_coords = [y_coords_large y_coords_small];
plot(x_coords(neuron_type == 1),y_coords(neuron_type ==1),'b^','MarkerSize',3, 'MarkerFaceColor', 'b');
hold on;
plot(x_coords(neuron_type ==-1),y_coords(neuron_type ==-1), 'ro','MarkerSize',3, 'MarkerFaceColor', 'r');

axis('manual');
for neur_shown = 1:ceil(number_of_neurons/2) %show pre-synapses for every 2nd neuron
    neur_index = ceil (rand(1)*number_of_neurons);
    if neuron_type(neur_index) == 1
        plotted_synapses_per_neuron = ceil(number_of_neurons/50); %show every 50th pre-synaptic connection
    else
        plotted_synapses_per_neuron = ceil(number_of_neurons/20);
    end
    for synapses_shown = 1:plotted_synapses_per_neuron 
        synapse_index = ceil(rand(1)*number_of_neurons);
        if W(neur_index,synapse_index) ~=0
            start_xy = [x_coords(synapse_index) y_coords(synapse_index)];
            end_xy = [x_coords(neur_index) y_coords(neur_index)];
            if W(neur_index,synapse_index) >= 1  %excitatory connection
                line ([start_xy(1) end_xy(1)], [start_xy(2) end_xy(2)], 'Color', 'blue','LineWidth',abs(W(neur_index,synapse_index))/2);               
            elseif W(neur_index,synapse_index) <= 1  %inhibitory connection
                line ([start_xy(1) end_xy(1)], [start_xy(2) end_xy(2)], 'Color', 'red','LineWidth',abs(W(neur_index,synapse_index))/2);
            end
        end
    end
end

hold off;

%% Simulate

    rates = zeros(size(W,1), length(time_vec), length(stimulation_amplitude_array));
    stimulation_index = 1;
    for stimulation_amplitude_index = 1:length(stimulation_amplitude_array)
        r = zeros(size(W,1), length(time_vec));
        stimulation = zeros(size(r));
        neuron_indices = 1:number_of_neurons;
        exc_neuron_indices = neuron_indices(neuron_type == 1);
        stimulation_start = 0; %seconds
        stimulation_end = t_end; %seconds
        total_stim_number = ceil(stim_freq*(stimulation_end - stimulation_start)); %number of random noise stimulations per neuron per simulation time
        for stim_neuron_index = 1:length(exc_neuron_indices)  %create Poisson noise events for each neuron
            for stimulation_no = 1:total_stim_number
                stimulation_time = ceil(rand(1,1)*(stimulation_end - stimulation_start)/dt+stimulation_start/dt); 
                stimulation(exc_neuron_indices(stim_neuron_index),stimulation_time) = stimulation_amplitude_array(stimulation_amplitude_index); %assign stim amplitude to stimulation times per neuron
            end
        end
        for i = 1:length(time_vec)-1
            input = W*r(:,i);
            r(neuron_type==1,i+1) = dt/tau_m*(-r(neuron_type==1,i) + response_function_arr(input(neuron_type==1), threshold, exc_slope, max_rate))+r(neuron_type==1,i)+stimulation(neuron_type==1,i);
            r(neuron_type==-1,i+1) = dt/tau_m*(-r(neuron_type==-1,i) + response_function_arr(input(neuron_type==-1), threshold, inh_slope, max_rate))+r(neuron_type==-1,i)+stimulation(neuron_type==-1,i);
        end
    disp(['For stimulation amplitude = ' num2str(stimulation_amplitude_array(stimulation_amplitude_index)) ', the max firing rate was ' num2str(max(max(r(1:number_of_neurons_large,:))))]);
    rates(:, :, stimulation_amplitude_index) = r; %collect firing rates for simulations with different numbers of stimulated neurons
    
    excitatory_rate = rates(neuron_type == 1, :, stimulation_amplitude_index);
    inhibitory_rate = rates(neuron_type == -1, :, stimulation_amplitude_index);
    figure;imagesc(excitatory_rate, [0 0.1]);
    %figure;imagesc(excitatory_rate, [0 20]);
    xlabel('Time (timestep)');
    ylabel('Excitatory neuron number');
    title(['Stimulation amplitude = ' num2str(stimulation_amplitude_array(stimulation_amplitude_index))]);

    figure;imagesc(inhibitory_rate);
    xlabel('Time(timestep)');
    ylabel('Inhibitory neuron number');
    title(['Stimulation amplitude = ' num2str(stimulation_amplitude_array(stimulation_amplitude_index))]);
    
    figure;plot(time_vec,sum(excitatory_rate(1:num_exc_neurons_large,:))/num_exc_neurons_large);
    hold on;
    plot(time_vec,sum(inhibitory_rate(1:num_inh_neurons_large,:))/num_inh_neurons_large);
    hold off;
    xlabel('Time(sec)');
    ylabel('Firing Rate');
    title('Large Network');
    
    figure;plot(time_vec,sum(excitatory_rate((num_exc_neurons_large+1):end,:))/num_exc_neurons_small);
    hold on;
    plot(time_vec,sum(inhibitory_rate((num_inh_neurons_large+1):end,:))/num_inh_neurons_small);
    hold off;
    xlabel('Time(sec)');
    ylabel('Firing Rate');
    title('Small Network');
    end
end
end
toc


%analysis of probability of E-E or E-I connection
for i = 1:number_of_neurons_large
   probs_large(i) = sum(W_large(i,:)>0)/number_of_neurons_large;
end

disp(['Average probability of E->E or E->I connection in Large network is ' num2str(mean(probs_large))]);

for i = 1:number_of_neurons_small
   probs_small(i) = sum(W_small(i,:)>0)/number_of_neurons_small;
end

disp(['Average probability of E->E or E->I connection in Small network is ' num2str(mean(probs_small))]);







