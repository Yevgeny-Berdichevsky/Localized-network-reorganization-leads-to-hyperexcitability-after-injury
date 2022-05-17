close all;
clear all;
%Code used to generate Figure 2 (D, F).  
%Different results may be produced every time the code is run due to stochastic nature of network generation.  
%Probability in Fig. 2(E) was
%determined with n = 7 per network size.  Network size can be changed by altering num_exc_neurons (defalt = 100)
%and num_inh_neurons (default = 10), while keeping
%num_exc_neurons/num_inh_neurons = 10
%%
%Parameter Definition
tau_m = 10e-3; %membrane time constant [sec]
t_end = 0.1; %end time for the stimulation [sec]
dt = 2e-5; %time step [sec]
threshold = 0.1; %input threshold for piece-wise linear neuron to fire
exc_slope = 0.01; %increase in firing rate with increasing input
inh_slope = exc_slope*4.3;
max_rate = 20; %maximum firing rate [Hz]
stochastic_input_scale = 0.01;
num_exc_neurons = 100; %number of excitatory neurons
num_inh_neurons = 10;  %number of inhibitory neurons
number_of_neurons = num_exc_neurons+num_inh_neurons; %maximum number of neurons to not run out of memory is 40,000
max_exc_input = 150;
max_inh_input = 10;

K_exc = 1;
K_inh = 1;

stimulation_amplitude = max_rate/400;

%%
%Create Network
radius = number_of_neurons/(2*pi);
dtheta = 2*pi/number_of_neurons;
neuron_theta = 0:dtheta:(2*pi-dtheta);
neuron_type = ones(1,length(neuron_theta));
for i=1:round(number_of_neurons/num_inh_neurons):length(neuron_theta)
    neuron_type(i) = -1;  %inhibitory neurons marked with '-1', excitatory with '1'
end

distance_between_neighbors = radius*dtheta;
lambda_exc = 2000*distance_between_neighbors;
lambda_inh = 16*distance_between_neighbors;


W = zeros(number_of_neurons, number_of_neurons);

%stochastic connectivity

for iter = 1:1000
    disp(['Network builder iteration = ', num2str(iter)]);
    for postsyn_neuron = 1:number_of_neurons
        set_of_presyn_neurons = datasample (1:number_of_neurons,max_exc_input/10);
        for index = 1:length(set_of_presyn_neurons)
            presyn_neuron = set_of_presyn_neurons(index);
            distance = abs((neuron_theta(postsyn_neuron) - neuron_theta(presyn_neuron))*radius);
            if neuron_type(presyn_neuron) == 1 && sum(W(postsyn_neuron,neuron_type == 1)) < max_exc_input
                if rand(1)<(K_exc*exp(-distance/lambda_exc)) && presyn_neuron ~= postsyn_neuron %do not make synapses to self
                    W(postsyn_neuron,presyn_neuron) = W(postsyn_neuron,presyn_neuron) + 1;
                end
            elseif neuron_type(presyn_neuron) == -1 && abs(sum(W(postsyn_neuron,neuron_type == -1))) < max_inh_input
                if  rand(1)<(K_inh*exp(-distance/lambda_inh)) && neuron_type(postsyn_neuron) ~= -1   %do not make inhibitory-inhbitory synapses
                    W(postsyn_neuron,presyn_neuron) = W(postsyn_neuron,presyn_neuron) - 1;
                end
            end
        end
    end
end

x_coords = radius*sin(neuron_theta);
y_coords = radius*cos(neuron_theta);

time_vec = 0:dt:t_end;

W = W - diag(diag(W)); %set neuron self-connectivity to zero

%%
%Run simulation
tic
stimulated_neuron_min = 1;
stimulated_neuron_max = 10;
rates = zeros(size(W,1), length(time_vec), stimulated_neuron_max - stimulated_neuron_min);
for stimulated_neuron_no = stimulated_neuron_min:stimulated_neuron_max
    r = zeros(size(W,1), length(time_vec));
    stimulation = zeros(size(r));
    neuron_indices = 1:number_of_neurons;
    exc_neuron_indices = neuron_indices(neuron_type == 1);
    stimulated_neurons = datasample (exc_neuron_indices, stimulated_neuron_no, 'Replace', false); %pick a number of excitatory neurons to stimulate
    stimulation_timestep = 100;
    stimulation(stimulated_neurons,stimulation_timestep) = stimulation_amplitude; 
    for i = 1:length(time_vec)-1
        input = W*r(:,i);
        r(neuron_type==1,i+1) = dt/tau_m*(-r(neuron_type==1,i) + response_function_arr(input(neuron_type==1), threshold, exc_slope, max_rate))+r(neuron_type==1,i)+stimulation(neuron_type==1,i);
        r(neuron_type==-1,i+1) = dt/tau_m*(-r(neuron_type==-1,i) + response_function_arr(input(neuron_type==-1), threshold, inh_slope, max_rate))+r(neuron_type==-1,i)+stimulation(neuron_type==-1,i);
    end
    disp(['For stimulated_neuron_no = ' num2str(stimulated_neuron_no) ', the max firing rate was ' num2str(max(max(r)))]);
    rates(:, :, stimulated_neuron_no) = r; %collect firing rates for simulations with different numbers of stimulated neurons
end

toc

%%
%Generate plots for 2 and 4 stimulated neurons.
for selected_no = 2:2:4
    excitatory_rate = rates (neuron_type == 1, :, selected_no);
    inhibitory_rate = rates (neuron_type == -1, :, selected_no);
    figure;imagesc(excitatory_rate);
    xlabel('Time (timestep)');
    ylabel('Excitatory Neuron Number');
    title(['Number of stimulated neurons = ' num2str(selected_no)]);
    exc_average = mean(excitatory_rate);
    inh_average = mean(inhibitory_rate);
    figure;plot(time_vec, exc_average);
    hold on;
    plot(time_vec, inh_average);
    xlabel('Time (sec)');
    ylabel('Firing Rate');
    legend('Excitatory Firing Rate', 'Inhibitory Firing Rate');
    title(['Number of stimulated neurons = ' num2str(selected_no)]);
    hold off;
end

%Generate system trajectories for 1 to 10 stimulated neurons
figure;hold on;
for selected_no = stimulated_neuron_min:stimulated_neuron_max
    excitatory_rate = rates(neuron_type == 1, :, selected_no);
    inhibitory_rate = rates(neuron_type == -1, :, selected_no);
    exc_average = mean(excitatory_rate);
    inh_average = mean(inhibitory_rate);
    plot(inh_average, exc_average);
end
xlabel('Inhibitory Firing Rate');
ylabel('Excitatory Firing Rate');
hold off;
xlim([0 6]);
ylim([0 1]);






