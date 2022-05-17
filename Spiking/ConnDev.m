function [Conn, syn_strength] = ...
    ConnDev(N,Loop_No,max_deg_out,max_deg_in,Conn_Probability,sig_mean)
% This function makes all the connections between neurons based on the set
% points and probability distributions 
% Input list:
% 'N' = Number of nodes
% 'Loop_No' = number of loops to run for creating connections
% 'max_deg_out' = matrix conatins the set points of sending connections for exc. and inh. neurons
% 'max_deg_in' = matrix conatins the set points of receiving connections for exc. and inh. neurons
% 'Conn_Probability' = matrix contains probability of connections between any two nodes
% based on their arc distances
% 'Conn' = Zero NbyN connectivity matrix
% 'Syn_strength' = Zero NbyN synaptic strength matrix
% 'sig_mean' = Mean value for sigmoid function that constrains the nodes ability
%to accumulate connections too rapidly

% Output list:
% 'Conn' = Filled NbyN connectivity matrix
% 'Syn_strength' = Filled NbyN synaptic strength matrix
Conn = zeros(N); % Connectivity matrix is a NbyN zero matrix to begin with
syn_strength = zeros(N); % Synaptic strength matrix is a NbyN zero matrix to begin with
connections = 0 : N;
Sig = 1./(1 + exp(0.005.*(connections'-sig_mean))); % Sigmoid function that ...
... dictates the probability of a neuron being chosen based on the number of connections it has
sigmoid_func = fittype('1./(1 + exp(-a.*(x-c)))');
inh_neuron_indx = find(max_deg_out < 300); % Indices of inh. neurons
exc_neuron_indx = find(max_deg_out > 141); % Indices of exc. neurons
max_degin_inhinh = 10; % Number of inhibitory connections that inhibitory neurons receive
max_degin_inhexc = 18; % Number of inhibitory connections that excitatory neurons receive
max_degin_excinh = 130; % Number of excitatory connections that inhibitory neurons receive
max_degin_excexc = 282; % Number of excitatory connections that inhibitory neurons receive
eta = 1000; %500
connections = 0 : N;
Pick_Prob = fit(connections',Sig, sigmoid_func,'start',[-0.005 sig_mean]); 
for hours = 1 : Loop_No
    True = sum(syn_strength,2) < max_deg_out;
    Check_out = find(True == 1);
    True_in = sum(syn_strength,1)' < max_deg_in;
    Check_in = find(True_in == 1);
    if isempty(Check_out) == 1 || isempty(Check_in) == 1
        break
    else
        syn_strength_dummy = syn_strength;
        syn_strength_dummy(inh_neuron_indx,:) = syn_strength_dummy(inh_neuron_indx,:) .* 1;
        no_axo = sum(syn_strength_dummy(Check_out,:),2);
        if length(Check_out) <= eta
            eta = length(Check_out);
        end
        if length(Check_out) == 1
            randomi = Check_out;
        else
            randomi = randsample(Check_out,eta,'true',Pick_Prob(no_axo));
            dist_prob = Conn_Probability(randomi,Check_in);
            syn_strength_dummy2 = syn_strength;
            syn_strength_dummy2(:,inh_neuron_indx) = syn_strength_dummy2(:,inh_neuron_indx) .* 3;
            no_den = sum(syn_strength_dummy2(:,Check_in));
            avail_incap = Pick_Prob(no_den);
            incap = repmat(avail_incap',size(dist_prob,1),1);
            tot_prob = dist_prob .* incap;
            for i = 1:size(tot_prob,1)
                if sum(syn_strength(randomi(i),:)) < max_deg_out(randomi(i))
                    if length(Check_in) == 1
                        j = Check_in;
                    else
                        j = randsample(Check_in,1,'true',tot_prob(i,:));
                        if sum(syn_strength(:,j)) < max_deg_in(j) 
                            check_inh_in = ismember(j,inh_neuron_indx);
                            check_inh_out = ismember(randomi(i),inh_neuron_indx);
                            all_inh_input = sum(syn_strength(inh_neuron_indx,j));
                            all_exc_input = sum(syn_strength(exc_neuron_indx,j));
                            if check_inh_out == 1 % Inhibitory neuron sending connection
                                if check_inh_in == 1 % Inhibitory neuron to inhibitory neuron (check for inh degree in of inh neuron)
                                    if all_inh_input < max_degin_inhinh
                                    Conn(randomi(i),j) = 1;
                                    syn_strength(randomi(i),j) = syn_strength(randomi(i),j) + 1;
                                    end
                                elseif check_inh_in == 0 % Inhibitory neuron to excitatory neuron (check for inh degree in of exc neuron)
                                    if all_inh_input < max_degin_inhexc
                                    Conn(randomi(i),j) = 1;
                                    syn_strength(randomi(i),j) = syn_strength(randomi(i),j) + 1; 
                                    end
                                end
                            else % Excitatory neuron sending connection
                                if check_inh_in == 1 % Excitatory onto inhibitory neuron
                                    if all_exc_input < max_degin_excinh
                                    Conn(randomi(i),j) = 1;
                                    syn_strength(randomi(i),j) = syn_strength(randomi(i),j) + 1; 
                                    end
                                elseif check_inh_in == 0 % Excitatory onto excitatory neuron
                                    if all_exc_input < max_degin_excexc
                                    Conn(randomi(i),j) = 1;
                                    syn_strength(randomi(i),j) = syn_strength(randomi(i),j) + 1; 
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end

