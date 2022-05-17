function [Conn_dev, syn_strength_dev] = ...
    ConnDevAftercut(N,Npost,Npre,Loop_No,...
    max_deg_out,max_deg_in,Conn_Probability,Conn,syn_strength, ...
    sig_mean,radius,X,Y,Z,max_inhin,max_excin,inh_neuron_indx,exc_neuron_indx)
% This function makes all the connections between neurons based on the set
% points and probability distributions after node removal
% Input list:
% 'N' = Number of nodes post node removal
% 'Npost' = Number of neurons missed indegree
% 'Npre' = Number of neurons missed outdegree
% 'Loop_No' = number of loops to run for creating compensatory connections
% 'max_deg_out' = matrix conatins the set points of sending connections for exc. and inh. neurons
% 'max_deg_in' = matrix conatins the set points of receiving connections for exc. and inh. neurons
% 'Conn_Probability' = matrix contains probability of connections between any two nodes
% based on their arc distances
% 'Conn' = Zero NbyN connectivity matrix
% 'Syn_strength' = Zero NbyN synaptic strength matrix
% 'sig_mean' = Mean value for sigmoid function that constrains the nodes ability
%to accumulate connections too rapidly
% 'X', 'Y', 'Z' = adjusted matrix of coordinates in the lesioned network
% 'inh_neuron_indx' = index of the inhibitory neurons in the lesioned
% network
% 'exc_neuron_indx' = index of the excitatory neurons in the lesioned
% network
% 'max_inhin' = matrix of number of inhibitory inputs in the lesioned network to each neuron
% 'max_excin' = matrix of number of excitatory inputs in the lesioned network to each neuron

% Output list:
% 'Conn_dev' = Filled NbyN connectivity matrix
% 'Syn_strength_dev' = Filled NbyN synaptic strength matrix

sigmoid_func = fittype('1./(1 + exp(-a.*(x-c)))');
connections = 0 : N;
eta = 1000;
% speedW = max_deg_in-sum(syn_strength)';
Sig = 1./(1 + exp(0.005.*(connections'-sig_mean)));
Pick_Prob = fit(connections',Sig, sigmoid_func,'start',[-0.005 sig_mean]); 
max_Z = max(Z(:));
Conn_dev = Conn;
syn_strength_dev = syn_strength;
% figure;colormap 'white';
% [xs,ys,zs] = sphere; surf(radius*xs,radius*ys,radius*zs);
% hold on;plot3(X,Y,Z,'.', 'MarkerSize', 20,'color', 'b'); 
% hold on;
for hours = 1 : Loop_No
    True_out = sum(syn_strength_dev(Npost,:),2) < max_deg_out(Npost);
    Check_out = find(True_out == 1);
    True_in = sum(syn_strength_dev(:,Npre))' < max_deg_in(Npre);
    Check_in = find(True_in == 1);
    if isempty(Check_out) == 1 || isempty(Check_in) == 1
        fprintf(num2str(hours));
        break
    else
        syn_strength_dummy = syn_strength_dev;
%         syn_strength_dummy(inh_neuron_indx,:) = syn_strength_dummy(inh_neuron_indx,:) .* 1.2;
        no_axo = sum(syn_strength_dummy(Npost(Check_out),:),2);
        if length(Check_out) <= eta
            eta = length(Check_out);
        end
        if length(Check_out) == 1
            randomi = Check_out;
        else
            randomi = randsample(Npost(Check_out),eta,'true',Pick_Prob(no_axo));
            dist_prob = Conn_Probability(randomi,Npre(Check_in));
            syn_strength_dummy2 = syn_strength_dev;
%             syn_strength_dummy2(:,inh_neuron_indx) = syn_strength_dummy2(:,inh_neuron_indx) .* 1.1;
            no_den = sum(syn_strength_dummy2(:,Npre(Check_in)));
            avail_incap = Pick_Prob(no_den);
            incap = repmat(avail_incap',size(dist_prob,1),1);
            tot_prob = dist_prob .* incap;
            for i = 1:size(tot_prob,1)
                if sum(syn_strength_dev(randomi(i),:)) < max_deg_out(randomi(i))
                    if length(Check_in) == 1
                        j = Npre(Check_in);
                    else
%                         weight = tot_prob(i,:).*speedW(Npre(Check_in))';
                        weight = tot_prob(i,:);
                        j = randsample(Npre(Check_in),1,'true',weight);
                        if sum(syn_strength_dev(:,j)) < max_deg_in(j) 
                            nodeI = [X(randomi(i)),Y(randomi(i)),Z(randomi(i))];
                            nodeII = [X(j),Y(j),Z(j)];
                            mid_node = (1/2*(nodeI+nodeII) / norm(1/2*(nodeI+nodeII))).*radius;
%                             mid_node_midI = (1/2*(nodeI+mid_node) / norm(1/2*(nodeI+mid_node))).*radius;
%                             mid_node_midII = (1/2*(nodeII+mid_node) / norm(1/2*(nodeII+mid_node))).*radius;
                            if (mid_node(3) < max_Z+20)
%                                 && mid_node_midI(3) < max_Z+20 && mid_node_midII(3) < max_Z+20)
                                check_inh_out = ismember(randomi(i),inh_neuron_indx);
                                all_inh_input = sum(syn_strength_dev(inh_neuron_indx,j));
                                all_exc_input = sum(syn_strength_dev(exc_neuron_indx,j));
                                if check_inh_out == 1 % Inhibitory neuron send connection
                                    if all_inh_input < max_inhin(j)
                                        Conn_dev(randomi(i),j) = 1;
                                        syn_strength_dev(randomi(i),j) = syn_strength_dev(randomi(i),j) + 1;
%                                         L1 = [nodeI;mid_node_midI];
%                                         L2 = [mid_node_midI;mid_node];
%                                         L3 = [mid_node;mid_node_midII];
%                                         L4 = [mid_node_midII;nodeII];
%                                         plot3(L1(:,1),L1(:,2),L1(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L2(:,1),L2(:,2),L2(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L3(:,1),L3(:,2),L3(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L4(:,1),L4(:,2),L4(:,3),'.-','LineWidth',1.5);hold on;
                                    end
                                else % Excitatory neuron send connection
                                    if all_exc_input < max_excin(j) 
                                        Conn_dev(randomi(i),j) = 1;
                                        syn_strength_dev(randomi(i),j) = syn_strength_dev(randomi(i),j) + 1;
%                                         L1 = [nodeI;mid_node_midI];
%                                         L2 = [mid_node_midI;mid_node];
%                                         L3 = [mid_node;mid_node_midII];
%                                         L4 = [mid_node_midII;nodeII];
%                                         plot3(L1(:,1),L1(:,2),L1(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L2(:,1),L2(:,2),L2(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L3(:,1),L3(:,2),L3(:,3),'.-','LineWidth',1.5);hold on;
%                                         plot3(L4(:,1),L4(:,2),L4(:,3),'.-','LineWidth',1.5);hold on;
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

