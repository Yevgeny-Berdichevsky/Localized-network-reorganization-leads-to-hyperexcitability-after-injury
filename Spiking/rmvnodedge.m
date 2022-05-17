function [N_new,Affected_prenodes,Affected_postnodes,max_deg_out_rmv,max_deg_in_rmv,...
    Conn_Probability_aftercut,Conn_rmv, syn_strength_rmv,...
    inh_neuron_indx_rmv,exc_neuron_indx_rmv,max_inhin_aftercut,max_excin_aftercut,Dist_rmv,X,Y,Z] = ...
    rmvnodedge(N, No_rmnodes,syn_strength_dev,Conn_dev,X,Y,Z,st_node,Conn_Probability_aftercut,...
    Dist,radius,max_deg_out,max_deg_in)
% This function removes the desired number of neurons set by the user at
% the inputs of main code and subsequently removes all the edges connected
% to the lost neurons as well as the connections passing through the lesion
% List of inputs:
% 'N' = number of neurons
% 'No_rmnodes' = number of nodes to remove
% 'syn_strength_dev' = matrix containing the synaptic strengths before node removal
% 'Conn_dev' = matrix containing the connectivity map before node removal
% 'X','Y','Z' = nodes coordinates before node removal
% 'st_node' = first starting node to remove (defines the location of lesion) 
% 'Conn_Probability_aftercut' = probability of connection based on distance under defined scenario (1 or 2 or 3?)
% 'Dist' = matrix of arc distances 
% 'radius' = radius of the sphere
% 'max_deg_out' = set points of sending connections  
% 'max_deg_in' = set points of receiving connections  
% List of outputs:
% 'N' = Number of remained neurons
% 'Affected_prenodes' = index of neurons lost outdegree
% 'Affected_postnodes' = index of neurons lost indegree
% 'max_deg_out_rmv' = adjusted matrix of setpoints
% 'max_deg_in_rmv' = adjusted matrix of setpoints
% 'Conn_Probability_aftercut' = adjusted matrix of probability of connections after cut
% 'Conn_rmv' = adjusted matrix of connectivity
% 'syn_strength_rmv' = adjusted matrix of synaptic strength
% 'inh_neuron_indx_rmv' = index inhibitory neurons in the lesioned network
% 'exc_neuron_indx_rmv' = index excitatory neurons in the lesioned network
% 'max_inhin_aftercut' = matrix of number of inhibitory inputs in the lesioned network to each neuron
% 'max_excin_aftercut' = matrix of number of excitatory inputs in the lesioned network to each neuron
% 'Dist_rmv' = adjusted matrix of the arc distances for the remainder of neurons

N_neuron = N; % Total number of neurons
N_new = N_neuron - No_rmnodes; % Number of remainder neurons 
syn_strength_rmv = syn_strength_dev; % New matrix for the synaptic strengths 
Conn_rmv = Conn_dev; % Connectivity map matrix
% Find nodes that lost connection directly from removed nodes
AA = Conn_dev(st_node:No_rmnodes+st_node-1,:);
[row,col] = find(AA==1);
pre_syn_lack = unique(col);
pre_syn_lack_rmv = pre_syn_lack - No_rmnodes+st_node-1;
pre_syn_lack_rmv(pre_syn_lack_rmv<=0) = []; % Neurons that lacks indegree by nodes removal
A = Conn_dev(:,st_node:No_rmnodes+st_node-1);
[row,col] = find(A==1);
post_syn_lack = unique(row);
post_syn_lack_rmv = post_syn_lack - No_rmnodes+st_node-1;
post_syn_lack_rmv(post_syn_lack_rmv<=0) = []; % Neurons that lacks outdegree by nodes removal
X(st_node:No_rmnodes+st_node-1) = []; % coordinates adjustment
Y(st_node:No_rmnodes+st_node-1) = [];
Z(st_node:No_rmnodes+st_node-1) = [];
Conn_Probability_aftercut(st_node:No_rmnodes+st_node-1,:) = [];  % adjustment for probability matrix
Conn_Probability_aftercut(:,st_node:No_rmnodes+st_node-1) = []; 
Dist_rmv = Dist;
Dist_rmv(st_node:No_rmnodes+st_node-1,:) = []; Dist_rmv(:,st_node:No_rmnodes+st_node-1) = [];
% Remove nodes and edges from connectivity map
Conn_rmv(st_node:No_rmnodes+st_node-1,:) = []; Conn_rmv(:,st_node:No_rmnodes+st_node-1) = [];
syn_strength_rmv(st_node:No_rmnodes+st_node-1,:) = []; syn_strength_rmv(:,st_node:No_rmnodes+st_node-1) = [];
% Remove edges pass through the lesion
max_Z = max(Z(:)); %Find Z of the highest neuron
max_Z_indx = find(Z == max_Z);
max_Z_indx(max_Z_indx == 0) = [];
edge_point = [X(max_Z_indx),Y(max_Z_indx),Z(max_Z_indx)]; %Coordinates of the neuron located at the highest Z, closest to the lesion area
% z_top = [0,0,1]; 
% lesion_angle = acos((dot(z_top,edge_point))/(norm(edge_point)));
Max_Dist = max(max(Dist));
%Plot the network after node removal and check for the en passant
%connections
figure;colormap 'white';
[xs,ys,zs] = sphere; surf(radius*xs,radius*ys,radius*zs);hold on;
hold on;plot3(X,Y,Z,'.', 'MarkerSize', 20,'color', 'b'); 
hold on;
colors = distinguishable_colors(500);
connections = 0;
% Find en passant connections to remove by finding 3 mid points
for i = 1:N_new
  Conn_no = find(Conn_rmv(i,:)>0);
  for ii = 1:size(Conn_no,2)
    nodeI = [X(i),Y(i),Z(i)];
    j = Conn_no(ii);
    nodeII = [X(j),Y(j),Z(j)];
    mid_node = (1/2*(nodeI+nodeII) / norm(1/2*(nodeI+nodeII))).*radius;
    mid_node_midI = (1/2*(nodeI+mid_node) / norm(1/2*(nodeI+mid_node))).*radius;
    mid_node_midII = (1/2*(nodeII+mid_node) / norm(1/2*(nodeII+mid_node))).*radius;
    if (mid_node(3) > max_Z || mid_node_midI(3) > max_Z || mid_node_midII(3) > max_Z)
        Conn_rmv(i,j) = 0;
        syn_strength_rmv(i,j) = 0;
        connections = connections + 1;
        affected_nodes1(i) = i;
        affected_nodes2(i) = j;
%         plot3(X(i),Y(i),Z(i),'.', 'MarkerSize', 20,'color', colors(connections,:)); hold on;
%         plot3(X(j),Y(j),Z(j),'.', 'MarkerSize', 20,'color', colors(connections,:)); hold on;
%         plot3(mid_node(1),mid_node(2),mid_node(3),'.', 'MarkerSize', 20,'color', 'g'); hold on;
%         plot3(mid_node_midI(1),mid_node_midI(2),mid_node_midI(3),'.', 'MarkerSize', 20,'color','y'); hold on;
%         plot3(mid_node_midII(1),mid_node_midII(2),mid_node_midII(3),'.', 'MarkerSize', 20,'color', 'y'); hold on;
%         L1 = [nodeI;mid_node_midI];
%         L2 = [mid_node_midI;mid_node];
%         L3 = [mid_node;mid_node_midII];
%         L4 = [mid_node_midII;nodeII];
%         plot3(L1(:,1),L1(:,2),L1(:,3),'.-','LineWidth',1.5,'color', colors(connections,:) );hold on;
%         plot3(L2(:,1),L2(:,2),L2(:,3),'.-','LineWidth',1.5,'color', colors(connections,:));hold on;
%         plot3(L3(:,1),L3(:,2),L3(:,3),'.-','LineWidth',1.5,'color', colors(connections,:));hold on;
%         plot3(L4(:,1),L4(:,2),L4(:,3),'.-','LineWidth',1.5,'color', colors(connections,:));hold on;
    end
  end
end

if No_rmnodes < 1000% if the number of nodes removed are less than 1000
    Affected_postnodes1 = [affected_nodes2';pre_syn_lack_rmv]; % Lost indegree
    Affected_postnodes = unique(Affected_postnodes1);
    Affected_postnodes(Affected_postnodes==0)=[];
    Affected_prenodes1 = [affected_nodes1';post_syn_lack_rmv]; % Lost outdegree
    Affected_prenodes = unique(Affected_prenodes1);
    Affected_prenodes(Affected_prenodes==0)=[];
elseif No_rmnodes >= 1000 % for 1000 and more removed nodes
    Affected_postnodes1 = [pre_syn_lack_rmv]; % Lost indegree
    Affected_postnodes = unique(Affected_postnodes1);
    Affected_postnodes(Affected_postnodes==0)=[];
    Affected_prenodes1 = [post_syn_lack_rmv]; % Lost outdegree
    Affected_prenodes = unique(Affected_prenodes1);
    Affected_prenodes(Affected_prenodes==0)=[];
end

% Remove synaptic strength of the lesion edge
% adjustment for max degree in and out vectors
inh_neuron_indx = find(max_deg_out < 300);
exc_neuron_indx = find(max_deg_out > 141);
max_deg_out_rmv = max_deg_out;
max_deg_in_rmv = max_deg_in;
max_deg_out_rmv(st_node:No_rmnodes+st_node-1) = [];
max_deg_in_rmv(st_node:No_rmnodes+st_node-1) = [];
inh_neuron_indx_rmv = find(max_deg_out_rmv < 300); %New Indices 
exc_neuron_indx_rmv = find(max_deg_out_rmv > 141);
all_inhin = sum(syn_strength_dev(inh_neuron_indx,:))'; % Sum of the number of remainder inhibitory inputs to each neuron in 2000x2000 matrix
all_inhin(st_node:No_rmnodes+st_node-1) = []; % Sum of the number of inhibitory inputs to each neuron after node removal
max_in_inhinh = all_inhin(inh_neuron_indx_rmv); % Sum of the number of inhibitory inputs to inhibitory neurons 
max_in_inhexc = all_inhin(exc_neuron_indx_rmv); % Sum of the number of inhibitory inputs to excitatory neurons
max_inhin_aftercut = zeros(N_new,1); % Sum of the number of inhibitory inputs to all neurons saved in a new matrix
max_inhin_aftercut(inh_neuron_indx_rmv,1) = max_in_inhinh; % //
max_inhin_aftercut(exc_neuron_indx_rmv,1) = max_in_inhexc; % //
all_excin = sum(syn_strength_dev(exc_neuron_indx,:))'; % Sum of the number of remainder excitatory inputs to each neuron in 2000x2000 matrix
all_excin(st_node:No_rmnodes+st_node-1) = []; % Sum of the number of excitatory inputs to each neuron after node removal
max_in_excinh = all_excin(inh_neuron_indx_rmv); % Sum of the number of excitatory inputs to inhibitory neurons after node removal
max_in_excexc = all_excin(exc_neuron_indx_rmv); % Sum of the number of inhibitory inputs to excitatory neurons after node removal
max_excin_aftercut = zeros(N_new,1); % Sum of the number of excitatory inputs to all neurons saved in a new matrix
max_excin_aftercut(inh_neuron_indx_rmv,1) = max_in_excinh; % //
max_excin_aftercut(exc_neuron_indx_rmv,1) = max_in_excexc; % //
max_deg_out_pre = sum(syn_strength_dev,2); % Sum of the number of sending connections from each neuron before 
max_deg_in_pre = sum(syn_strength_dev);  % Sum of the number of sending connections from each neuron
max_deg_out_rmv = max_deg_out_pre;
max_deg_out_rmv(st_node:No_rmnodes+st_node-1) = [];
max_deg_in_rmv = max_deg_in_pre';
max_deg_in_rmv(st_node:No_rmnodes+st_node-1) = [];
end