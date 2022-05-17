function [Conn_Probability, Dist, Conn_Probability_aftercut] = ...
    ArcDistCalc(N,X,Y,Z,exc_neuron_indx,inh_neuron_indx,radius,scenario)
% This function finds arc distances between any two nodes and creates
% distribution probabilities of connections. Sigmoid distributions for
% excitatory and inhibitory neurons are different as are in the realistic
% biophysical condition. Inhibitory neurons make more local connections
% whereas the excitatory neurons send longer axons.

% Input list:
% 'N' = number of neurons
% 'X', 'Y', 'Z' = cartesian coordiantes of nodes/neurons
% 'exc_neuron_indx' = index of excitatory neurons in the connectivity matrix
% 'inh_neuron_indx' = index of inhibitory neurons in the connectivity matrix
% 'radius' = radius of the sphere
% 'sigma_exc' = mean value for sigma distribution of distance wighted
% probability of connections for excitatory neurons
% 'sigma_inh' = mean value for sigma distribution of distance wighted
% probability of connections for inhibitory neurons
% 'sigma_exc_aftercut' = mean value for sigma distribution of distance wighted
% probability of connections for excitatory neurons after node removal
% 'sigma_inh_aftercut' = mean value for sigma distribution of distance wighted
% probability of connections for inhibitory neurons after node removal

% Output list:
% 'Conn_Probability' = Matrix indicating probability of connections between
% any two nodes
% 'Dist' = Arc distances between any two nodes
% 'Conn_Probability_aftercut' = Matrix indicating probability of connections between
% any two nodes after node removal

for I = 1 : N
    for II = 1 : N
        nodeI = [X(I),Y(I),Z(I)];
        nodeII = [X(II),Y(II),Z(II)];
        Dist(I,II) = radius * atan2(norm(cross(nodeI,nodeII)),dot(nodeI,nodeII));
    end
end
Dist = fix(Dist); 
Max_Dist = max(max(Dist)); % Maximum arc lenght between neurons in the defined network

sigma_exc = (1/4)* Max_Dist; % Excitatory 'mean' value in the sigmoid distribution for distance weighted probability of connection 
sigma_inh = (1/5)* Max_Dist; % inhibitory 'mean' value in the sigmoid distribution for distance weighted probability of connection 
if scenario == 1
    sigma_exc_aftercut = (1/25)* Max_Dist; % (Maximum constraint) For weighted probability of connection after cut or node removal 
    sigma_inh_aftercut = (1/25)* Max_Dist; % (Maximum constraint) For weighted probability of connection after cut or node removal 
elseif scenario == 2
    sigma_exc_aftercut = (1/12)* Max_Dist; % (Intermediate constraint) For weighted probability of connection after cut or node removal 
    sigma_inh_aftercut = (1/12)* Max_Dist; % (Intermediate constraint) For weighted probability of connection after cut or node removal
elseif scenario == 3
    sigma_exc_aftercut = (1/4)* Max_Dist; % (Same as before) For weighted probability of connection after cut or node removal 
    sigma_inh_aftercut = (1/5)* Max_Dist; % (Same as before) For weighted probability of connection after cut or node removal 
end
Dist_vector = unique(sort(Dist(:)));
Dist_vector(Dist_vector == 0) = []; % All the possible distances
Sig_Dist_exc = 0.03+0.7./(1 + exp(0.005.*(Dist_vector-sigma_exc))); % Defining a sigmoid distribution for distances for excitatory neurons
sigmoid = fittype('0.03+0.7./(1 + exp(-a.*(x-c)))');
Func_Sig_exc = fit(Dist_vector,Sig_Dist_exc, sigmoid,'start',[-5e-3 sigma_exc]); % Find the best fit 
Sig_Dist_exc_aftercut = 1./(1 + exp(0.01.*(Dist_vector-sigma_exc_aftercut))); % Defining a sigmoid distribution for distances for excitatory neurons
sigmoid_aftercut = fittype('1./(1 + exp(-a.*(x-c)))');
Func_Sig_exc_aftercut = fit(Dist_vector,Sig_Dist_exc_aftercut, sigmoid_aftercut,'start',[-1e-2 sigma_exc_aftercut]); % Find the best fit 
% title('Probability of Conection Vs Distance for Excitatory Neurons');box off; 
Conn_Probability_exc = reshape(Func_Sig_exc(Dist),[N N]); % Make the N by N matrix for probability 
Conn_Probability_exc = Conn_Probability_exc - diag(diag(Conn_Probability_exc)); % Zero out the self connections
% Next 3 lines if we want formation of self connections
% random_weight_exc = max(Sig_Dist_exc).*zeros(N,1);
% random_weight = max(Gauss_Dist).*rand(N_new,1);
% Zero_Dist_weight_exc = diag(random_weight_exc);
% Conn_Probability_exc = Conn_Probability_exc + Zero_Dist_weight_exc;
Conn_Probability_exc_aftercut = reshape(Func_Sig_exc_aftercut(Dist),[N N]); % Make the N by N matrix for probability 
Conn_Probability_exc_aftercut = Conn_Probability_exc_aftercut - diag(diag(Conn_Probability_exc_aftercut)); % Zero out the self connections
Sig_Dist_inh = 1./(1 + exp(0.007.*(Dist_vector-sigma_inh))); % Defining a sigmoid distribution for distances for inhibitory neurons
sigmoid_inh = fittype('1./(1 + exp(-a.*(x-c)))');
Func_Sig_inh = fit(Dist_vector,Sig_Dist_inh, sigmoid_inh,'start',[-7e-3 sigma_inh]); % Find the best fit of Gaussian distribution
Sig_Dist_inh_aftercut = 1./(1 + exp(0.01.*(Dist_vector-sigma_inh_aftercut))); % Defining a sigmoid distribution for distances for excitatory neurons
sigmoid_inh_aftercut = fittype('1./(1 + exp(-a.*(x-c)))');
Func_Sig_inh_aftercut = fit(Dist_vector,Sig_Dist_inh_aftercut, sigmoid_inh_aftercut,'start',[-1e-2 sigma_inh_aftercut]); % Find the best fit 
% figure('color','white'); plot(Dist_vector,Sig_Dist_inh, 'LineWidth', 3); 
% set(gca, 'fontsize', 16); set(gca,'TickDir','out','fontsize', 16, 'linewidth', 1.2);
% xlabel('Distance(um)'); ylabel('Probability');
% title('Probability of Conection Vs Distance for Inhibitory Neurons');box off; 
Conn_Probability_inh = reshape(Func_Sig_inh(Dist),[N N]); % Make the N by N matrix for probability 
Conn_Probability_inh = Conn_Probability_inh - diag(diag(Conn_Probability_inh)); % Zero out the self connections
Conn_Probability_inh_aftercut = reshape(Func_Sig_inh_aftercut(Dist),[N N]); % Make the N by N matrix for probability 
Conn_Probability_inh_aftercut = Conn_Probability_inh_aftercut - diag(diag(Conn_Probability_inh_aftercut)); % Zero out the self connections
% Next 3 lines if we want formation of self connections
% random_weight_exc = max(Gauss_Dist_inh).*zeros(N,1);
% Zero_Dist_weight_exc = diag(random_weight_exc);
% Conn_Probability_inh = Conn_Probability_inh + Zero_Dist_weight_exc;
Conn_Probability_exc(inh_neuron_indx,:) = 0;
Conn_Probability_inh(exc_neuron_indx,:) = 0;
Conn_Probability = Conn_Probability_exc + Conn_Probability_inh;
Conn_Probability_exc_aftercut(inh_neuron_indx,:) = 0;
Conn_Probability_aftercut = Conn_Probability_exc_aftercut + Conn_Probability_inh_aftercut;


%% For plotting all scenarios
sigma_exc_aftercut1 = (1/25)* Max_Dist; 
sigma_exc_aftercut2 = (1/12)* Max_Dist; % For weighted probability of connection under intermediate condition
Sig_Dist_exc_aftercut1 = 1./(1 + exp(0.01.*(Dist_vector-sigma_exc_aftercut1))); % Defining a sigmoid distribution for distances for excitatory neurons
Sig_Dist_exc_aftercut2 = 1./(1 + exp(0.01.*(Dist_vector-sigma_exc_aftercut2))); % Defining a sigmoid distribution for distances for excitatory neurons

figure('color','white'); plot(Dist_vector,Sig_Dist_exc, 'k','LineWidth', 3); 
set(gca, 'fontsize', 16); set(gca,'TickDir','out','fontsize', 24, 'linewidth', 0.75);
xlabel('Distance(um)'); ylabel('Connection Weight'); hold on;
plot(Dist_vector,Sig_Dist_inh, 'r','LineWidth', 3); 
set(gca, 'fontsize', 16); set(gca,'TickDir','out','fontsize', 24, 'linewidth', 0.75);
% xlabel('Distance(um)'); ylabel('Probability');
hold on;
plot(Dist_vector,Sig_Dist_exc_aftercut2, 'b','LineWidth', 3); 
set(gca, 'fontsize', 16); set(gca,'TickDir','out','fontsize', 24, 'linewidth', 0.75);
% xlabel('Distance(um)'); ylabel('Probability');
hold on;
plot(Dist_vector,Sig_Dist_exc_aftercut1, 'g','LineWidth', 3); 
set(gca, 'fontsize', 16); set(gca,'TickDir','out','fontsize', 24, 'linewidth', 0.75);
% xlabel('Distance(um)'); ylabel('Probability');
legend('Developing exc.','Developing inh.','Intermediate','Most confined');
box off; legend boxoff;


end