% Return average
%
% Input:
% Data - data matrix (signals * time. Signals are 1-dimensional)
% label - cluster labels for each signal (determine clusters)
% T - filter parameter for scattering transform: time resolution
% Q - filter parameter for scattering transform: frequency resolution. Take a^j where Q = 1/ log_2(a)
% M_opt - parameters for scattering transform
%
% Output:
% A - matrix with scattering coefficients over all signals
% S_mean - average of scattering coefficients for all signals over all clusters
%
function [A, S_mean] = compute_scattering_multiple_signals(data, label, T, Q, M_opt, scat_post_process, output_file_name)

if(~exist('output_file_name', 'var')) % if no output file name is given
    output_file_name= 'transformToR'; % default file name
end

if(ischar(data)) % allow to load data from name
    data = load(data);
    label = data.label;
    data = data.data;
end

[M, N] = size(data);
scat_opt.M = M_opt;
filt_opt = default_filter_options('audio', T); % always use 'audio'?
filt_opt.Q = Q;
filt_opt.J = T_to_J(T, filt_opt);
Wop = wavelet_factory_1d(N, filt_opt, scat_opt);
S = cell(M,1);
S_matrix = cell(M,1);

lab = unique(label);
K = length(lab); % number of clusters
S_mean = cell(K,1);
clusters_size = accumarray(label,1); % size of each cluster

for i=1:M % loop on sigmals
    if(mod(i, 10) == 0)
        run_signal_number = i
    end
    S{i} = scat(data(i,:)', Wop);
    num_vals = length(S{i}{1}.signal{1}); % number of values for each node (group) of scattering coefficients
    
    S_vec = format_scat(S{i}); % transfer to one matrix
    if(i == 1) % first time
        num_layers = length(S{i});
        num_nodes_per_layer = zeros(num_layers, 1);
        for j=1:K
            for layer = 1:num_layers % separate to layers
                num_nodes_per_layer(layer) = size(S{i}{layer}.meta.j, 2);
                S_mean{j}{layer} = zeros(num_nodes_per_layer(layer), num_vals);
            end
        end
        cum_num_nodes_per_layer = [0 cumsum(num_nodes_per_layer)'];
        num_scat_coefficients = numel(S_vec); % num_vals * sum(num_nodes_per_layer);
        A = zeros(M, num_scat_coefficients);
    end
    switch lower(scat_post_process)
        case 'none'
        case 'scattergram'
            S_vec(2:(num_nodes_per_layer(2)+1),:) = scattergram_layer(S{i}{1+1},[]);
            ctr = num_nodes_per_layer(2)+2;
            for j=1:num_nodes_per_layer(2) % number of nodes in first layer
                run_j = j
                tmp = scattergram_layer(S{i}{1+2},j-1);
                cur_num_nodes = sum(S{i}{3}.meta.j(1,:) == j-1);
                S_vec(ctr:(ctr+cur_num_nodes-1),:) = tmp(end-cur_num_nodes+1:end,:);
                ctr = ctr+cur_num_nodes;
            end
            
        case 'scattergram_normalized'
            logrenorm_S = log_scat(renorm_scat(S{i}));
            
            % Display the transformed coefficients.
            S_vec(2:(num_nodes_per_layer(2)+1),:)  = scattergram_layer(logrenorm_S{1+1},[]);
            ctr = num_nodes_per_layer(2)+2;
            for j=1:num_nodes_per_layer(2) % number of nodes in first layer
                tmp = scattergram_layer(logrenorm_S{1+2},j-1);
                cur_num_nodes = sum(S{i}{3}.meta.j(1,:) == j-1);
                S_vec(ctr:(ctr+cur_num_nodes-1),:) = tmp(end-cur_num_nodes+1:end,:);
                ctr = ctr+cur_num_nodes;
            end
            
    end % switch post-processing format
    cluster_ind = find(lab == label(i));
    for layer = 1:num_layers % separate to layers
        S_mean{cluster_ind}{layer} = S_mean{cluster_ind}{layer} + S_vec((cum_num_nodes_per_layer(layer)+1):cum_num_nodes_per_layer(layer+1),:);
    end
    A(i,:) = S_vec(:);
end

for j=1:K % Normalize
    for layer = 1:num_layers % separate to layers
        S_mean{j}{layer} =  S_mean{j}{layer} / clusters_size(j); % normalize % zeros(size(S_mean1{1}));
    end
end

save(output_file_name, 'A', 'S_mean'); % save output in .mat format.
dlmwrite(file_name_to_txt(output_file_name), A, 'delimiter', '\t', 'precision', 3); % save also in csv/txt format  - to be used for clustering




