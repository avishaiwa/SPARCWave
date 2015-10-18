% Run on three datasets
scat_root_dir = 'C:\Users\user\Dropbox\wavelet'; % change according to local path
figures_dir = fullfile(scat_root_dir, 'docs/figs');
datasets_dir = fullfile(scat_root_dir, 'Data');

datasets_vec = {'Berkeley', 'Phenome', 'Wheat'};
input_files_vec = {'berkeley_data', 'phoneme_data','wheat_data'}; % fill weith file names
output_files_vec = {'scatter_berkeley_data', 'scatter_phoneme_data', 'scatter_wheat_data'}; % fill weith file names
labels_vec = {{'Boys', 'Girls'}, {'sh', 'dcl', 'iy', 'aa', 'ao'}, {'High', 'Low'}};

% Parameters:
T_vec = [32, 32, 32];  % time resolution for computing scattering coefficients
Q_mat = [2 1; 8 1; 8 1];  % frequency resolution :1 / log_2(a) where coefficients are a^j
M_opt = [2, 2, 2]; % options for scattering (??)

run_flag = 0; % 0 : only plot, 1 - re-run scattering computation 

plot_opt = []; % Plotting flags: 
plot_opt.plot_flag = 1;  % plot scattering matrices using imagesc
plot_opt.plot_horizontal = 1; % 0 - plot each cluster vertically, 1 - plot each cluster horizontaly 
plot_opt.plot_bw = 0; % 0 - plot in gray-scale (black&while), 1 - plot in color 
plot_opt.normalize_colormaps = 1; % 0 - set all layers on same scale, 1 - use different maps for each scale
plot_opt.plot_weights = 1; % 1 - plot also clustering weights, 0 - plot just cluster averages 
plot_opt.weights_groups = 0; % 1 - take weights from group-sparse clustering, 0 - take weights from sparse clustering 
plot_opt.show_labels = 1; % 1 - show cluster labels, 0 - don't show 

scat_post_process_vec = {'none', 'scattergram', 'scattergram_normalized'};

for s=1:1 % length(scat_post_process_vec)
    % loop on all datasets:
    for d=1:length(datasets_vec)
        data_struct = load(fullfile(datasets_dir, input_files_vec{d})); % load data
        
        weights_file = fullfile(datasets_dir, [strrep(input_files_vec{d}, 'data', 'weights') repmat('_group', 1, plot_opt.weights_groups) '.csv'])
        plot_opt.weights = loadcellfile(weights_file, [], ',');
        if(plot_opt.weights_groups)
            plot_opt.weights = cell2num(plot_opt.weights(2:end,2:end));
        else
            plot_opt.weights = cell2num(plot_opt.weights(2:end,:));
        end
        
        if(run_flag)
            [A, S_mean] = compute_scattering_multiple_signals(data_struct.data, data_struct.label, T_vec(d), Q_mat(d,:), M_opt(d), ...
                scat_post_process_vec{s}, fullfile(datasets_dir, [output_files_vec{d} '_' scat_post_process_vec{s}])); % run scattering and save coefficient to file
        else % load A, S_mean from file
            load(fullfile(datasets_dir, [output_files_vec{d} '_' scat_post_process_vec{s}])); 
        end
        
        if(plot_opt.plot_flag)
            %    figure; imagesc(S_mean{1}) % plot scattering coefficients
            scattering_imagesc(S_mean, labels_vec{d}, plot_opt, ...
                fullfile(figures_dir, [output_files_vec{d}  '_' scat_post_process_vec{s} '_layers'])); % here do 3 image scales
        end
    end
end

