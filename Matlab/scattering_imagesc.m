% Show scattering coefficients layer by layer
% Input:
% S - cell array where S{i} is a matrix with the coefficients of the i-1-th layer
% labels_vec - label for each cluster
% plot_opt - structure with several parameters fr plotting:
%           plot_flag - how to display images:
%            0 - one image for each cluster and each layer
%            1 - one image for all clusters for each layer
% output_image_file - where to save figure (optional)
%
% Output: (none)
%
function X = scattering_imagesc( S, labels_vec, plot_opt, output_image_file)

X = [];

num_clusters = length(S);
num_layers = length(S{1});
num_nodes_in_layer = zeros(num_layers, 1);
if(~exist('labels_vec', 'var') || isempty(labels_vec)) % labels for clusters
    labels_vec = num2str_cell(num2cell((1:num_clusters)));
end
if(length(labels_vec) == num_clusters)
    labels_vec{num_clusters+1} = 'w'; % 'Weights';
end
if(plot_opt.plot_bw)
    sep_line_color = 'red';
else
    sep_line_color = 'black';
end


max_val = -99999999999.9; min_val = -max_val;
for c=1:num_clusters % Get min and max values
    max_val = max(max_val, max(max_cell(S{c})));
    min_val = min(min_val, min(min_cell(S{c})));
    
    for i=1:num_layers
        if(c == 1)
            [num_nodes_in_layer(i), num_points] = size(S{c}{i});
        end
    end
end
cum_num_nodes_in_layer = [0 cumsum(num_nodes_in_layer)'];
S_img_struct = cell(num_clusters, 1);

figure; hold on;
%cluster_sizes = zeros(num_clusters, 1);
for i=1:num_layers
    S_img = [];
    
    for c=1:num_clusters
        if((i==1) && (plot_opt.plot_horizontal))
            S_img_struct{c} = [];
        end
        plot_now=0;
        switch plot_opt.plot_flag
            case 0 % one plot for each clusterXlayer configuration
                subplot(num_plot_layers, num_clusters, c+(i-1)*num_clusters); hold on; plot_now=1;
                S_img = S{c}{i};
                if(i == num_layers)
                    xlabel(labels_vec{c});
                end
            case 1 % one plot for all clusters
                if(plot_opt.plot_horizontal == 0)
                    S_img = [S_img S{c}{i}];
                    if(c == num_clusters)
                        subplot(num_plot_layers, 1, i); hold on; plot_now=1;
                    end
                else % plot horizontally
                    S_img_struct{c} = [S_img_struct{c} S{c}{i}']; % make separate plot for each cluster
                    %%                    if(c == num_clusters)
                    %%                        S_img = reshape(S_img, num_nodes_in_layer(i)*num_clusters, num_points);
                    %%                        subplot(1, num_layers, i); hold on;  plot_now=1;
                    %%                    end
                    if(i == num_layers)
                        subplot_ax = subplot(num_clusters+plot_opt.plot_weights, 1, c); hold on;  plot_now=1;
                        S_img = S_img_struct{c}; % copy to image
                    end
                end % if plot horizontally
        end % switch plot type
        
        if(plot_now)
            if(plot_opt.plot_weights) % plot also weights
                if(plot_opt.plot_horizontal==0)
                    S_img = [S_img' plot_opt.weights((cum_num_nodes_in_layer(i)+1):cum_num_nodes_in_layer(i+1),:)']';
                end
            end
            imagesc(S_img);
            if(plot_opt.plot_bw)
                colormap(1-colormap(gray));
            end
            %            if(i == 1) % first layer
            %                ColorMap = get(gcf,'Colormap');
            %            end
            %            else % next layers
            if((c == 1) && plot_opt.plot_horizontal) % first cluster
                ColorMap = get(gcf,'Colormap');
            end
            freezeColors;
            colormap(subplot_ax, ColorMap);
            if(~plot_opt.normalize_colormaps)
                caxis manual; caxis([min_val, max_val]);
                %                    colormap(gcf, ColorMap);
            end
            %            end
            if( ((i == num_layers) || plot_opt.normalize_colormaps) && (~plot_opt.plot_horizontal) ) % first layer
                colorbar;
            end
        end % if plot_now
        
    end % loop on clusters
    
    for c=1:(num_clusters+plot_opt.plot_horizontal*plot_opt.plot_weights) % loop on clusters again
        if(c == 1)
            if(plot_opt.plot_horizontal == 0)
                ylabel(['Scat-Layer ' num2str(i)]);
            else % plot horizontal
                layer_str = []; % create a long label with all layers
                for i_layer=1:(num_layers) % plot lines - loop again
                    layer_str = [layer_str ' layer ' num2str(i_layer-1)];
                    layer_str = [layer_str repmat(' ', 1, num_nodes_in_layer(i_layer))];
                end
            end
        end % if c == 1
        
        switch plot_opt.plot_flag
            case 1 % plot lines separating clusters
                if(plot_opt.plot_horizontal == 0)
                    if(c < num_clusters)
                        line([ c*num_points+0.5,  c*num_points+0.5], ...
                            [0.5, num_nodes_in_layer(i)+0.5], 'color', sep_line_color, 'linewidth', 3);
                        ylim([0.5, num_nodes_in_layer(i)+0.5]);
                        xlim([0.5, num_clusters*num_points+0.5]);
                    end
                else % plot horiznontally
                    %%                    if(c < num_clusters+plot_opt.plot_weights)
                    %%                        line([0.5, num_points+0.5], ...
                    %%                            [c*num_nodes_in_layer(i)+0.5,  c*num_nodes_in_layer(i)+0.5], 'color', sep_line_color, 'linewidth', 3);
                    %%                        ylim([0.5, (num_clusters+plot_opt.plot_weights)*num_nodes_in_layer(i)+0.5]);
                    %5                        xlim([0.5, num_points+0.5]);
                    
                    %%       line([0.5, num_points+0.5], ...
                    %%           [c*num_nodes_in_layer(i)+0.5,  c*num_nodes_in_layer(i)+0.5], 'color', sep_line_color, 'linewidth', 3);
                    if(i == num_layers)
                        cur_ax = subplot(num_clusters+plot_opt.plot_weights, 1, c); hold on;
                        if(c == num_clusters+1) % plot_opt.plot_weights)
                            imagesc(plot_opt.weights');  colormap(cur_ax, 1-colormap(gray)); % set(gcf, 'Colormap', 1-colormap(gray)); % colormap(1-colormap(gray));
                            freezeColors;
                            
                            set(gca, 'xtick', 0.5 * ( 1+ cum_num_nodes_in_layer(1:end-1) + cum_num_nodes_in_layer(2:end) ));
                            set(gca, 'xticklabel', {'0', '1', '2'}, 'fontsize', 14, 'fontweight', 'bold');
                            xlabel('Scattering Layers', 'fontsize', 14, 'fontweight', 'bold');
                            
                            y_tick = get(gca, 'ytick');
                            if(max(y_tick) ~= round(max(y_tick)))
                                y_tick = y_tick(2:2:end-1);
                                set(gca, 'ytick', y_tick);
                            end
                            
                            
                        else % remove ticks, non-last cluster
                            set(gca, 'xtick', []);
                            set(gca,'FontSize',14, 'fontweight', 'bold')
                        end
                        if(plot_opt.show_labels)
                            ylabel(labels_vec{c});
                        end
                        xlim([0.5, sum(num_nodes_in_layer)+0.5]);
                        ylim([0.5, num_points+0.5]);
                        for i_layer=1:(num_layers-1+plot_opt.plot_weights) % plot lines - loop again
                            line([cum_num_nodes_in_layer(i_layer+1)+0.5,  cum_num_nodes_in_layer(i_layer+1)+0.5], ...
                                [0.5, num_points+0.5], 'color', sep_line_color, 'linewidth', 3);
                        end
                        y_tick = get(gca, 'ytick');
                        if(max(y_tick) ~= round(max(y_tick)))
                            y_tick = y_tick(2:2:end-1);
                            set(gca, 'ytick', y_tick);
                        end
                        
                        if(c == 1)
                            %                        title(layer_str);
                        end
                    end % last layer
                end % if not plot horitontal
        end % switch plot type
    end % loop on clusters again
end % loop on layers

for c=1:num_clusters+plot_opt.plot_weights
    subplot_ax(c) = subplot(num_clusters+plot_opt.plot_weights, 1, c);
    if(c == 1)
        title('(b)', 'fontsize', 14, 'fontweight', 'bold'); 
    end
    y_tick = get(subplot_ax(c) , 'ytick');
%    if(max(y_tick) ~= round(max(y_tick)))
    if(max(y_tick <= 3))
        y_tick = [1 2]; % y_tick(2:2:end-1);
        set(subplot_ax(c) , 'ytick', y_tick);
    end
    
    %    colormap(subplot_ax, ColorMap);
end




if(plot_opt.plot_horizontal == 1)
    colormap(ColorMap);
    axes_pos = [0.965-0.025*num_clusters, 1.055+0.0925*num_clusters, 1.03*num_clusters, -0.132+0.835*num_clusters] / (num_clusters+1)
    axes('Position', axes_pos, 'Visible', 'off');                 caxis manual; caxis([min_val, max_val]);
    colorbar('fontsize', 14); % set colorbar for all points
    cbfreeze;
    
    colormap(1-colormap(gray));
    axes_pos_weights = [0.965-0.025*num_clusters, 0.12+0.11*num_clusters, 1.03*num_clusters, 0.665-0.0125*num_clusters ] / (num_clusters+1)
    axes('Position', axes_pos_weights , 'Visible', 'off');  caxis manual; caxis([min(plot_opt.weights(:)), max(plot_opt.weights(:))]);
    colorbar('fontsize', 14); % set colorbar for weights
    cbfreeze;
end % plot horizontal
% These work for 2 cluster:
%axes_pos = [0.99, 1.437-0.1*num_clusters, 1.03*num_clusters, -0.2+0.87*num_clusters] / (num_clusters+1)
%axes('Position', axes_pos, 'Visible', 'off'); colorbar; % set colorbar for all points
%freezeColors;
%
%colormap(1-colormap(gray));
%axes_pos_weights = [0.99, 0.14+0.1*num_clusters, 1.03*num_clusters, 0.2+0.22*num_clusters ] / (num_clusters+1)



if(exist('output_image_file', 'var') && (~isempty(output_image_file))) % save image to file
    my_saveas(gcf, output_image_file, {'epsc', 'pdf', 'jpg'}); %
end
