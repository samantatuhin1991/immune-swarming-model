% === Average combined_cell_dynamics over 1..100 folders and plot ===
% Expects each file: combined_cell_dynamics.dat with at least:
%   col 1: time (min)
%   col 2: total immune cells in wound
%
% Saves:
%   avg_combined_cell_dynamics.dat  (time, mean_n_w)
%   avg_combined_cell_dynamics.csv  (Time, Mean_n_w)
%
% Plot:
%   time vs averaged total immune cells in wound (one curve), large fonts, grid off.

root_dir    = pwd;
num_folders = 100;

time_ref = [];
nw_mat   = [];

tol = 1e-10;
col = 0;

for i = 1:num_folders
    folder_name = fullfile(root_dir, num2str(i));
    data_file   = fullfile(folder_name, 'combined_cell_dynamics.dat');

    if isfile(data_file)
        data = load(data_file);  % expect [time, total_cells_in_wound, ...]
        if size(data, 2) < 2
            fprintf(2, 'Warning: %s has <2 columns. Skipping.\n', data_file);
            continue;
        end

        t  = data(:, 1);
        nw = data(:, 2);   % total immune cells in wound

        % Ensure t is strictly monotonic for interp1
        if any(diff(t) <= 0)
            [t, uniqIdx] = unique(t, 'stable');
            nw = nw(uniqIdx);
        end

        % Initialize reference time grid from the first valid file
        if isempty(time_ref)
            time_ref = t(:);
            N        = numel(time_ref);
            nw_mat   = NaN(N, num_folders);
        end

        % If time grids differ, interpolate onto time_ref
        if numel(t) ~= numel(time_ref) || any(abs(t(:) - time_ref) > tol)
            nw_i = interp1(t, nw, time_ref, 'linear', NaN);
        else
            nw_i = nw(:);
        end

        col = col + 1;
        nw_mat(:, col) = nw_i;
    else
        fprintf('File not found in %s, skipping...\n', folder_name);
    end
end

if col == 0
    error('No valid combined_cell_dynamics.dat files found in 1..%d.', num_folders);
end

% Keep only filled columns
nw_mat = nw_mat(:, 1:col);

% Average over folders (omit NaNs)
mean_nw = mean(nw_mat, 2, 'omitnan');

% ---- SAVE averaged data (DAT + CSV) ----
avg_out = [time_ref, mean_nw];
out_dat = fullfile(root_dir, 'avg_combined_cell_dynamics_alpha_35.dat');   % space-delimited
out_csv = fullfile(root_dir, 'avg_combined_cell_dynamics_alpha_35.csv');   % with header

writematrix(avg_out, out_dat, 'Delimiter', ' ');
writetable(table(time_ref, mean_nw, ...
    'VariableNames', {'Time','Mean_n_w'}), out_csv);

fprintf('Saved averaged data to:\n  %s\n  %s\n', out_dat, out_csv);

% ---- PLOT averaged total immune cells in wound ----
fig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);
hold on;

p1 = plot(time_ref, mean_nw, '-', 'LineWidth', 3, 'Color', [0 0 0]);

xlabel('Time (min)', 'FontSize', 22, 'FontWeight', 'bold');
ylabel('Total immune cells in wound', 'FontSize', 22, 'FontWeight', 'bold');

legend(p1, {'\langle n_{\mathrm{wound}}(t) \rangle'}, ...
       'FontSize', 18, 'Location', 'best');

set(gca, 'FontSize', 20, 'LineWidth', 1.8, ...
         'Box', 'on', 'XGrid', 'off', 'YGrid', 'off');

hold off;

% ---- SAVE FIGURE ----
png_file = fullfile(root_dir, 'avg_combined_cell_dynamics.png');
pdf_file = fullfile(root_dir, 'avg_combined_cell_dynamics.pdf');
fig_file = fullfile(root_dir, 'avg_combined_cell_dynamics.fig');

exportgraphics(fig, png_file, 'Resolution', 300);  % high-quality PNG
exportgraphics(fig, pdf_file);                     % vector PDF
savefig(fig, fig_file);                            % MATLAB FIG

fprintf('Figure saved as:\n  %s\n  %s\n  %s\n', png_file, pdf_file, fig_file);
