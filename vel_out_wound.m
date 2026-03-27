% Root directory
root_dir = pwd;

% Number of folders
num_folders = 200;

% Collect data as tables (safer for alignment)
Ts = cell(num_folders,1);
valid = false(num_folders,1);

for i = 1:num_folders
    folder_name = fullfile(root_dir, num2str(i));
    data_file  = fullfile(folder_name, 'average_velocity_outside_wound.dat');

    if isfile(data_file)
        M = readmatrix(data_file);       % expects 2 cols: time, avg velocity
        if size(M,2) < 2
            fprintf('Bad format in %s (expected 2 columns), skipping...\n', data_file);
            continue
        end
        T = array2table(M, 'VariableNames', {'Time','Vel'});
        Ts{i} = T;
        valid(i) = true;
    else
        fprintf('File not found in %s, skipping...\n', folder_name);
    end
end

% If nothing found, bail out gracefully
if ~any(valid)
    error('No input files found. Nothing to average.');
end

% Outer-join all tables on Time to align data (handles unequal lengths)
Tall = Ts{find(valid,1)};
Tall.Properties.VariableNames{2} = sprintf('Vel_%d', find(valid,1));
for i = find(valid')  % iterate valid indices
    if i == find(valid,1), continue; end
    Ti = Ts{i};
    Ti.Properties.VariableNames{2} = sprintf('Vel_%d', i);
    Tall = outerjoin(Tall, Ti, 'Keys','Time', 'MergeKeys', true);
end

% Convert to numeric array, compute mean across velocity columns (omit NaN)
time_col = Tall.Time;
vel_mat  = Tall{:, setdiff(Tall.Properties.VariableNames, {'Time'})};
avg_velocity = mean(vel_mat, 2, 'omitnan');

% Save averaged data (time vs avg_velocity) safely
avg_data = [time_col, avg_velocity];
out_txt  = fullfile(root_dir, 'avg_velocity_vs_time_gamma1_240.dat');
out_csv  = fullfile(root_dir, 'avg_velocity_vs_time_gamma1_240.csv');

% Text DAT without header
writematrix(avg_data, out_txt, 'Delimiter', ' ');
% CSV with header
Tavg = table(time_col, avg_velocity, 'VariableNames', {'Time','AvgVelocity'});
writetable(Tavg, out_csv);

fprintf('Averaged data saved to:\n  %s\n  %s\n', out_txt, out_csv);

% Plot
figure('Units','normalized','Position',[0.05 0.05 0.85 0.8]);
plot(time_col, avg_velocity, 'k-', 'LineWidth', 3);
xlabel('Time (min)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Average Velocity Outside the Wound [\mum/min]', 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'FontSize', 18, 'LineWidth', 1.5);
grid off; box on;

% SAVE the figure (this was the missing piece)
png_file = fullfile(root_dir, 'avg_velocity_vs_time.png');
pdf_file = fullfile(root_dir, 'avg_velocity_vs_time.pdf');
fig_file = fullfile(root_dir, 'avg_velocity_vs_time.fig');

exportgraphics(gcf, png_file, 'Resolution', 300);  % high-quality PNG
exportgraphics(gcf, pdf_file);                     % vector PDF
savefig(gcf, fig_file);                            % MATLAB figure

fprintf('Figure saved as:\n  %s\n  %s\n  %s\n', png_file, pdf_file, fig_file);
