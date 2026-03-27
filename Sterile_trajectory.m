%Author: T.S.
%  100 minutes, 
%   Length: micrometers [µm]; Time: minutes [min]; Velocity: [µm/min]; Counts: cells
%

clear; clc; close all;
%rng(42);  % reproducibility
rng('shuffle'); 

%% ---------------------- ORIGINAL NUMERIC PARAMETERS (as given) -----------
R_model      = 5;        % [µm]
rw_model     = 20;       % [µm]
rmax_model   = 150;      % [µm]
c0_model     = 0.1;      % chemo baseline (arb. units)
r0_model     = 35;       % [µm]
cs_model     = 1;        % Hill half-saturation
v0_model     = 0.75;     % speed scale (per-step; will convert)
alpha_model  = 35;       % secretion scaling (kept exactly)
gamma0_model = 0.06;     % stopping rate (per-step; will convert)
gamma1_model = 12;       % modulation (per-step; will convert)
num_cells    = 120;
t_steps      = 5000;     % desired number of steps
num_bins     = 20;
eps_val      = 100;      % stabilizer in gamma() denominator
reservoir_ratio = 0.95;  % keep ≥95% of initial outermost-bin count

%% ---------------------- TIME: 10000 steps = 100 minutes ------------------
% dt_phys = 150 / t_steps;          % [min/step] = 0.01
% time    = (0:t_steps) * dt_phys;  % [min], length = 10001
% T       = numel(time);

dt_phys   = 0.05;                        % [min/step]
total_time = 100;                        % [min]
t_steps    = round(total_time / dt_phys);% ≈ 2000 steps
time       = (0:t_steps) * dt_phys;      % [min], goes 0 → 100
T          = numel(time);

%% ---------------------- Unit-bridge to physical (per-minute) ------------
% Preserve per-step effects: v*dt and gamma*dt
v0     =  v0_model     / dt_phys;   % [µm/min]
gamma0 =  gamma0_model / dt_phys;   % [1/min]
gamma1 =  gamma1_model / dt_phys;   % [1/min]

% Geometry / chemo (already in µm, same numeric values)
R    = R_model;
rw   = rw_model;
rmax = rmax_model;
c0   = c0_model;
r0   = r0_model;
cs   = cs_model;
alpha= alpha_model;

%% ---------------------- Initialization -----------------------------------
% Random positions in [rw, rmax], uniform in area
theta   = 2*pi*rand(1, num_cells);
r_cell  = sqrt(rand(1, num_cells) * (rmax^2 - (rw+0.001)^2) + (rw+0.001)^2);
x_cells = r_cell .* cos(theta);
y_cells = r_cell .* sin(theta);
is_in_wound = false(1, num_cells);

% Trajectory tracking buffers and arrival times
trajX = cell(1, num_cells);    % x(t) history per cell
trajY = cell(1, num_cells);    % y(t) history per cell
trajT = cell(1, num_cells);    % time history per cell [min]
arrival_time = inf(1, num_cells);  % first time (min) inside wound

% Binning
r_bin_edges = linspace(rw, rmax, num_bins + 1);
bin_areas   = pi * (r_bin_edges(2:end).^2 - r_bin_edges(1:end-1).^2);
bin_centers = 0.5*(r_bin_edges(1:end-1) + r_bin_edges(2:end));

% Initial outermost-bin count (correct definition)
r_init = sqrt(x_cells.^2 + y_cells.^2);
outer_mask_init    = (r_init >= r_bin_edges(end-1) & r_init < r_bin_edges(end));
outer_init_count   = sum(outer_mask_init);
outer_init_area    = bin_areas(end);
outer_init_density = outer_init_count / outer_init_area;          %#ok<NASGU>
target_outer_count = ceil(reservoir_ratio * outer_init_count);    % target count

%% ---------------------- Chemo & velocity (STABILIZED) --------------------
% Raw chemo (can be negative / huge if used directly)
chemo_raw = @(r, c0v, ns) (c0v + alpha * ns) .* exp(-r / r0);

% Clamp chemo to a physical, non-negative range to avoid blowups:
c_max = 1e4;  % reasonable upper bound; adjust if needed
chemo_concentration = @(r, c0v, ns) ...
    min(max(chemo_raw(r, c0v, ns), 0), c_max);

% Hill function gets only non-negative chemo values (0 ≤ c ≤ c_max)
hill_function       = @(c) c ./ (c + cs);

% Velocity profile; with the clamped chemo, |hill_function(...)| ≤ 1
velocity_profile    = @(r, c0v, ns) v0 * ( ...
                            hill_function(chemo_concentration(r + R, c0v, ns)) - ...
                            hill_function(chemo_concentration(max(r - R, 0), c0v, ns)) );

%% ---------------------- Storage -----------------------------------------
% population size changes (due to reservoir) → use cell arrays for per-time velocities
velocities_list        = cell(T, 1);    % 1 x N_t signed radial [µm/min]
va_rw_list             = cell(T, 1);    % 1 x N_t [µm/min] at r=rw
cell_density           = zeros(num_bins, T);   % [cells/µm^2]
mean_velocity          = zeros(num_bins, T);   % [µm/min]
overall_mean_velocity  = zeros(1, T);         % [µm/min]
rho_rw                 = zeros(1, T);         % [cells/µm^2]
gamma_vec              = zeros(1, T);         % [1/min]
n_s                    = zeros(1, T);         % [cells]
n_ss                   = zeros(1, T);         % [cells]
n_w                    = zeros(1, T);         % [cells]
flux                   = zeros(1, T);         % [cells/min]
num_cells_in_wound     = zeros(1, T);         % [cells]
density_outside_wound  = zeros(T, 1);         % [cells/µm^2]

% Data log for cells inside wound
fid = fopen('cells_in_wound.data', 'w');

%% ---------------------- Video setup (robust + time overlay) --------------
make_video   = true;     % set false to skip writing video
video_stride = 100;      % capture every 100 steps
if make_video
    v = VideoWriter('cell_dynamics.avi');
    v.FrameRate = 30;
    open(v);
end
fig = figure('Position',[100,100,1280,960],'Resize','off','Color','w');
ax  = axes('Parent',fig,'Position',[0,0,1,1]); hold(ax,'on');
axis(ax,'equal'); xlim(ax,[-rmax rmax]); ylim(ax,[-rmax rmax]);

%% ===================== High-quality snapshot setup =======================
% 50 evenly spaced times from 0→100 min, saved at 600 DPI (PNG) and optional PDF.
snap_count    = 50;
snap_times    = linspace(0, total_time, snap_count);   % [min]
snap_dpi      = 600;                                   % 300–1200 suggested
snap_save_pdf = true;                                  % also save vector PDF
snap_dir      = sprintf('snaps_%ddpi', snap_dpi);
if ~exist(snap_dir,'dir'), mkdir(snap_dir); end

% Dedicated clean figure (no legend/southoutside text) for publication snaps
snapFig = figure('Color','w','Units','inches','Position',[1,1,7,7], 'Visible','off');
snapAx  = axes('Parent',snapFig,'Position',[0.12 0.12 0.76 0.76]); hold(snapAx,'on');
axis(snapAx,'equal'); xlim(snapAx,[-rmax rmax]); ylim(snapAx,[-rmax rmax]);
xlabel(snapAx,'X [\mum]'); ylabel(snapAx,'Y [\mum]');
box(snapAx,'on'); grid(snapAx,'off');

% Precomputed circle
th_circ = linspace(0, 2*pi, 360);

% Next snapshot index to hit
nextSnapIdx = 1;

%% ---------------------- Main loop ----------------------------------------
for t = 1:T
    % State and distances
    num_cells_in_wound(t) = sum(is_in_wound);
    if mod(t, 500) == 1
        fprintf('Time: %.2f minutes, Cells in Wound: %d\n', time(t), num_cells_in_wound(t));
    end
    r_now = sqrt(x_cells.^2 + y_cells.^2);

    % Chemo coupling uses previous n_s
    n_s_prev = n_s(max(t-1,1));

    % --- per-time containers (match current N) ---
    N       = numel(x_cells);
    v_t     = zeros(1, N);             % signed radial speeds at time t
    varw_t  = zeros(1, N);             % speeds evaluated at r=rw (same for all cells)

    % --- Update cells (explicit Euler on radial direction, STABILIZED) ---
    for i = 1:N
        r_i = r_now(i);

        % raw velocity from chemo gradient
        v_i = velocity_profile(r_i, c0, n_s_prev);  % [µm/min]

        % defensive checks against NaN/Inf and runaway magnitudes
        if ~isfinite(v_i)
            v_i = 0;
        end
        % enforce |v_i| ≤ v0 to be safe
        v_i = max(min(v_i, v0), -v0);

        v_t(i)    = v_i;
        varw_t(i) = velocity_profile(rw, c0, n_s_prev);

        if r_i > 0
            ux = x_cells(i)/r_i; 
            uy = y_cells(i)/r_i;
            x_cells(i) = x_cells(i) + v_i * ux * dt_phys;
            y_cells(i) = y_cells(i) + v_i * uy * dt_phys;

            % Hard projection into the domain if a cell steps outside rmax
            r_new = sqrt(x_cells(i)^2 + y_cells(i)^2);
            if r_new > rmax
                scale = rmax / r_new;
                x_cells(i) = x_cells(i) * scale;
                y_cells(i) = y_cells(i) * scale;
            end
        end

        if r_i <= rw && ~is_in_wound(i)
            is_in_wound(i) = true;
        end
    end

    % --- Save cell data for those in wound (AFTER motion this step) ---
    r_now_after = sqrt(x_cells.^2 + y_cells.^2);
    for i = 1:N
        if is_in_wound(i)
            r_i = max(r_now_after(i), 1e-12);
            vx  = v_t(i) * x_cells(i) / r_i;
            vy  = v_t(i) * y_cells(i) / r_i;
            fprintf(fid, 'Time: %.2f, X: %.2f, Y: %.2f, Vx: %.3f, Vy: %.3f\n', ...
                time(t), x_cells(i), y_cells(i), vx, vy);
        end
    end

    % ================= REPLENISH OUTERMOST BIN (reservoir) =================
    % Keep outermost-bin count ≥ target_outer_count
    r_now_after = sqrt(x_cells.^2 + y_cells.^2);  % recompute after motion
    outer_mask  = (r_now_after >= r_bin_edges(end-1) & r_now_after < r_bin_edges(end));
    outer_count = sum(outer_mask);
    to_add = max(0, target_outer_count - outer_count);

    if to_add > 0
        % Sample new cells uniformly in area within the outermost bin
        r_lo = r_bin_edges(end-1); 
        r_hi = r_bin_edges(end);
        U     = rand(1, to_add);
        r_new = sqrt(r_lo^2 + (r_hi^2 - r_lo^2).*U);
        th    = 2*pi*rand(1, to_add);
        x_new = r_new .* cos(th);
        y_new = r_new .* sin(th);

        % Append to state
        x_cells     = [x_cells, x_new];
        y_cells     = [y_cells, y_new];
        is_in_wound = [is_in_wound, false(1, to_add)];

        % Compute and append velocities for these new cells at the SAME time step
        v_new   = arrayfun(@(r) velocity_profile(r, c0, n_s_prev), r_new);
        % defensive caps for new cells as well
        v_new(~isfinite(v_new)) = 0;
        v_new = max(min(v_new, v0), -v0);

        v_t     = [v_t, v_new]; %#ok<AGROW>
        varw_t  = [varw_t, repmat(velocity_profile(rw, c0, n_s_prev), 1, to_add)]; %#ok<AGROW>

        % Update r_now_after after adding new cells
        r_now_after = sqrt(x_cells.^2 + y_cells.^2);
    end
    % ======================================================================

    % Expand trajectory arrays if new reservoir cells were added
    currentN   = numel(x_cells);
    oldN_arr   = numel(arrival_time);
    if oldN_arr < currentN
        arrival_time(oldN_arr+1:currentN) = inf;
    end
    if numel(trajX) < currentN
        trajX{currentN} = [];
        trajY{currentN} = [];
        trajT{currentN} = [];
    end

    % Append this step's position & time for all cells
    for i = 1:currentN
        if isempty(trajX{i})
            % first time this cell appears
            trajX{i} = x_cells(i);
            trajY{i} = y_cells(i);
            trajT{i} = time(t);
        else
            trajX{i}(end+1) = x_cells(i);
            trajY{i}(end+1) = y_cells(i);
            trajT{i}(end+1) = time(t);
        end
    end

    % Record first arrival times
    for i = 1:currentN
        if is_in_wound(i) && isinf(arrival_time(i))
            arrival_time(i) = time(t);  % minutes
        end
    end

    % --- Near-wound density (use AFTER replenishment) ---
    delta_r   = 2*R;
    ring_mask = (r_now_after >= rw & r_now_after < rw + delta_r);
    area_near = pi * ((rw + delta_r)^2 - rw^2);
    rho_rw(t) = sum(ring_mask) / area_near;

    % --- Bin densities & mean |v| (AFTER replenishment) ---
    total_velocity_sum   = 0;
    total_particle_count = 0;
    for b = 1:num_bins
        bin_indices = (r_now_after >= r_bin_edges(b) & r_now_after < r_bin_edges(b + 1));
        cell_count  = sum(bin_indices);
        cell_density(b, t) = cell_count / bin_areas(b);
        if cell_count > 0
            v_mean = mean(abs(v_t(bin_indices)));
            mean_velocity(b, t) = v_mean;
            total_velocity_sum  = total_velocity_sum + sum(abs(v_t(bin_indices)));
            total_particle_count = total_particle_count + cell_count;
        end
    end
    overall_mean_velocity(t) = total_velocity_sum / max(1, total_particle_count);

    % --- Density just outside wound (first bin) ---
    density_outside_wound(t) = sum(r_now_after >= r_bin_edges(1) & r_now_after < r_bin_edges(2)) / bin_areas(1);
    if t == 1
        density_at_bin_time_zero = density_outside_wound(1); %#ok<NASGU>
    end

    % --- Flux into wound (cells/min) ---
    va_rw_mean = mean(abs(varw_t));
    flux(t)    = 2 * pi * rw * rho_rw(t) * va_rw_mean;

    % --- Phenotype dynamics inside wound (SEMI-IMPLICIT, STABLE) ---
    if t == 1
        gamma_vec(t) = gamma0;
    else
        gamma_vec(t) = gamma0 + gamma1 * ( n_ss(t-1) / ( n_ss(t-1) + eps_val ) );
    end

    gamma_t  = gamma_vec(t);
    ns_prev  = n_s_prev;
    nss_prev = n_ss(max(t-1,1));

    % Semi-implicit update for n_s:
    %   n_s^{t+1} = n_s^t + dt*(flux - gamma * n_s^{t+1})
    % => (1 + dt*gamma)*n_s^{t+1} = n_s^t + dt*flux
    n_s_new = (ns_prev + dt_phys * flux(t)) / (1 + dt_phys * gamma_t);
    n_s_new = max(n_s_new, 0);   % enforce non-negativity

    % Update n_ss with the new n_s
    n_ss_new = nss_prev + dt_phys * gamma_t * n_s_new;
    n_ss_new = max(n_ss_new, 0);

    n_s(t)  = n_s_new;
    n_ss(t) = n_ss_new;
    n_w(t)  = n_s_new + n_ss_new;

    % --- Store per-time velocity vectors (variable length) ---
    velocities_list{t} = v_t;
    va_rw_list{t}      = varw_t;

    % -------------------- Video (clear colors + legend + time overlay) ----
    if make_video && mod(t, video_stride)==0
        cla(ax);

        % circles
        th = linspace(0, 2*pi, 360);
        % outer boundary (green dashed)
        h_outer = plot(ax, rmax*cos(th), rmax*sin(th), 'g--', 'LineWidth', 1.2); hold(ax,'on');
        % bin edges (light gray, excluding rw and rmax)
        for rr = r_bin_edges(2:end-1)
            plot(ax, rr*cos(th), rr*sin(th), 'Color', [0.7 0.7 0.7], 'LineStyle','--', 'LineWidth', 0.8);
        end
        % wound boundary (black solid)
        h_wound = plot(ax, rw*cos(th), rw*sin(th), 'k-', 'LineWidth', 1.8);

        % cells
        h_in  = plot(ax, x_cells(is_in_wound),  y_cells(is_in_wound),  'ro', 'MarkerSize', 5, 'MarkerFaceColor','r');
        h_out = plot(ax, x_cells(~is_in_wound), y_cells(~is_in_wound), 'bo', 'MarkerSize', 3);

        axis(ax,'equal'); xlim(ax,[-rmax rmax]); ylim(ax,[-rmax rmax]);
        xlabel(ax,'X [\mum]'); ylabel(ax,'Y [\mum]');

        legend(ax, [h_wound, h_outer, h_in, h_out], ...
               {'Wound boundary (r = r_w)', 'Outer boundary (r = r_{max})', ...
                'Cells in wound', 'Cells outside wound'}, ...
               'Location','southoutside','NumColumns',2,'Box','off');

        % Time overlay
        txt_time = sprintf('Time = %.2f min', time(t));
        text(ax, -0.95*rmax, 0.92*rmax, txt_time, ...
            'FontSize',12, 'FontWeight','bold', ...
            'BackgroundColor','w', 'Margin',4);

        drawnow;
        frame = getframe(fig);
        if ~isempty(frame.cdata)
            writeVideo(v, frame);
        end
    end

    % ================== Save high-quality snapshots =======================
    while nextSnapIdx <= snap_count && time(t) + 1e-9 >= snap_times(nextSnapIdx)
        % Split sets
        maskW  = is_in_wound;
        maskNW = ~is_in_wound;

        % Render a clean snapshot (no legend/writeups)
        render_clean_snapshot(snapAx, rmax, rw, th_circ, ...
                              x_cells(maskW), y_cells(maskW), ...
                              x_cells(maskNW), y_cells(maskNW));

        baseName = sprintf('snap_%03d_%06.2fmin', nextSnapIdx, snap_times(nextSnapIdx));
        pngPath  = fullfile(snap_dir, [baseName '.png']);
        try
            exportgraphics(snapAx, pngPath, 'Resolution', snap_dpi, 'BackgroundColor','white', 'ContentType','image');
        catch
            set(snapFig,'PaperPositionMode','auto');
            print(snapFig, pngPath, '-dpng', ['-r' num2str(snap_dpi)]);
        end
        if snap_save_pdf
            pdfPath = fullfile(snap_dir, [baseName '.pdf']);
            try
                exportgraphics(snapAx, pdfPath, 'ContentType','vector', 'BackgroundColor','white');
            catch
                print(snapFig, pdfPath, '-dpdf', '-bestfit');
            end
        end
        fprintf('Saved snapshot %d/%d at t=%.2f min -> %s\n', nextSnapIdx, snap_count, snap_times(nextSnapIdx), pngPath);
        nextSnapIdx = nextSnapIdx + 1;
    end
    % ======================================================================
end

% Close resources
fclose(fid);
if make_video, close(v); end

%% ---------------------- Post-processing & plots (ALL FIGURES) ------------
% Time-averaged mean |v| per bin
average_velocity_bins = zeros(1, num_bins);
for b = 1:num_bins
    vv = mean_velocity(b, :);
    average_velocity_bins(b) = mean(vv(~isnan(vv)));
end

% ---- Initial density bar plot ----
figure;
bar(bin_centers, cell_density(:, 1));
xlabel('Distance from Wound Center r [\mum]');
ylabel('Cell Density at t = 0 [cells/\mum^2]');
title('Initial Cell Density Across Annular Rings'); grid on;

% ---- Heatmap: Mean |velocity| [µm/min] ----
figure('Units','inches','Position',[1,1,8,6]);
imagesc(time, bin_centers, mean_velocity); colorbar;
xlabel('Time [min]'); ylabel('Distance from Wound Center r [\mum]');
title('Mean |Velocity| of Cells in Annular Rings [\mum/min]');
set(gca,'YDir','normal','FontSize',12,'LineWidth',1.5);
[X, Y] = meshgrid(time, bin_centers);
data_mv = [X(:), Y(:), mean_velocity(:)];
save('mean_velocity_heatmap.dat', 'data_mv', '-ascii', '-double');
print(gcf,'mean_velocity_heatmap.png','-dpng','-r300');

% ---- Mean velocity vs distance over time (each column is a time) ----
figure;
plot(r_bin_edges(1:end-1), mean_velocity, 'v-', 'LineWidth', 1);
xlabel('Distance from Wound Center r [\mum]');
ylabel('Mean |Velocity| [\mum/min]');
title('Mean Velocity per Bin over Time'); grid on;
data_mean_velocity = [r_bin_edges(1:end-1)', mean_velocity];
save('mean_velocity_vs_distance.dat','data_mean_velocity','-ascii');
saveas(gcf,'mean_velocity_distance.png');

% ---- Time-averaged |v| vs distance ----
figure;
plot(r_bin_edges(1:end-1), average_velocity_bins, 'o-', 'LineWidth', 2);
xlabel('Distance from Wound Center r [\mum]');
ylabel('Time-avg Mean |Velocity| [\mum/min]');
title('Time-averaged Mean |Velocity| vs Distance'); grid on;
data_averaged_velocity = [r_bin_edges(1:end-1)', average_velocity_bins'];
save('averaged_velocity_vs_distance.dat','data_averaged_velocity','-ascii');
saveas(gcf,'averaged_velocity_distance.png');

% ---- Overall mean |v| vs time (all cells) ----
figure;
plot(time, overall_mean_velocity, 'v-', 'LineWidth', 2);
xlabel('Time [min]'); ylabel('Overall Mean |Velocity| [\mum/min]');
title('Overall Mean |Velocity| of Cells'); grid on;
data_overall_mean_v = [time(:), overall_mean_velocity(:)];
save('average_velocity_outside_wound.dat','data_overall_mean_v','-ascii');
saveas(gcf,'average_velocity_outside_wound_all_cells.png');

% ---- Cell density over time per ring ----
figure; hold on;
for b = 1:num_bins
    plot(time, cell_density(b,:), 'DisplayName', sprintf('Bin %d', b));
end
xlabel('Time [min]'); ylabel('Cell Density [cells/\mum^2]');
title('Cell Density in Annular Rings Over Time'); legend; grid on; colormap(jet(num_bins));
data_cell_density = [time; cell_density]';
save('cell_density_over_time.dat','data_cell_density','-ascii');
saveas(gcf,'Cells_density_time.png');

% ---- Mean velocity over time per ring ----
figure; hold on;
for b = 1:num_bins
    plot(time, mean_velocity(b,:), 'DisplayName', sprintf('Bin %d', b));
end
xlabel('Time [min]'); ylabel('Mean |Velocity| [\mum/min]');
title('Mean |Velocity| in Annular Rings Over Time'); legend; grid on; colormap(jet(num_bins));
data_mean_velocity_bins = [time; mean_velocity]';
save('mean_velocity_in_bins_over_time.dat','data_mean_velocity_bins','-ascii');
saveas(gcf,'Mean_velocity_time.png');

% ---- Cell density heatmap ----
figure;
imagesc(time, bin_centers, cell_density); colorbar;
xlabel('Time [min]'); ylabel('Distance from Wound Center r [\mum]');
title('Cell Density in Annular Rings [cells/\mum^2]');
set(gca,'YDir','normal');
[TT, RR] = meshgrid(time, bin_centers);
data_cd = [TT(:), RR(:), cell_density(:)];
save('cell_density_heatmap.dat','data_cd','-ascii');
saveas(gcf,'Cell_denisty_heatmap.png');

% ---- Density just outside wound (first ring) ----
figure;
plot(time, density_outside_wound, 'k-', 'LineWidth', 2); hold on;
yline(density_outside_wound(1), '--', 'Density at t=0', 'LabelHorizontalAlignment','left');
xlabel('Time [min]'); ylabel('\rho [cells/\mum^2]');
title('Density in Bin Just Outside the Wound'); grid on;
data_density_outside = [time(:), density_outside_wound(:)];
save('density_outside_wound.dat','data_density_outside','-ascii');
saveas(gcf,'Cell_denisty_justoutside_wound_time.png');

% ---- Influx diagnostics ----
cumulative_influx  = cumtrapz(time, flux);
total_cells_influx = trapz(time, flux); %#ok<NASGU>

% Combined plot (cells in wound + cumulative influx)
figure;
plot(time, num_cells_in_wound, 'r-', 'LineWidth', 2, 'DisplayName', 'Cells in Wound'); hold on;
plot(time, cumulative_influx,  'b--', 'LineWidth', 1.5, 'DisplayName', 'Cumulative Influx');
xlabel('Time [min]'); ylabel('Value');
title('Cell Dynamics and Influx in Wound Healing'); legend('show'); grid on;
saveas(gcf,'Combined_Cell_Dynamics_Plot.png');

% Save combined data
data_combined = [time(:), num_cells_in_wound(:), cumulative_influx(:)];
save('combined_cell_dynamics.dat','data_combined','-ascii');

% ---- Log-log plot (cells in wound & cumulative influx) ----
figure;
loglog(time, max(num_cells_in_wound,1e-12), 'r-', 'LineWidth', 2, 'DisplayName', 'Cells in Wound'); hold on;
loglog(time, max(cumulative_influx,1e-12),  'b--', 'LineWidth', 1.5, 'DisplayName', 'Cumulative Influx');
xlabel('Time [min]'); ylabel('Value');
title('Log-Log Plot of Cell Dynamics and Influx'); legend('show'); grid on;
data_combined_ll = [time(:), num_cells_in_wound(:), cumulative_influx(:)];
save('loglog_combined_cell_dynamics.dat','data_combined_ll','-ascii');
saveas(gcf,'LogLog_Cell_Dynamics_Plot.png');

% ---- Overall absolute mean velocity vs time (from velocities_list) ----
average_velocity = zeros(1, T);
for t = 1:T
    vt = velocities_list{t};
    if ~isempty(vt)
        average_velocity(t) = mean(abs(vt));
    end
end
figure;
plot(time, average_velocity, 'r-', 'LineWidth', 2, 'DisplayName', 'Absolute Mean Velocity');
xlabel('Time [min]'); ylabel('Absolute Mean Velocity [\mum/min]');
title('Overall Absolute Mean Velocity vs Time'); grid on; legend('show');
saveas(gcf,'Absolute_Mean_Velocity_Plot_Linear.png');
save('absolute_mean_velocity_linear.dat','-ascii','average_velocity'); % data values alone
dlmwrite('absolute_mean_velocity_linear.dat',[time(:), average_velocity(:)],'delimiter',' '); % (time,value)

% ---- Log-scale average velocity ----
figure;
loglog(time, max(average_velocity,1e-12), 'b-', 'LineWidth', 2, 'DisplayName', 'Average Mean Velocity');
xlabel('Time [min]'); ylabel('Average Mean Velocity [\mum/min]');
title('Log Scale Overall Average Mean Velocity vs Time'); grid on; legend('show');
saveas(gcf,'Overall_average_mean_velocity_plot.png');
dlmwrite('overall_absolute_mean_velocity_log_scale.dat',[time(:), average_velocity(:)],'delimiter',' ');

% ---- va_rw statistics (from va_rw_list) ----
mean_va_rw = zeros(1, T);
for t = 1:T
    vt = va_rw_list{t};
    if ~isempty(vt)
        mean_va_rw(t) = mean(abs(vt));
    end
end
figure;
plot(time, mean_va_rw, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean Abs(va\_rw)');
xlabel('Time [min]'); ylabel('Mean Abs(va\_rw) [\mum/min]');
title('Mean Absolute Value at Wound Ring vs Time'); grid on; legend('show');
saveas(gcf,'Mean_Abs_va_rw_Plot_Linear.png');
dlmwrite('mean_abs_va_rw_linear.dat',[time(:), mean_va_rw(:)],'delimiter',' ');

% ---- Log-log va_rw ----
figure;
loglog(time, max(mean_va_rw,1e-12), 'r-', 'LineWidth', 2, 'DisplayName', 'Mean Abs(va\_rw)');
xlabel('Time [min]'); ylabel('Mean Abs(va\_rw) [\mum/min]');
title('Log-Log: Mean Abs(va\_rw) vs Time'); grid on; legend('show');
saveas(gcf,'Mean_Abs_va_rw_Plot_LogLog.png');
dlmwrite('mean_abs_va_rw_loglog.dat',[time(:), mean_va_rw(:)],'delimiter',' ');

% ---- Gamma and populations ----
figure;
plot(time, gamma_vec, 'b', 'LineWidth', 2);
xlabel('Time [min]'); ylabel('\gamma(t) [1/min]');
title('\gamma vs Time'); grid on;
saveas(gcf,'gamma_vs_time_final.png');
dlmwrite('time_vs_gamma_final.dat',[time(:), gamma_vec(:)],'delimiter',' ');

figure;
plot(time, n_s, 'b', 'LineWidth', 2); hold on;
plot(time, n_ss,'r', 'LineWidth', 2);
plot(time, n_w, 'g', 'LineWidth', 2); hold off;
legend('n_s (Secreting Cells)','n_{ss} (Stopped Cells)','n_w (Total Cells)','Location','best');
xlabel('Time [min]'); ylabel('Cell Count [cells]');
title('Immune Cell Dynamics in the Wound'); grid on;
saveas(gcf,'immune_cell_dynamics_final.png');
dlmwrite('immune_cell_dynamics_final.dat',[time(:), n_s(:), n_ss(:), n_w(:)],'delimiter',' ');

%% ===================== SAVE FULL TRAJECTORIES (ALL CELLS) =====================
% File format:
%   cell_id  step_index  time_min  x_um  y_um  arrival_time_min
% For cells that never reach the wound, arrival_time_min = -1.

all_traj_file = 'all_cell_trajectories.dat';
fid_all = fopen(all_traj_file, 'w');
fprintf(fid_all, '# cell_id step_index time_min x_um y_um arrival_time_min\n');

num_cells_final = numel(trajX);
for i = 1:num_cells_final
    x_i = trajX{i};
    y_i = trajY{i};
    t_i = trajT{i};
    L   = numel(x_i);
    if isempty(t_i)
        continue;
    end
    if isinf(arrival_time(i))
        a_i = -1;   % never arrived
    else
        a_i = arrival_time(i);
    end
    for k = 1:L
        fprintf(fid_all, '%d %d %.6f %.6f %.6f %.6f\n', ...
            i, k, t_i(k), x_i(k), y_i(k), a_i);
    end
end

fclose(fid_all);
fprintf('Saved ALL cell trajectories to %s\n', all_traj_file);

%% ===================== SAVE ARRIVER TRAJECTORIES (SUBSET) =====================
Tcut = 100;  % [min] fixed time horizon (full simulation time)
arrived_idx = find(arrival_time <= Tcut);

traj_filename = 'arriver_trajectories.dat';
fid_traj = fopen(traj_filename,'w');
fprintf(fid_traj, '# cell_id step_index time_min x_um y_um arrival_time_min\n');

for kk = 1:numel(arrived_idx)
    i  = arrived_idx(kk);
    ti = trajT{i};
    xi = trajX{i};
    yi = trajY{i};
    L  = numel(ti);
    for k = 1:L
        if ti(k) > Tcut + 1e-9
            break;  % only keep up to Tcut
        end
        fprintf(fid_traj, '%d %d %.6f %.6f %.6f %.6f\n', ...
            i, k, ti(k), xi(k), yi(k), arrival_time(i));
    end
end

fclose(fid_traj);
fprintf('Saved arriver trajectories to %s\n', traj_filename);


% Helper to map arrival time to a colormap index (earlier=dark, later=bright)
cmap = parula(256);                          % perceptually uniform
if isempty(arrived_idx)
    tmin = 0; tmax = 1;
else
    tmin = min(arrival_time(arrived_idx));
    tmax = max(arrival_time(arrived_idx));
    if ~isfinite(tmin), tmin = 0; end
    if ~isfinite(tmax), tmax = tmin + 1; end
end
mapColor = @(ta) cmap( max(1, min(256, ...
    round(1 + (ta - tmin) * (256-1) / max(1e-12, (tmax - tmin))) )), :);

% (A) 2D Trajectory overlay (0..Tcut) for all arrivers
figure('Color','w','Units','inches','Position',[1,1,6.5,6.0]); hold on;
th = linspace(0, 2*pi, 360);
plot(rmax*cos(th), rmax*sin(th), 'k--', 'LineWidth', 1.0); % outer boundary
plot(rw*cos(th),   rw*sin(th),   'k-',  'LineWidth', 1.8); % wound boundary

for k = 1:numel(arrived_idx)
    i   = arrived_idx(k);
    ti  = trajT{i};
    xi  = trajX{i};
    yi  = trajY{i};
    mask = ti <= Tcut + 1e-9;
    if any(mask)
        ci = mapColor(arrival_time(i));
        plot(xi(mask), yi(mask), '-', 'LineWidth', 0.8, 'Color', ci);
    end
end

axis equal; xlim([-rmax rmax]); ylim([-rmax rmax]);
xlabel('X [\mum]'); ylabel('Y [\mum]');
title(sprintf('Trajectories of cells that reached the wound by T = %d min', Tcut));
grid on;

% Colorbar explaining color = arrival time to wound
cb = colorbar; colormap(cmap);
if tmax > tmin
    caxis([tmin tmax]);
else
    caxis([tmin tmin+1]); % degenerate case if all arrive same time
end
ylabel(cb, 'Arrival time to wound [min]');
text(-0.96*rmax, 0.92*rmax, ...
     {'Color encodes first-arrival time to the wound:', ...
      sprintf('earlier (dark) \x2192 later (bright), %.1f–%.1f min', tmin, tmax)}, ...
     'FontSize',9,'BackgroundColor','w','Margin',4);

saveas(gcf, sprintf('trajectories_arrivers_full_%dmin.png', Tcut));

% (B) r(t) vs time for all arrivers (same colors, using per-cell times)
figure('Color','w','Units','inches','Position',[1,1,7.0,5.0]); hold on;
for k = 1:numel(arrived_idx)
    i   = arrived_idx(k);
    ti  = trajT{i};
    xi  = trajX{i};
    yi  = trajY{i};
    mask = ti <= Tcut + 1e-9;
    if any(mask)
        r_i = sqrt(xi(mask).^2 + yi(mask).^2);
        ci  = mapColor(arrival_time(i));
        plot(ti(mask), r_i, '-', 'LineWidth', 0.8, 'Color', ci);
    end
end
yline(rw,'k-','LineWidth',1.2,'DisplayName','Wound radius r_w');
xlabel('Time [min]'); ylabel('r(t) [\mum]');
title(sprintf('Radius vs Time for cells that reached wound by %d min', Tcut));
grid on; cb = colorbar; colormap(cmap); 
if tmax>tmin, caxis([tmin tmax]); else, caxis([tmin tmin+1]); end
ylabel(cb,'Arrival time to wound [min]');
saveas(gcf, sprintf('r_vs_time_arrivers_%dmin.png', Tcut));

% (C) x(t) and y(t) vs time (same colors) — two aligned panels
figure('Color','w','Units','inches','Position',[1,1,7.5,7.0]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
for k = 1:numel(arrived_idx)
    i   = arrived_idx(k);
    ti  = trajT{i};
    xi  = trajX{i};
    mask = ti <= Tcut + 1e-9;
    if any(mask)
        ci = mapColor(arrival_time(i));
        plot(ti(mask), xi(mask), '-', 'LineWidth', 0.8, 'Color', ci);
    end
end
ylabel('x(t) [\mum]'); title('x-coordinate vs Time (arrivers)'); grid on;

nexttile; hold on;
for k = 1:numel(arrived_idx)
    i   = arrived_idx(k);
    ti  = trajT{i};
    yi  = trajY{i};
    mask = ti <= Tcut + 1e-9;
    if any(mask)
        ci = mapColor(arrival_time(i));
        plot(ti(mask), yi(mask), '-', 'LineWidth', 0.8, 'Color', ci);
    end
end
xlabel('Time [min]'); ylabel('y(t) [\mum]'); title('y-coordinate vs Time (arrivers)'); grid on;

% one shared colorbar for both subplots
h = colorbar('Position',[0.92 0.11 0.02 0.78]); colormap(cmap);
if tmax>tmin, caxis([tmin tmax]); else, caxis([tmin tmin+1]); end
ylabel(h,'Arrival time to wound [min]');
saveas(gcf, sprintf('xy_vs_time_arrivers_%dmin.png', Tcut));


disp('Simulation complete with reservoir replenishment, trajectory saving, overlays, and time-series panels (publication-ready).');
disp(['High-quality snapshots saved in folder: ', snap_dir]);

%% ======================= Local helper (bottom) ===========================
function render_clean_snapshot(axH, rmax, rw, th_circ, xw, yw, xout, yout)
% Draw a clean publication-quality frame on axH without legend/writeups.
    cla(axH);
    plot(axH, rmax*cos(th_circ), rmax*sin(th_circ), 'k--', 'LineWidth', 1.0); hold(axH,'on');
    plot(axH, rw*cos(th_circ),   rw*sin(th_circ),   'k-',  'LineWidth', 1.5);
    plot(axH, xout, yout, 'bo', 'MarkerSize', 3, 'MarkerFaceColor','none');
    plot(axH, xw,   yw,   'ro', 'MarkerSize', 5, 'MarkerFaceColor','r');
    axis(axH,'equal'); xlim(axH,[-rmax rmax]); ylim(axH,[-rmax rmax]);
    xlabel(axH,'X [\mum]'); ylabel(axH,'Y [\mum]');
    set(axH,'FontSize',11,'LineWidth',1.2);
    box(axH,'on'); grid(axH,'off');
end
