function smpj_plot_mov_corr(file_list)
% smpj_plot_mov_corr Plot motion correction parameters from SPM realignment
%
% Usage:
%   smpj_plot_mov_corr(file_list)
%
% Description:
%   This function plots the motion correction parameters for a series of
%   runs specified by the input `file_list`. Vertical dashed lines are added
%   to indicate the boundaries between runs.
%
% Inputs:
%   file_list: A cell array of strings, where each string is the path to a
%              file containing motion correction parameters from SPM. These
%              files typically have a prefix of `rp_` and contain six columns:
%              three for translations (in mm) and three for rotations (in radians).
%
% Plot Details:
%   - The first subplot shows translations in the X, Y, and Z directions.
%   - The second subplot shows rotations around the X, Y, and Z axes (pitch, roll, and yaw).
%
% Example:
%   file_list = {'path/to/rp_file1.txt', 'path/to/rp_file2.txt'};
%   smpj_plot_mov_corr(file_list);


mov = []; % To store concatenated movement data
verticalLinePositions = [];
currentPosition = 0; % Track the current position on the x-axis

for f = 1:length(file_list)
    movRun = dlmread(file_list{f}); % Load alignment file
    mov = [mov; movRun]; % Concatenate the file data

    nVol = length(movRun); % Length of run in volumes
    currentPosition = currentPosition + nVol; % Update the position to the end of the current segment
    
    % Store the position for the vertical line at the end of each run, except the last
    if f < length(file_list)
        verticalLinePositions(end+1) = currentPosition;
    end
end

% Now, plot the vertical lines at the stored positions
figure('Position', [100, 100, 800, 600]);
subplot(211);
plot(mov(:, 1:3)); % Plot the concatenated movement data for XYZ movements
ylabel('Translation (mm)');
hold on; % Allows plotting additional lines on top

% Loop through the stored positions and plot a dashed vertical line at each
for i = 1:length(verticalLinePositions)
    xline(verticalLinePositions(i), '--', 'Color', 'k'); % Dashed line, change 'k' to preferred color
end

yline(0, '-', 'Color', 'k')
xlim([1, length(mov)])
legend('x', 'y', 'z')

hold off

subplot(212);
plot(mov(:, 4:end)); % Plot the concatenated movement data for XYZ movements
ylabel('Rotation (rad)');
hold on; % Allows plotting additional lines on top

% Loop through the stored positions and plot a dashed vertical line at each
for i = 1:length(verticalLinePositions)
    xline(verticalLinePositions(i), '--', 'Color', 'k'); % Dashed line, change 'k' to preferred color
end

yline(0, '-', 'Color', 'k')
xlim([1, length(mov)])
legend('pitch', 'roll', 'jaw')

hold off

sgtitle('Motion correction parameters')

