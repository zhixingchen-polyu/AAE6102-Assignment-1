% Extract user's latitude, longitude, and height
latitudes = navSolutionsCT.usrPosLLH(:, 1);  % Latitude
longitudes = navSolutionsCT.usrPosLLH(:, 2); % Longitude
heights = navSolutionsCT.usrPosLLH(:, 3);     % Height

% Create user trajectory plot and save it
figure; % Create a new figure window
plot3(longitudes, latitudes, heights, 'b-', 'LineWidth', 0.5); % 3D plot, trajectory as blue line
hold on;

scatter3(longitudes, latitudes, heights, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5); % Mark each position with red points
grid on;

% Add title and axis labels
title('User Trajectory in LLH Coordinates');
xlabel('Longitude (degrees)');
ylabel('Latitude (degrees)');
zlabel('Height (meters)');

% Adjust axis limits
%xlim([min(longitudes), max(longitudes)]);
%ylim([min(latitudes), max(latitudes)]);
%zlim([min(heights), max(heights)]);
legend('Trajectory', 'User Position');

% Add real position marker (drawn last to ensure it's on top)
realLatitude = 22.328444770087565; % Real latitude
realLongitude = 114.1713630049711;  % Real longitude
realHeight = 0;  % Assume real height is 0 (or adjust as necessary)
scatter3(realLongitude, realLatitude, realHeight, 50, 'k', 'filled', 'DisplayName', 'Real Position'); % Mark real position with black point

% Update legend to include real position
legend('Trajectory', 'User Position', 'Real Position');

% Save user trajectory plot
saveas(gcf, 'User_Trajectory.png'); % Save as PNG file
close(gcf); % Close the current figure window

% Extract velocity data
velocities = navSolutionsCT.usrVelENU; % Get all velocities

% Create user velocity plot and save it
figure; % Create a new figure window
hold on;
plot(velocities(1, :), 'g', 'DisplayName', 'East Velocity'); % East velocity
plot(velocities(2, :), 'r', 'DisplayName', 'North Velocity'); % North velocity
plot(velocities(3, :), 'b', 'DisplayName', 'Up Velocity'); % Up velocity
grid on;
title('User Velocity in ENU Coordinates');
xlabel('Samples');
ylabel('Velocity (m/s)');
legend show;
hold off;

% Save user velocity plot
saveas(gcf, 'User_Velocity.png'); % Save as PNG file
close(gcf); % Close the current figure window