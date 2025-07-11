%% XBeach advanced course: non-hydrostatic modelling
% Modeling of nature-based flood protection solutions with XBeach-NH
% v1.0  de Ridder     Nov-18
% v2.0 Midkiff        Oct-24
clear all
close all

%% 0. Input
sim         = 'C:\Users\middl\OneDrive\Documents\Desktop\Biloxi\Models\Whole\';

% load model
bed     = load([sim 'GIS_import_1.txt']);
x       = load([sim 'x.grd']);
y       = load([sim 'y.grd']);

% Resize bed to match the dimensions of x and y (2965 x 3408)
bed_resized = interp2(bed, linspace(1, size(bed, 2), size(x, 2)), linspace(1, size(bed, 1), size(x, 1))');

% Now plot with the original x, y, and resized bed
figure;
surf(x, y, bed_resized);  % Plot the bed as a 3D surface (just to impress the people who don't know what I'm doing)
shading interp;
colorbar;
xlabel('x-axis');
ylabel('y-axis');
zlabel('Elevation (inverted)');
title('2D Bathymetry Plot');
axis tight;
view(3);

figure;
pcolor(x, y, bed_resized);  % Plot the bed
shading interp;
colorbar;
xlabel('x-axis');
ylabel('y-axis');
zlabel('Elevation');
title('2D Bathymetry Plot');

min_plant_elevation = -1;
max_plant_elevation = 2/3;
plantMax = (bed_resized >= min_plant_elevation) & (bed_resized <= max_plant_elevation); % Selects all of the areas with elevations that are appropriate for mangroves
figure; pcolor(plantMax); shading interp; % Plots all of those locations to setup the next command

prompt = "x_min: ";
x_min = input(prompt);

prompt = "x_max: ";
x_max = input(prompt);

prompt = "y_min: ";
y_min = input(prompt);

prompt = "y_max: ";
y_max = input(prompt);

if x_min < min(x(:)) || x_max > max(x(:)) || y_min < min(y(:)) || y_max > max(y(:))
    error('Input bounds are outside the range of the grid.');
end

x_locations = (x >= x_min) & (x <= x_max); % Allows us to set the minimum and maximum x values for the bounds of our vegetation
y_locations = (y >= y_min) & (y <= y_max); % Allows us to set the minimum and maximum y values for the bounds of our vegetation
rect_mask = x_locations & y_locations; % Uses the two previous lines to make a rectangular crop of our final vegetation area

ind = plantMax & rect_mask; % All of the vegetation locations

figure;
pcolor(x, y, bed_resized);  % Plot the bed as a 2D color map
shading interp;             % Smooth out the colors for better visualization
colorbar;                   % Add a color bar to show the scale of bed elevation
xlabel('x-axis');
ylabel('y-axis');
title('Cropped Region on 2D Bathymetry Plot');

% Hold on to the current plot so the rectangle can be added
hold on;

% Add a rectangle to show the cropping area (based on your x_min, x_max, y_min, y_max)
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], 'EdgeColor', 'r', 'LineWidth', 2);

% Optional: Adjust axis scaling if needed
axis tight;  % Adjust axes to fit the data tightly

figure;
subplot(1, 2, 1);
imagesc(plantMax);
title('plantMax');

subplot(1, 2, 2);
imagesc(rect_mask);
title('rect_mask');

pointPlant = sum(plantMax(:))
pointMask = sum(rect_mask(:))
pointTotal = sum(ind(:))

%% 1. Vegetation
% create vegmap
veg = bed_resized *0;
veg(ind) = 1;

figure;
imagesc(veg);
title('Vegetation Map');

figure;
pcolor(x, y, bed_resized);        % Plot the bathymetry
shading interp;
hold on;

pcolor(x, y, veg);  % Plot the vegetation
shading interp;
alpha(0.7);  % Adjust transparency for the vegetation layer

xlabel('x-axis');
ylabel('y-axis');
title('2D Bathymetry with Vegetation Map');

assert(size(x, 1) == size(y, 1) && size(y, 2) == size(bed_resized, 2), 'Grid dimensions do not match!'); % Lets you know if your grid is funky

%% 2. Write the veg to a vegmap file
% Open the file for writing
fileID = fopen([sim '\vegmap.txt'], 'w');

% Loop through each row of the veg matrix
for i = 1:size(veg, 1)
    % Write the current row, space-separated
    fprintf(fileID, '%d ', veg(i, :));  % Write entire row in one go
    fprintf(fileID, '\n');              % Add a newline at the end of each row
end

fprintf('Number of vegetation points: %d\n', sum(veg(:)));

% Close the file
fclose(fileID);

