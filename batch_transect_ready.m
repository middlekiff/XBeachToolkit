clear all
close all
clc

% Adds OET
addpath 'C:\Users\middl\OneDrive\Documents\Desktop\CEDS DEMO - Copy\OET'; % Provides additional functions
data = 'C:\Users\middl\OneDrive\Documents\Desktop\CEDS DEMO - Copy\'; % Where the supplemental data is
cd(data) % Sets current working directory to data folder. This is where the new folders will be nested


%% Data Import
% Sets model directory and output
bathy = readmatrix([data 'breakwater.csv']);
x1 = bathy(:, 3);
y1 = bathy(:, 4);
z1 = bathy(:, 2);

% Making sure the grid aligns to the original data
xq = unique(x1);
yq = unique(y1);
[Xq, Yq] = meshgrid(xq, yq);

% Interpolate Z values onto grid
Zq = griddata(x1, y1, z1, Xq, Yq, 'cubic');


%% Cleaning the elevation data
% Remove the "%" signs at the start of lines 30-43 if you need to trim
% elevations. Otherwise the following code will be ignored.

%min = -12
%max = 0
%mask = (Zq < min | Zq > max); % Find cells outside the valid range


%% Calculating new values for outliers or NaN cells
% Preallocate local means matrix
%local_means = NaN(size(Zq));
% Replace outliers with local means
%Zq(mask) = local_means(mask);
% Verify no NaNs remain
%if any(isnan(Zq(:)))
%    Zq = fillmissing(Zq, 'nearest', 'EndValues', 'nearest');
%end


%% Parameter Setup
bathy = Zq/3.28; % Converting elevation from feet to meters- you won't always use this
Zq = bathy;

% Plots the bathymetry
figure; pcolor(bathy); shading interp; colorbar; title('measured bathymetry'); 
figure; surf(bathy); shading interp; colorbar; title('measured bathymetry [m]'); 


%% Preparing Simulation Transects
% Start at the water side, end just past the intervention

transects = ["A","B","C"];
transect_coords = struct( ...
    "A", struct("x1", 1350, "x2", 950, "y1", 400, "y2", 400), ...
    "B", struct("x1", 1300, "x2", 900, "y1", 700, "y2", 700), ...
    "C", struct("x1", 1250, "x2", 850, "y1", 1000, "y2", 1000));


%% Plot all of the transects together
% This is for later reference
figure;
h = pcolor(bathy);
shading interp;
colorbar;
title('Measured Bathymetry');
hold on;

% Force the background (pcolor) to be skipped in the legend
set(h, 'HandleVisibility', 'off');

legendHandles = gobjects(length(transects), 1); % preallocate
legendLabels = strings(length(transects), 1);

for k = 1:length(transects)
    name = transects(k);
    coords = transect_coords.(name);

    % Generate evenly spaced points along the transect
    numberOfPoints = 100;
    xPoints = linspace(coords.x1, coords.x2, numberOfPoints);
    yPoints = linspace(coords.y1, coords.y2, numberOfPoints);
    
    % Plot the transect and store the handle
    legendHandles(k) = plot([coords.x1 coords.x2], [coords.y1 coords.y2], ...
        'LineWidth', 1.5);
    legendLabels(k) = [char(name)];
end

legend(legendHandles, legendLabels);
saveas(gcf, fullfile(data, ['transect_map.png']));


%% Loop Over Each Transect
for k = 1:length(transects)
    name = transects(k);
    coords = transect_coords.(name);

    % Generate evenly spaced points along the transect
    numberOfPoints = 100;
    xPoints = linspace(coords.x1, coords.x2, numberOfPoints);
    yPoints = linspace(coords.y1, coords.y2, numberOfPoints);

    % Interpolate Z values from surface
    bathymetryValues = interp2(Zq, xPoints, yPoints);

    % Define output directory
    outdir = fullfile(pwd, ['transect_' char(name)]);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % Saves a png for reference
    figure;
    plot(1:numberOfPoints, bathymetryValues);
    title(['Transect ' char(name)]);
    xlabel('Distance');
    ylabel('Elevation (m)');
    saveas(gcf, fullfile(outdir, ['transect_' char(name) '_profile.png']));

    %% Wave and Grid Settings
    Hm0 = 4.3;
    Tp = 5.66;
    mainang = 270; % This is VERY important. The angle is the direction waves come from. 0 = Starts North, Moves South

    zs0 = 0; % storm surge level in meters
    dxmin = 1;
    nx = length(xq);
    dx = 1;
    dy = 1;

    x = (0:nx-1)*dx;
    y = 0; % 1D profile
    z = bathymetryValues;

    % Create grid
    [xgr, zgr] = xb_grid_xgrid(x, z, 'dxmin', dxmin, 'Tm', Tp/1.2, 'wl', zs0, 'nonh', 0, 'ppwl', 20);
    bathymetry = xb_grid_add('x', x, 'z', z, 'posdwn', 1);

    % Load/create jonswap wave spectrum file
    % This does apply other paramters but you hardly ever need to change
    % them. If you do, edit the other script.
    run('create_jonswap.m');

    % Create parameter settings
    pars = xb_generate_settings('xori',0,'yori',0, ...
        'wavemodel','surfbeat','morphology',0,'sedtrans',0, ...
        'wbctype','parametric','bcfile','jonswap.txt', ...
        'outputformat','netcdf', ...
        'nglobal',{'zs','u','hh','H','zb'},'tintg',60, ...
        'npoints',{'100 1','550 1'},'npointvar',{'zs','u','zb'},'tintp',10, ...
        'nmeanvar',{'zs','u','H'},'tintm',60, ...
        'order',1, ...
        'zs0',zs0, ...
        'dy',dy,...
        'rt', 7200,'tstart',0,'tstop',7200);

    % Merge input structure
    xbm_si = xs_join(bathymetry, pars);

    % Write params file
    xb_write_input(outdir, xbm_si);

    % Copy XBeach files
    copyfile(fullfile(pwd, 'xbeach.exe'), outdir);
    copyfile(fullfile(pwd, 'netcdf.dll'), outdir);

    % Run the simulation
    cd(outdir);
    system('.\xbeach.exe &');  % or remove & if you want them to run sequentially
    cd(data);  % return to base directory

end
