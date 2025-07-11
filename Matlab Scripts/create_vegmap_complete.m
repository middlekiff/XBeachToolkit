%% XBeach advanced course: non-hydrostatic modelling
% Modeling of nature-based flood protection solutions with XBeach-NH
% v1.0  de Ridder     Nov-18
clear all
close all

%% 0. Input
sim         = 'C:\Users\middl\OneDrive\Documents\Desktop\Matlab Designer Toolbox\transect_A\';

% load model
bed     = load([sim 'bed.dep']);
x       = load([sim 'x.grd']);

figure;
plot(x,bed*-1)

prompt = "vegStart: ";
vegStart = input(prompt);

prompt = "vegEnd: ";
vegEnd = input(prompt);

%% 1. Vegetation

% create vegmap
veg = bed *0;
ind = find(x>vegStart &x<vegEnd);
veg(ind) = 1;

% plot
figure;
plot(x,bed*-1+veg);hold on;
plot(x,bed*-1);hold on;
legend('veg','zb')

%% 2. Write the veg to a vegmap file
fileID = fopen([sim '\vegmap.txt'],'w');
for ii=1:length(veg)
    fprintf(fileID, [num2str(veg(ii)) ' ']);
end
fprintf(fileID, '\n');
for ii=1:length(veg)
    fprintf(fileID, [num2str(veg(ii)) ' ']);
end
fclose(fileID);