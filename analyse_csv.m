%% Loading the data
data1 = readmatrix('1a.csv');
data2 = readmatrix('1b.csv');
data3 = readmatrix('1c.csv');
data4 = readmatrix('2a.csv');
data5 = readmatrix('2b.csv');
data6 = readmatrix('2c.csv');
data7 = readmatrix('3a.csv');
data8 = readmatrix('3b.csv');
data9 = readmatrix('3c.csv');

%% Sorting the data
x1 = data1(:,1);
z1 = data1(:,3);

x2 = data2(:,1);
z2 = data2(:,3);

x3 = data3(:,1);
z3 = data3(:,3);

x4 = data4(:,1);
z4 = data4(:,3);

x5 = data5(:,1);
z5 = data5(:,3);

x6 = data6(:,1);
z6 = data6(:,3);

x7 = data7(:,1);
z7 = data7(:,3);

x8 = data8(:,1);
z8 = data8(:,3);

x9 = data9(:,1);
z9 = data9(:,3);

%% Plot individual results
figure;
plot(x1,z1);
hold on;
plot(x2,z2);
plot(x3,z3);
title('Breakwater A Influence on Water Height');
    xlabel('Distance Across Transect');
    ylabel('Water Height (m)');
    legend({'Transect 1','Transect 2','Transect 3'})
    hold off;

figure;
plot(x4,z4);
hold on;
plot(x5,z5);
plot(x6,z6);
title('Breakwater B Influence on Water Height');
    xlabel('Distance Across Transect');
    ylabel('Water Height (m)');
    legend({'Transect 1','Transect 2','Transect 3'})
    hold off;

figure;
plot(x7,z7);
hold on;
plot(x8,z8);
plot(x9,z9);
title('Breakwater C Influence on Water Height');
    xlabel('Distance Across Transect');
    ylabel('Water Height (m)');
    legend({'Transect 1','Transect 2','Transect 3'})
    hold off;

%% Calculate and plot mean results
amean = (data1 + data2 + data3)/3;
ameanx = amean(:,1);
ameanz = amean(:,3);
figure;
plot(ameanx,ameanz,'r-');
hold on;

bmean = (data4 + data5 + data6)/3;
bmeanx = bmean(:,1);
bmeanz = bmean(:,3);
plot(bmeanx,bmeanz,'g-');

cmean = (data7 + data8 + data9)/3;
cmeanx = cmean(:,1);
cmeanz = cmean(:,3);
plot(cmeanx,cmeanz,'b-');
hold off;
title('Mean Water Height by Breakwater Type');
    xlabel('Distance Across Transect');
    ylabel('Water Height (m)');
    legend({'Type A','Type B','Type C'})


