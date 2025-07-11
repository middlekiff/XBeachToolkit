
data1 = readmatrix('1a.csv');
data2 = readmatrix('1b.csv');
data3 = readmatrix('1c.csv');
data4 = readmatrix('2a.csv');
data5 = readmatrix('2b.csv');
data6 = readmatrix('2c.csv');
data7 = readmatrix('3a.csv');
data8 = readmatrix('3b.csv');
data9 = readmatrix('3c.csv');

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

h0 = z1(290) - z1(210)
h1 = z2(290) - z2(210)
h2 = z3(290) - z3(210)
h3 = z4(290) - z4(210)
h4 = z5(290) - z5(210)
h5 = z6(290) - z6(210)
h6 = z7(290) - z7(210)

 figure('NumberTitle', 'off', 'Name', 'Hurricane Rita Storm Surge');

 subplot(7,1,1); plot(x1,z1,'r-','LineWidth',2);
 grid on; title('Hurricane Rita Storm Surge','No Mangroves');

 subplot(7,1,2); plot(x2,z2,'y-','LineWidth',2);
 grid on; title('10m of Mangroves');
 subplot(7,1,3); plot(x3,z3,'g-','LineWidth',2);
 grid on; title('20m of Mangroves');

 subplot(7,1,4); plot(x4,z4,'c-','LineWidth',2);
 ylabel('Water level'); grid on; title('30m of Mangroves');

 subplot(7,1,5); plot(x5,z5,'b-','LineWidth',2);
 grid on; title('50m of Mangroves');

 subplot(7,1,6); plot(x6,z6,'m-','LineWidth',2);
 grid on; title('75m of Mangroves');

 subplot(7,1,7); plot(x7,z7,'k-','LineWidth',2); xlabel('Distance across transect'); 
 grid on; title('100m of Mangroves');

 figure; title('Hurricane Ivan Storm Surge'); 
 plot(x1,z1,'m-','LineWidth',2); hold on; % originally red
 plot(x2,z2,'c-','LineWidth',2); % originally yellow
 % plot(x3,z3,'g-','LineWidth',2) % originally green
 % plot(x4,z4,'c-','LineWidth',2); % originally cyan
 plot(x5,z5,'g-','LineWidth',2); % originally blue
 % plot(x6,z6,'m-','LineWidth',2); % originally magenta
 plot(x7,z7,'r-','LineWidth',2) % originally black

r0 = (-1)*(h0)*(1/160)*100
r1 = (-1)*(h1)*(1/160)*100
r2 = (-1)*(h2)*(1/160)*100
r3 = (-1)*(h3)*(1/160)*100
r4 = (-1)*(h4)*(1/160)*100
r5 = (-1)*(h5)*(1/160)*100
r6 = (-1)*(h6)*(1/160)*100
