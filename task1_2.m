clear all, close all
%%
outdata = importdata('/Users/bradphilipbar/Documents/ME_561/Task1/OUTPUT.dat');
% Calculating the exact solution
nu = 1.568e-5;
H = 1;
Re = 1000;
UB = Re*nu/(2*H);
P0 = 3*nu*UB/H^2;
y=linspace(0,1,10000);
Ux = -P0/(2*nu) * (y.^2 - H^2);
Umax = P0*H^2/(2*nu);
Unorm = Ux/Umax;
% plot vel
plot(y,Unorm,'r-',outdata(:,1),outdata(:,2),'k.','LineWidth',2,'MarkerSize',12)
legend('exact','numerical')
ylabel('U_x / U_{max}')
xlabel('y/H')
figure
% plot res
resdata = importdata('/Users/bradphilipbar/Documents/ME_561/Task1/RES.dat');
semilogy(resdata(:,1),resdata(:,2),'k.-','LineWidth',2,'MarkerSize',15)
ylabel('Residual')
xlabel('Time')




