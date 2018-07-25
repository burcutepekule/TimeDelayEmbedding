%% DEMONSTRATE HENONS MAP FNN %%
clear all;close all;clc;
sampleSize  = 1e3; % number of samples, IC included
x(1)=0;
y(1)=0;
a=1.4;
b=0.3;
for i=2:sampleSize
    x(i)=1-1.4*(x(i-1)^2)+y(i-1);
    y(i)=b*x(i-1);
end
%% EMBEDDING FOR m=1 and m=2, NEIGHBORS ARE DIFFERENT
pickNNs  = [128 824 566];
screensize = get(0,'ScreenSize');
sz         = [600 500];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, ~]    = tight_subplot(1,2,[.1 .1],[.1 .1],[.1 .1]);
set(ax(1),'position',[0.130 0.35 0.83 0.6])
set(ax(2),'position',[0.130 0.13 0.83 0.2])
axes(ax(1))
plot(x(1:end-1),x(2:end),'.','Color',[0.5 0.5 0.5],'MarkerSize',4)
hold on;
scatter(x(pickNNs(1)), x(pickNNs(1)+1), 50, 'b','filled')
scatter(x(pickNNs(2)), x(pickNNs(2)+1), 50, [0.1961, 0.8039, 0.1961],'filled')
scatter(x(pickNNs(3)), x(pickNNs(3)+1), 50,'r','filled')
set(gca,'XTickLabel',[]);
% xlabel('$x[n]$','interpreter','latex')
ylabel('$x[n-1]$','interpreter','latex')
grid on;
axes(ax(2))
set(0,'defaultAxesFontSize',22)
plot(x,ones(1,sampleSize),'.','Color',[0.5 0.5 0.5],'MarkerSize',4)
hold on;
scatter(x(pickNNs(1)), 1, 50, 'b','filled')
scatter(x(pickNNs(2)), 1, 50, [0.1961, 0.8039, 0.1961],'filled')
scatter(x(pickNNs(3)), 1, 50, 'r','filled')
set(gca,'YTickLabel',[]);
xlabel('$x[n]$','interpreter','latex','FontSize',22)
grid on;
set(gca,'fontsize',22)

%% TRY m=3, NEIGHBORS STAY THE SAME
screensize = get(0,'ScreenSize');
sz         = [600 500];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, ~]    = tight_subplot(1,1,[.55 .55],[.15 .05],[.05 .05]);
set(ax(1),'position',[0.13 0.1 0.85 0.85])
axes(ax(1))
plot3(x(1:end-2),x(2:end-1),x(3:end),'.','Color',[0.5 0.5 0.5],'MarkerSize',4)
hold on;
scatter3(x(pickNNs(1)), x(pickNNs(1)+1), x(pickNNs(1)+2), 50, 'b','filled')
scatter3(x(pickNNs(2)), x(pickNNs(2)+1), x(pickNNs(2)+2), 50, [0.1961, 0.8039, 0.1961],'filled')
scatter3(x(pickNNs(3)), x(pickNNs(3)+1), x(pickNNs(3)+2), 50,'r','filled')
xlabel('$x[n-2]$','interpreter','latex')
ylabel('$x[n-1]$','interpreter','latex')
zlabel('$x[n]$','interpreter','latex')
grid on;
set(gca,'fontsize',22)

