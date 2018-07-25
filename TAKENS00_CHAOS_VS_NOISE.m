%% SIMULATE A MACKEY-GLASS SYSTEM DATA GENERATION %%
clear all;close all;clc;
set(0,'defaultAxesFontSize',20)
screensize = get(0,'ScreenSize');
sz         = [500 350];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, ~]    = tight_subplot(1,1,[.55 .55],[.15 .05],[.05 .05]);
set(ax(1),'position',[0.1 0.12 0.85 0.8])
axes(ax(1))
sampleSize  = 1e3-1; % number of samples, IC excluded
dt          = 1; % time step 
a           = 0.2;     % value for a in eq (1)
b           = 0.1;     % value for b in eq (1)
tauMc       = 30;		% delay constant in eq (1)
x0          = 0.9;		% initial condition: x(t=0)=x0
varSys      = 0.0;     % system noise
[data,~]    = mckyGlss(a,b,tauMc,x0,dt,sampleSize,varSys,'euler');
subsamp     = 10;
idxBeg      = 50;
idxEnd      = 1000;
idxs        = idxBeg:subsamp:idxEnd;
% plot((idxBeg:idxEnd)./10,data(idxBeg:idxEnd)-0.8,'r','linewidth',1)
% hold on;
scatter(idxs./10,data(idxs)-0.8,'k','filled')
grid on;
%%
clear all;close all;clc;
set(0,'defaultAxesFontSize',20)
screensize = get(0,'ScreenSize');
sz         = [500 350];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, ~]    = tight_subplot(1,1,[.55 .55],[.15 .05],[.05 .05]);
set(ax(1),'position',[0.1 0.12 0.85 0.8])
axes(ax(1))
rng(7)
data=sqrt(0.06).*randn(1,1e2);
scatter(1:length(data),data,'k','filled')
grid on;
