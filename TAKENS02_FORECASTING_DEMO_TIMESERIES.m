%% OPTION 1 : SIMULATE A LORENZ ATTRACTOR FOR DATA GENERATION %%
clear all;close all;clc;
sampleSize  = 5*1e3-1; % number of samples, IC excluded
dt          = 0.01; % time step 
[data,y,z,t]= eulerLorenz( 10, 8/3, 28, [0.1 0 0 ], sampleSize, dt); % simulate Lorenz System with Euler method
% D    = 2.06;   % Attractor dimension of Lorenz System
% m    = ceil(2*D+1); % theoretical minimum embedding dimension (m) requirement 
%% OPTION 2 : SIMULATE A MACKEY-GLASS SYSTEM DATA GENERATION %%
clear all;close all;clc;
sampleSize  = 5*1e3-1; % number of samples, IC excluded
dt          = 0.1; % time step 
a           = 0.2;     % value for a in eq (1)
b           = 0.1;     % value for b in eq (1)
tauMc       = 30;		% delay constant in eq (1)
x0          = 0.9;		% initial condition: x(t=0)=x0
varSys      = 0.0;     % system noise
[data,~]    = mckyGlss(a,b,tauMc,x0,dt,sampleSize,varSys,'euler');
%% OPTION 3 : LOAD ONE OF THE DATASETS %%
clear all;close all;clc;
load('sunSpots.mat');data=data'; %Sunspot time series 
% load('minDailyTemp');data=data'; %Daily Temp. time series
sampleSize  = length(data)-1;
%% FORECASTING %%
tau         = 5; % time delay 
xObs        = data; % observed signal -> one can add noise here to see what happens
xKeep       = data; % keep noise free signal
m    = 3; % embedding dimension
T    = 1; % Prediction time -> T=1 means predict the very next point in the time series
Nf   = floor((sampleSize+1)*.4); % Size of the fitting set
Nt   = sampleSize-Nf+1; % Size of the testing set
Ns   = 1; % Sampling rate -> Ns=1 uses all points in the time series
iSet = Nf:Ns:Nf+Nt-T; % Index set of the testing set
numOfNearNeigh = 20; %Number of nearest neighbours for regression
err            = zeros(1,length(iSet)); % Allocate memory for error vector
forecast       = zeros(1,length(iSet)); % Allocate memory for forecast results
jSet           = 1+(m-1)*tau:Nf-T; %Index set of the fittiing set 
for iIdx=1:length(iSet) %For every point in the testing set ...
    sprintf('i : %d/%d',iIdx,length(iSet))
    i               = iSet(iIdx); 
    testDelayVector = xObs(i:-tau:i-(m-1)*tau); % construct the delay vector of x[i]
    d_ij            = zeros(1,length(jSet)); % allocate memory for distance vector
    for jIdx=1:length(jSet)
        j              = jSet(jIdx);
        fitDelayVector = xObs(j:-tau:j-(m-1)*tau); % construct the delay vector of x[j]
        d_ij(jIdx)     = norm(testDelayVector-fitDelayVector); %distances |x[i]-x[j]|
    end
    [~,distIdxs] = sort(d_ij,'ascend'); % sort the distances
    jVecSorted   = jSet(distIdxs);
    kNeigIdxs    = jVecSorted(1:numOfNearNeigh); % indexes of the nearest neighbours
    futureVec    = zeros(numOfNearNeigh,1);   % allocate memory for future values of the neighbours 
    XMat         = zeros(numOfNearNeigh,m+1); % allocate memory for padded delayed neighbour vectors
    for l=1:numOfNearNeigh
        j_l           = kNeigIdxs(l); % index of the neighbour
        futureVec(l,1)= xObs(j_l+T);  % future value (t+T) of the neighbour vector
        XMat(l,:)     = [1 xObs(j_l:-tau:j_l-(m-1)*tau)]; % delayed neighbour vectors, padded with 1s for regression 
    end
     alphas               = regress(futureVec,XMat); % regression with built-in function "regress"
%      alphas               = (transpose(XMat)*XMat)\transpose(XMat)*futureVec; % regression, matrix multiplication form
    forecast(iIdx)       = [1 testDelayVector]*alphas; % forecasted values 
    err(iIdx)            = abs(forecast(iIdx)-xKeep(i+T)); % error value for the forecast
end
errRMS  = sqrt(sum(err.^2))/std(xObs); % Normalized RMS forecasting error
%% PLOTTING %%
close all;
set(0,'defaultAxesFontSize',20)
screensize = get(0,'ScreenSize');
sz         = [800 600];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, ~]    = tight_subplot(1,2,[.55 .55],[.15 .05],[.05 .05]);
set(ax(1),'position',[0.12 0.60 0.84 0.38])
set(ax(2),'position',[0.12 0.10 0.84 0.38])
axes(ax(1))
plot(1:sampleSize+1,xKeep,'k','linewidth',2)
hold on;
plot((iSet+T),forecast,'ro--','MarkerSize',2)
grid on;
axis tight
xlabel('$n$','interpreter','latex')
ylabel('$x[n]$','interpreter','latex')
leg = legend('$x[n]$','$\hat{x}[n]$');
set(leg,'interpreter','latex')
axes(ax(2))
plot([1:Nf iSet+T],[zeros(1,Nf) err],'b','linewidth',2)
grid on;
axis tight
xlabel('$n$','interpreter','latex')
ylabel('$\Big|x[n]-\hat{x}[n]\Big|$','interpreter','latex')
