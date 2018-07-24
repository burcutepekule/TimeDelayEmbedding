%% DEMONSTRATE ONE POINT %%
clear all;close all;clc;
sampleSize  = 5*1e3-1; % number of samples, IC excluded
dt          = 0.01; % time step 
[x,y,z,t]   = eulerLorenz( 10, 8/3, 28, [0.1 0 0 ], sampleSize, dt); % simulate Lorenz System with Euler method
tau         = 5; % time delay 
xObs        = x; % observed signal -> one can add noise here to see what happens
xKeep       = x; % keep noise free signal
% D    = 2.06   % Attractor dimension of Lorenz System
% m    = 2*D+1; % theoretical minimum embedding dimension (m) requirement 
m    = 3; % For illustration purposes -> to plot in 3D
T    = 1; % Prediction time -> T=1 means predict the very next point in the time series
Nf   = floor((sampleSize+1)*.4); % Size of the fitting set
Nt   = sampleSize-Nf+1; % Size of the testing set
Ns   = 1; % Sampling rate -> Ns=1 uses all points in the time series
iSet = Nf:Ns:Nf+Nt-T; % Index set of the testing set
numOfNearNeigh = 20; %Number of nearest neighbours for regression
jSet           = 1+(m-1)*tau:Nf-T; %Index set of the fittiing set
i              = iSet(1); % Demonstrate for i=iSet(1)=Nf;
testDelayVector= xObs(i:-tau:i-(m-1)*tau); % construct the delay vector of x[i]
d_ij           = zeros(1,length(jSet)); % allocate memory for distance vector
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
%     alphas               = regress(futureVec,XMat); % regression with built-in function "regress"
alphas               = (transpose(XMat)*XMat)\transpose(XMat)*futureVec; % regression, matrix multiplication form
forecastValue        = [1 testDelayVector]*alphas; % forecasted values
forecastDelayVector  = [forecastValue xObs(i+T-tau:-tau:i-(m-1)*tau)]; % delay vector with the forecasted value to project on the attractor
futureDelayVector    = xObs(i+T:-tau:i-(m-1)*tau); % delay vector with the actual future value to project on the attractor
%% PLOTTING %
close all;clc;
set(0,'defaultAxesFontSize',20)
screensize = get(0,'ScreenSize');
sz         = [600 400];
xpos       = ceil((screensize(3)-sz(1))/2); % center the figure on the
ypos       = ceil((screensize(4)-sz(2))/2); % center the figure on the
h          = figure('Position', [xpos , ypos, sz(1), sz(2)]);
[ax, pos]  = tight_subplot(1,1,[.55 .55],[.15 .05],[.05 .05]);
cmap = summer(numOfNearNeigh); % colormap to plot nearest neighbours
gray = 0.7.*ones(1,3); % color to plot the attractor
set(ax(1),'position',[0.15 0.150 0.8 0.7])
axis(ax(1));
xReconstruct = [x(1:sampleSize-2*tau)' x(tau+1:sampleSize-tau)' x(2*tau+1:sampleSize)']; % reconstructed attractor from x(t)
plot3(xReconstruct(:,1),xReconstruct(:,2),xReconstruct(:,3),'color',gray)
hold on;
h0 = scatter3(testDelayVector(3),testDelayVector(2),testDelayVector(1),50,'k','filled'); % delayed test vector on the attractor
h1 = scatter3(futureDelayVector(3),futureDelayVector(2),futureDelayVector(1),50,[0.7 0.7 0.7],'filled'); % true value of the delayed future vector
for nIdx=1:numOfNearNeigh
    scatter3(XMat(nIdx,4),XMat(nIdx,3),XMat(nIdx,2),50,cmap(nIdx,:),'filled'); % plot the nearest neighbours on the attractor
end
h2   = scatter3(forecastDelayVector(3),forecastDelayVector(2),forecastDelayVector(1),90,'r'); % delayed forecast vector on the attractor
% plotting details : to focus on the part of the attractor to visualize the
% points of interest -> change if index i is different 
xmin = -20; ymin = -20; zmin = -20;
xmax = -10; ymax = -10; zmax = -10;
xlim([xmin xmax])
ylim([ymin ymax])
zlim([zmin zmax])
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
xlabel('$x(t)$','interpreter','latex')
ylabel('$x(t-\tau)$','interpreter','latex')
zlabel('$x(t-2\tau)$','interpreter','latex')
view([-14,25]) % rotates the 3D  plot
grid on;
leg = legend([h0 h1 h2],{'$\mathbf{x}(t)$','$\mathbf{x}(t+T)$','$\hat{\mathbf{x}}(t+T)$'},'location','northeastoutside');
set(leg,'interpreter','latex');
align_axislabel([0 8 0],h) % aligns the axis labels for a better 3D vision
