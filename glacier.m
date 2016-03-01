%--------------------------------
% Simulation of an eroding valley under a glacier
%--------------------------------
% Franklin Hinckley
% 14 February 2016
%--------------------------------
%
%--------------------------------

%% Clean up
clearvars
close all
clc

%% Constants [MacGregor et al]
g = 9.81;   % standard gravity [m/s^2]
A = 2.1e-16; % Arrhenius Constant [Pa^-3*yr^-1]
gamma = 0.01; % net mass balance gradient [m/yr/m]
rhoI = 917; % density of ice [kg/m^3]
rhoW = 1000; % density of water [kg/m^3]
usl0 = 0.0012; % sliding coefficient [m/yr/Pa]
edot0 = 0.0001; % erosion coefficient []
bmax = 3; % maximum accumulaton [m/yr]

%% Set up initial valley shape
dx = 1000; % [m]
len = 200*1000; % [m]
x = 0:dx:len;

SR = 20/1000;  % [m/km]

zMax = 3000; % [m]
zR = zMax - SR*x;

%% Time array
dt = 1/365/2; % [yr]
tSim = 1000; % [yr]
t = 0:dt:tSim;

P_ELA = 500; % period for variations in ELA [yr]
P_W = 1; % period for variations in water table [yr]

%% Initial glacier shape
H = zeros(size(x));

%% Allocate output
% Specify how often to save
saveInd = 250;

% Allocate output
zR_S = zeros(length(x),floor(length(t)/saveInd));
H_S = zeros(length(x),floor(length(t)/saveInd));
usl_S = zeros(length(x)-1,floor(length(t)/saveInd));
Q_S = zeros(length(x)+1,floor(length(t)/saveInd));
ELA_S = zeros(floor(length(t)/saveInd),1);
V = zeros(floor(length(t)/saveInd),1);

% Assign initial values
zR_S(:,1) = zR;
H_S(:,1) = H;
Q = zeros(length(x)+1,1);

%% Main loop
% Initialize counter for saving output
jj = 1;

% Loop
for ii = 1:length(t)
    % Compute glacier surface elevation 
    z = zR + H;
        
    % Compute current ELA
    ELA = 2500;% + 500*sin((2*pi/P_ELA)*t(ii));
    
    % Evaluate net accumulation/ablation 
    b = gamma*(z - ELA);
    
    % Cap accumulation
    b(b > bmax) = bmax;
    
    % Compute ice surface slope
    dzdx = diff(z)/dx;
    
    % Get box-centered heights
    Hm = (H(1:end-1) + H(2:end))/2;
    
    % Compute basal shear stress
    dzRdx = diff(zR)/dx;
    tauB = rhoI*g*Hm.*dzRdx;
    
    % Compute water table level
    Hwtl = 75 + 5*sin((2*pi/P_W)*t(ii));
    Hw = Hm - Hwtl;
    Hw = Hw.*(Hw > 0);
    
    % Compute sliding speed  
    Ne = rhoI*g*Hm - rhoW*g*Hw;
    usl = zeros(size(Hm));
    usl(Hm > 0) = (usl0*tauB(Hm > 0).^2)./Ne(Hm > 0);
    
    % Compute flux
    Q(2:end-1) = usl.*Hm - A*((rhoI*g*dzdx).^3).*((Hm.^5)/5);
    
    % Error trap (stops simulation if there is a numeric crash)
    if any(isnan(Q))
        error('NaN in Q')
    end
    
    % Compute flux gradient 
    dQdx = diff(Q)/dx; 
    
    % Compute thickness rate 
    dHdt = b - dQdx';
    
    % Compute erosion rate
    eDot = edot0 * usl;

    % Update glacier thickness
    H = H + dHdt*dt;
    
    % Remove negative thickness
    H = H.*(H > 0);

    % Update valley 
    zR = zR - [0 eDot*dt];
    
    % Check if save point
    if mod(ii,saveInd) == 1
        % Compute ice volume and save
        V(jj) = trapz(x,H);

        % Save output 
        ELA_S(jj)   = ELA;
        zR_S(:,jj)  = zR; % rock elevation [m]
        H_S(:,jj)   = H; % glacier thickness [m]
        Q_S(:,jj)   = Q'; % discharge [m^3/yr]
        usl_S(:,jj) = usl';
        
        % Increment counter
        jj = jj + 1;
    end
    
    % Update progress bar
    if mod(ii,100) == 1
        progressbar(ii/length(t))
    end
end

% Clean up progress bar
progressbar(1)

%% Plots
% Animation
figure
M = [];
for ii = 1:floor(length(t)/saveInd)
    % Find ELA position
    ELAind = find(zR_S(:,ii)+H_S(:,ii) < ELA_S(ii),1,'first');
    ELApos = x(ELAind)/1000;
    
    % Plot glacier
    subplot(2,1,1)
    title('Glacier','Fontsize',14)
    % Glacier
    fill([zeros(size(x)) x/1000],...
        [zeros(size(x)) zR_S(:,ii)'],'k')
    hold on
    % Bed
    fill([x/1000 fliplr(x/1000)],...
        [zR_S(:,ii)' fliplr(zR_S(:,ii)'+H_S(:,ii)')],'c')
    % ELA elevation
    plot([min(x/1000) max(x/1000)],[ELA_S(ii) ELA_S(ii)],'--b')
    % ELA position marker
    plot([ELApos ELApos],[0 4200],'--b')
    hold off
    % Write current ELA on plot
    tH = text('String',['ELA: ' num2str(ELA_S(ii),4) ' m']);
    tH.Position = [150 3750];
    tH.Color = 'b';
    % Write current date on plot
    tH = text('String',['Time: ' num2str(t(ii)*saveInd,4) ' yr']);
    tH.Position = [150 500];
    % Set limits and axes
    ylim([0 4200])
    ylabel('Elevation [m]')
    xlabel('Position [km]')
    
    % Plot flux
    subplot(2,1,2)
    % Flux
    plot(x/1000,Q_S(2:end,ii))
    hold on
    % ELA position marker
    plot([ELApos ELApos],[0 1e6],'--b')
    hold off
    ylim([0 1e6])
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    
    % Save frame
    %M = [M getframe(gcf)];
    pause(0.01)    
end

% Make movie
% v = VideoWriter('glacierWater.m4v','MPEG-4');
% open(v)
% for ii = 1:length(M)
%     writeVideo(v,M(ii))
% end
% close(v)

% Peak sliding speed
figure
plot(t(1:saveInd:end),max(usl_S))
ylabel('Peak Sliding Speed [m/s]')
xlabel('Time [yr]')
title('Sliding Speed')

% Volume vs time
Veq = (1-(1/exp(1)))*V(end); % 1 - (1/e) of equilibrium volume
eqInd = find(V > Veq,1,'first');
teq = t(eqInd); % time to Veq
figure
hold on
plot(t(1:saveInd:end),V)
%plot([t(1) t(end)],[Veq Veq],'--k')
%plot([teq teq],[0 10e7],'--k')
hold off
xlabel('Time [yr]')
ylabel('Ice Volume [m^3]')
ylim([0 1e8])
%tH = text('String',['t_{eq}: ' num2str(teq,4) ' yr']);
%tH.Position = [teq+50 1e7];
