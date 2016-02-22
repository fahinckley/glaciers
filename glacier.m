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

%% Set up initial valley shape
dx = 1000; % [m]
len = 200*1000; % [m]
x = 0:dx:len;

SR = 20/1000;  % [m/km]

zMax = 3000; % [m]
zR = zMax - SR*x;

%% Time array
dt = 1/52/4; % [yr]
tSim = 2500; % [yr]
t = 0:dt:tSim;

P = 2500; % period for variations in ELA [yr]

%% Initial glacier shape
H = zeros(size(x));

%% Main loop
zR_S = zeros(length(x),length(t));
H_S = zeros(length(x),length(t));
Q_S = zeros(length(x)+1,length(t));
ELA_S = zeros(length(t),1);
V = zeros(length(t),1);
zR_S(:,1) = zR;
H_S(:,1) = H;
for ii = 1:length(t)
    % Compute glacier surface elevation [X]
    z = zR + H;
        
    % Compute current ELA
    ELA = 2500;% + 500*sin((2*pi/P)*t(ii));
    ELA_S(ii) = ELA;
    
    % Evaluate net accumulation/ablation [X]
    b = gamma*(z - ELA);
    
    % Compute ice surface slope
    dzdx = diff(z)/dx;
    
    % Get box-centered heights
    Hm = zeros(size(dzdx)); % box-centered surface height [m]
    for jj = 1:length(dzdx)
        Hm(jj) = mean(H(jj:jj+1));
    end
    
    % Compute basal shear stress
    dzRdx = diff(zR)/dx;
    tauB = rhoI*g*Hm.*dzRdx;
    
    % Compute sliding speed
    usl = zeros(size(Hm));
    for jj = 1:length(Hm)
        if Hm(jj) > 0
            Hw = Hm(jj) - 75;
            Ne = rhoI*g*Hm(jj) - rhoW*g*Hw;
            usl(jj) = (usl0*tauB(jj)^2)/Ne;
        end
    end
    
    % Compute flux
    Q = usl.*Hm - A*((rhoI*g*dzdx).^3).*((Hm.^5)/5);
    %Q = usl.*H(1:end-1) + A*((rho*g*dzdx).^3).*((H(1:end-1).^5)/5);
    
    % Error trap (stops simulation if there is a numeric crash)
    if any(isnan(Q))
        error('NaN in Q')
    end
    
    % Pad flux array [X]
    Q = [0 Q 0];

    % Compute flux gradient [X]
    dQdx = diff(Q)/dx; 
    
    % Compute thickness rate [X] 
    dHdt = b - dQdx;
    
    % Compute erosion rate
    eDot = edot0 * usl;

    % Update glacier thickness [X]
    H = H + dHdt*dt;
    
    % Remove negative thickness [X]
    H = H.*(H > 0);

    % Update valley 
    zR = zR - [0 eDot*dt];
    
    % Compute ice volume
    V(ii) = trapz(x,H);
    
    % Save output [X]
    zR_S(:,ii) = zR; % rock elevation [m]
    H_S(:,ii) = H; % glacier thickness [m]
    Q_S(:,ii) = Q'; % discharge [m^3/yr]
    
    % Update progress bar [X]
    if mod(ii,100)
        progressbar(ii/length(t))
    end

end

% Clean up progress bar
progressbar(1)

%% Plots
% Animation
indPlot = 250;
figure
M = [];
for ii = 1:indPlot:length(t)
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
    plot([min(x/1000) max(x/1000)],[2500 2500],'--b')
    % ELA position marker
    plot([ELApos ELApos],[0 4200],'--b')
    hold off
    % Write current ELA on plot
    tH = text('String',['ELA: ' num2str(ELA_S(ii),4) ' m']);
    tH.Position = [150 3750];
    tH.Color = 'b';
    % Write current date on plot
    tH = text('String',['Time: ' num2str(t(ii),4) ' yr']);
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
    plot([ELApos ELApos],[0 5e5],'--b')
    hold off
    ylim([0 5e5])
    ylabel('Flux [m^3/yr]')
    xlabel('Position [km]')
    
    % Save frame
    M = [M getframe(gcf)];
    pause(0.01)    
end

% Make movie
M = M(1:indPlot:end);
movie2avi(M,'glacierErosion','fps',24)

% Volume vs time
Veq = (1-(1/exp(1)))*V(end); % 1 - (1/e) of equilibrium volume
eqInd = find(V > Veq,1,'first');
teq = t(eqInd); % time to Veq
figure
hold on
plot(t,V)
%plot([t(1) t(end)],[Veq Veq],'--k')
%plot([teq teq],[0 10e7],'--k')
hold off
xlabel('Time [yr]')
ylabel('Ice Volume [m^3]')
ylim([0 10e7])
%tH = text('String',['t_{eq}: ' num2str(teq,4) ' yr']);
%tH.Position = [teq+50 1e7];
