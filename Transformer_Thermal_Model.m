%% ==========================================================
% IEEE Std. C57.91-2011 Recursive Transformer Thermal Model
% Calculates:
%   1. Hot-Spot Temperature (thetaHS)
%   2. Accelerated Aging Factor (FAA)
%
% Input:
%   Pg.Data  -> Transformer loading (kVA)
%
% Output:
%   thetaHS -> Hot-spot temperature (°C)
%   FAA     -> Accelerated aging factor
%% ==========================================================

% clear
clc

%% Transformer Parameters

P_rated = 250e3;          % Rated transformer capacity (VA)

R = 5;                    % Load-loss ratio
n = 0.8;                  % Oil exponent
m = 0.8;                  % Winding exponent

dThetaTO_R = 55;          % Rated top-oil rise (°C)
dThetaH_R  = 30;          % Rated hot-spot rise (°C)

tauTO = 180;              % Top-oil time constant (minutes)
tauW  = 4;                % Winding time constant (minutes)

dt = 15;                  % 15-minute interval (96 samples/day)

%% Ambient Temperature Profile (Hourly)

thetaA_hourly = 1.75 * [...
15.0 14.2 13.6 13.3 13.2 ...
13.2 13.4 13.6 13.9 14.6 ...
15.4 16.2 18.0 21.0 23.0 ...
24.0 25.0 25.0 24.0 23.0 ...
21.0 20.0 18.0 17.5];

%% Convert Hourly Temperature to 96 Samples (15-min Resolution)

thetaA = interp1(1:24,...
                 thetaA_hourly,...
                 linspace(1,24,96),...
                 'linear');

%% Transformer Loading

DT_load = fillmissing(Pg.Data,'linear');

idx = round(linspace(1,length(DT_load),96));

DT_load = abs(DT_load(idx));

%% If Pg.Data is Active Power (W), Convert to Apparent Power
% Uncomment if necessary

powerFactor = 0.95;
DT_load = DT_load/powerFactor;

%% Preallocate

N = length(DT_load);

thetaHS = zeros(1,N);

FAA = zeros(1,N);

dThetaTO = zeros(1,N);

dThetaH = zeros(1,N);

%% Initial Conditions (Based on Initial Load)

K0 = DT_load(1)/P_rated;

dThetaTO(1) = dThetaTO_R*((K0^2*R + 1)/(R + 1))^n;

dThetaH(1) = dThetaH_R*(K0^(2*m));

thetaHS(1) = thetaA(1) + dThetaTO(1) + dThetaH(1);

FAA(1) = exp((15000/383) - (15000/(thetaHS(1)+273)));

%% Recursive IEEE Thermal Model

for k = 2:N

    K = DT_load(k)/P_rated;

    %% Ultimate Top-Oil Rise

    dThetaTO_U = dThetaTO_R*((K^2*R + 1)/(R + 1))^n;

    %% Recursive Top-Oil Temperature

    dThetaTO(k) = dThetaTO(k-1) + ...
        (dThetaTO_U - dThetaTO(k-1))...
        *(1-exp(-dt/tauTO));

    %% Ultimate Hot-Spot Rise

    dThetaH_U = dThetaH_R*(K^(2*m));

    %% Recursive Hot-Spot Rise

    dThetaH(k) = dThetaH(k-1) + ...
        (dThetaH_U - dThetaH(k-1))...
        *(1-exp(-dt/tauW));

    %% Hot-Spot Temperature

    thetaHS(k) = thetaA(k) + ...
                 dThetaTO(k) + ...
                 dThetaH(k);

    %% Accelerated Aging Factor

    FAA(k) = exp((15000/383) ...
          - (15000/(thetaHS(k)+273)));

end

    %% Cumulative Aging Factor

    cum_aging = cumsum(FAA);

%% Plot

figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 3.2, 2.4]);

ax = axes;
hold(ax,'on');

y_limits = [0 100];

% ---------- Background regions ----------
% (scaled to 96 points instead of 86400 sec)
scale = 96/86400;

% Off-peak
patch([0 21600 21600 0]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [0.85 0.92 1], 'EdgeColor','none','FaceAlpha',0.35);

patch([79201 86400 86400 79201]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [0.85 0.92 1], 'EdgeColor','none','FaceAlpha',0.35);

% Mid-peak
patch([21601 32400 32400 21601]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [1 1 0.85], 'EdgeColor','none','FaceAlpha',0.35);

patch([43081 64800 64800 43081]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [1 1 0.85], 'EdgeColor','none','FaceAlpha',0.35);

% Peak
patch([32401 43080 43080 32401]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [1 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.35);

patch([64801 79200 79200 64801]*scale, [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
      [1 0.85 0.85], 'EdgeColor','none','FaceAlpha',0.35);

x = 1:96;

% ---------- Plot ----------
yyaxis left
plot(x, thetaHS,'b','LineWidth',1.2)
ylabel('Temperature (°C)','Color','b','FontSize',8)
ax.YColor = 'b';

yyaxis right
plot(x, FAA,'r','LineWidth',1.2)
hold on
plot(x, cum_aging,'k','LineWidth',1.2)   % cumulative aging in black
ylabel('Accelerated Aging Factor','Color','r','FontSize',8)
ax.YColor = 'r';

% ---------- Axes ----------
xlabel('Hours','FontSize',8)

xlim([1 96])
xticks(linspace(1,96,25))   % 0–24 hours
xticklabels(string(0:24))
xtickangle(90)

set(gca,...
    'FontSize',7,...
    'LineWidth',0.8,...
    'Box','on',...
    'TickDir','out',...
    'TickLength',[0.015 0.015]);

% ---------- Layout ----------
set(gca,'Position',[0.13 0.22 0.72 0.7]);  % reduce width

hold off

% ---------- Export ----------
set(gcf,'Color','w');
set(gca,'LooseInset',get(gca,'TightInset'));

exportgraphics(gcf,'New_DTaging_C2_S2_Proposed.png',...
    'Resolution',400,...
    'BackgroundColor','white');