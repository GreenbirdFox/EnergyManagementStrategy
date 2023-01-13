clc;
clear
test_data = xlsread('Test_data.xlsx');
t = test_data(:,1);      % [s] Time
v = test_data(:,2);      % [mph] Velocity
v = v*1.61;              % Velocity mph to km/h
P_dem = cal_P(v);        % Power demand
ts = 1;                  % Time step
N = length(t);           % Length of time vector

% Input data plots
subplot(2,1,1);
plot(t,P_dem)
xlabel('Time [s]'); 
ylabel('Power [kW]'); 
title('Power demand');  
subplot(2,1,2);
plot(t,v);
xlabel('Time [s]'); 
ylabel('Vehicle Speed [km/h]'); 
title('Desired velocity');

f1_wt_en = 0.001;       % No. of grams consumed per unit energy consumption
Pe_max = 30000;         % [W] Maximum engine power
Pb_max = 15000;         % [W] Maximum battery power
Q_batt = 18900;         % [As] Battery capacity
U_oc = 320;             % [V] Open circuit voltage of the battery
SOC_min = 0.2;          % Lower SOC limit
SOC_max = 0.9;          % Upper SOC limit

SOC_grid = linspace(SOC_min,SOC_max,100);
ns = length(SOC_grid);

V = zeros(ns, N);       % Value function
V(:, N) = 0;            % Boundary condition
% Backward recursion / Iterate throught time vector
for i = N-1:-1:1
    % Iterate through SOC grid
    for j = 1:ns
        % Constraints
        lb = max([(((SOC_max-SOC_grid(j))*Q_batt*U_oc)/-ts), -Pb_max,  ...
            P_dem(i)-Pe_max]);  % lower bound P_batt
        ub = min([(((SOC_min-SOC_grid(j))*Q_batt*U_oc)/-ts), Pb_max,  ...
            P_dem(i)]);         % Upper bound P_batt
        % 250 Feasible actions
        P_batt_grid = linspace(lb,ub,250);
        P_eng = P_dem(i)-P_batt_grid;
        % State-Transition functions
        if SOC_grid(j)>=SOC_max-0.567
            % Decision variable a=0 was added
            SOC_next = SOC_grid(j);
        elseif (SOC_min<=SOC_grid(j)) && (SOC_grid(j)<=SOC_max-0.567)
            % Decision variable b was added
            SOC_next = SOC_grid(j)+0.294; 
        end
        % Decision policy
        V_nxt = interp1(SOC_grid, V(:,i+1), SOC_next);
        % Value function/cost-to-go
        c2g = (ts*f1_wt_en*P_eng)./1;
        [V(j,i), k] = min([c2g + V_nxt]);
        u_opt(j,i) = P_batt_grid(k);
    end
end

% Try with initial SOC 80%
[Pb_08, Pe_08, FC_08, SOC_08] = RUNHEV(.8,N,SOC_grid,u_opt,P_dem);
% Objective-Contribution Functions
% Calculations are in RUNHEV function
Consumption = zeros(1, N);
for i = 1:N-1
    if FC_08(i) >= 0
        Consumption = [Consumption(i) cumsum(FC_08(1:i))];
    else
        FC_08(i) = 0;
        Consumption = [Consumption(i) cumsum(FC_08(1:i))];
    end
    Consumption(1) = 0;
end

% SOC figure
figure;
plot(SOC_08);
hold on;
plot(linspace(SOC_min,SOC_min, length(t)));
plot(linspace(SOC_max,SOC_max, length(t)));
xlabel('Time [s]'); 
ylabel('SOC [%]'); 
title('SOC');
legend('SOC 0.8', 'SOC_min', 'SOC_max');

% Consumption figure
figure;
plot(Consumption);
xlabel('Time [s]'); 
ylabel('Cumulation Fuel Comsumption');
title('Comsumption');

function p_dem = cal_P(v)
    m=10059;    
    g=9.8;
    f=0.01;
    CD=0.7;
    A=8.6;
    delta=1.1;
    coef=0.9;
    n=length(v); 
    for i=1:n-1
       dv(i)=(v(i+1)-v(i))/3.6;
    end
    dvdt=[0 dv]';
    % Calculate power
    p=v/3600/coef.*(m*g*f+CD*A*v.^2/21.15+delta*m*dvdt); 
    p_dem=p*1000;       % P_dem kW to W
end


function [Pb_act, Pe_act, FC_act, SOC_act]=RUNHEV(SOC0,N,SOCgrd,u_opt,Pd)
    SOC_act(1) = SOC0;
    SOC_grid = SOCgrd;
    P_dem = Pd;
    U_oc = 320;
    ts = 1;
    Q_batt = 18900;
    fl_wt_en = 0.001;
    % Forward loop
    for i = 1:N-1
        Pb_act(i) = interp1(SOC_grid, u_opt(:,i), SOC_act(i));
        Pe_act(i) = P_dem(i) - Pb_act(i);
        % Assume Efficiency = 100%
        % Objective-Contribution Functions
        FC_act(i) = (ts*fl_wt_en*Pe_act(i))./1;
        SOC_act(i+1) = SOC_act(i) - ((ts*Pe_act(i))/(Q_batt*U_oc));
    end
end
