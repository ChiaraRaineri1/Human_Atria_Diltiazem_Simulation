clear all
close all
clc

% define vector of constants
Const = [1 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45 0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0.0];

% define empty vectors for the variables we want to find
Vmax_v = [];
Vrest_v = [];
Vdiff_v = [];
APD90_v = [];

% define constant settings
settings.BCL = 800;         % period of stimulation (ms)
settings.TSim = 24000;      % simulation time (ms)
settings.Dur_stim = 1;      % duration of the stimulation (ms)
settings.Amp_stim = -42;    % 1.5 * amplitude of threshold stimulation (mV)
settings.NumStim = 30;      % number of stimulations

% cycle to find the variables
for C = Const
    
    settings.storeLast = 30;     % 30 APs

    [StateControl, TiControl, currents]=maleckar_main(settings, C);
    figure(1)
    title('APs')
    plot(TiControl, StateControl(:,25), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('V [mV]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    
    settings.storeLast = 1;     % 30th AP

    [StateControl,TiControl, currents]=maleckar_main(settings, C);
    L = length(StateControl(:,25));
    
    % calculate variables
    
    Vmax = max(StateControl(:,25));

    Vrest = min(StateControl(:,25));

    Vdiff = Vmax - Vrest;

    V90 = Vmax - 0.9*Vdiff;
    
    Vmax_index = find(StateControl(:,25) == Vmax);
    t0 = TiControl(Vmax_index);
    
    idx = find_index(StateControl((Vmax_index+1):L,25), V90);
    idx90 = idx + Vmax_index;
    t90 = TiControl(idx90);
    
    APD90 = t90 - t0;
   
    % plot the last AP to see how it changes w.r.t. Const
    
    figure(2)
    title('last AP')
    plot(TiControl, StateControl(:,25), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('V [mV]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    % plot the currents to see how they change w.r.t. Const
    
    figure(3)
    title('i_C_a_L')
    plot(TiControl, currents(:,1), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(4)
    title('i_K_1')
    plot(TiControl, currents(:,2), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(5)
    title('i_r_e_l')
    plot(TiControl, currents(:,3), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(6)
    title('i_t')
    plot(TiControl, currents(:,4), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(7)
    title('i_K_r')
    plot(TiControl, currents(:,5), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(8)
    title('i_K_s')
    plot(TiControl, currents(:,6), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(9)
    title('i_N_a')
    plot(TiControl, currents(:,7), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    figure(10)
    title('i_K_u_r')
    plot(TiControl, currents(:,8), 'DisplayName', strcat('Const= ', num2str(C)));
    xlabel('time [ms]');
    ylabel('I [mA/mF]');
    grid on
    hold on
    legend('show')
    
    mycolors = jet(length(Const));
    ax = gca; 
    ax.ColorOrder = mycolors;
    
    % save variables we found
    
    Vmax_v = [Vmax_v Vmax];
    Vrest_v = [Vrest_v Vrest];
    Vdiff_v = [Vdiff_v Vdiff];
    APD90_v = [APD90_v APD90];

end


%%
% plot Vrest w.r.t. Const
figure(11)
plot(Const, Vrest_v, 'k*:', 'LineWidth', 2);
set(gca, 'xdir', 'Reverse');
title('V_r_e_s_t')
xlabel('Const [-]');
ylabel('V [mV]');

% plot Vmax w.r.t. Const
figure(12)
plot(Const, Vmax_v,'k*:', 'LineWidth', 2);
set(gca, 'xdir', 'Reverse');
title('V_m_a_x')
xlabel('Const [-]');
ylabel('V [mV]');

% plot Vdiff w.r.t. Const
figure(13)
plot(Const, Vdiff_v,'k*:', 'LineWidth', 2);
set(gca, 'xdir', 'Reverse');
title('V_d_i_f_f')
xlabel('Const [-]');
ylabel('V [mV]');

% plot APD90 w.r.t. Const
figure(14)
plot(Const, APD90_v, 'k*:', 'LineWidth', 2);
set(gca, 'xdir', 'Reverse');
title('APD_9_0')
xlabel('Const [-]');
ylabel('APD_9_0 [ms]');




