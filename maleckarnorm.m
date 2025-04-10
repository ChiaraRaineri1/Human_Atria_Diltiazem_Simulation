%===============================================================================
% CellML file:   C:\Documents and Settings\usuario\Escritorio\Modelos Cellml Modificados\maleckar_normalizado2.cellml
% CellML model:  Maleckar
% Date and time: 26/03/2012 at 13:14:51
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2012 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================

function dY = maleckarnorm(time, Y, flag, settings, Ns, Const)
time_real=(Ns-1)*settings.BCL+time;

%dummy=X(7);
%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [0.630471, 0.646226, 0.45453, 0.002665, 0.43071, 0.001098, 0.948202, 0.000014, 0.998578, 0.998561, 1.814418, 5.588239, 130.019282, 0.004661, 0.000054, 1.38222, 0.026604, 0.012843, 0.190077, 0.714719, 7.1e-5, 6.5e-5, 129.502075, 8.488527, -73.941851, 0.875262, 0.870692, 0.003325, 0.000371, 0.966869];
% 
% YNames = {'Ca_rel', 'Ca_up', 'F1', 'F2', 'O_Calse', 'r', 's', 'd_L', 'f_L1', 'f_L2', 'Ca_c', 'K_c', 'Na_c', 'n', 'pa', 'O', 'O_C', 'O_TC', 'O_TMgC', 'O_TMgMg', 'Ca_d', 'Ca_i', 'K_i', 'Na_i', 'V', 'h1', 'h2', 'm', 'a_ur', 'i_ur'};
% YUnits = {'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'millimolar', 'millimolar', 'millimolar', 'millimolar', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless'};
% YComponents = {'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_handling_by_the_SR', 'Ca_independent_transient_outward_K_current_r_gate', 'Ca_independent_transient_outward_K_current_s_gate', 'L_type_Ca_channel_d_L_gate', 'L_type_Ca_channel_f_L1_gate', 'L_type_Ca_channel_f_L2_gate', 'cleft_space_ion_concentrations', 'cleft_space_ion_concentrations', 'cleft_space_ion_concentrations', 'delayed_rectifier_K_currents_n_gate', 'delayed_rectifier_K_currents_pa_gate', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_Ca_buffering', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'intracellular_ion_concentrations', 'membrane', 'sodium_current_h1_gate', 'sodium_current_h2_gate', 'sodium_current_m_gate', 'ultra_rapid_K_current_aur_gate', 'ultra_rapid_K_current_iur_gate'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Ca_rel (millimolar) (in Ca_handling_by_the_SR)
% 2: Ca_up (millimolar) (in Ca_handling_by_the_SR)
% 3: F1 (dimensionless) (in Ca_handling_by_the_SR)
% 4: F2 (dimensionless) (in Ca_handling_by_the_SR)
% 5: O_Calse (dimensionless) (in Ca_handling_by_the_SR)
% 6: r (dimensionless) (in Ca_independent_transient_outward_K_current_r_gate)
% 7: s (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
% 8: d_L (dimensionless) (in L_type_Ca_channel_d_L_gate)
% 9: f_L1 (dimensionless) (in L_type_Ca_channel_f_L1_gate)
% 10: f_L2 (dimensionless) (in L_type_Ca_channel_f_L2_gate)
% 11: Ca_c (millimolar) (in cleft_space_ion_concentrations)
% 12: K_c (millimolar) (in cleft_space_ion_concentrations)
% 13: Na_c (millimolar) (in cleft_space_ion_concentrations)
% 14: n (dimensionless) (in delayed_rectifier_K_currents_n_gate)
% 15: pa (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
% 16: O (dimensionless) (in intracellular_Ca_buffering)
% 17: O_C (dimensionless) (in intracellular_Ca_buffering)
% 18: O_TC (dimensionless) (in intracellular_Ca_buffering)
% 19: O_TMgC (dimensionless) (in intracellular_Ca_buffering)
% 20: O_TMgMg (dimensionless) (in intracellular_Ca_buffering)
% 21: Ca_d (millimolar) (in intracellular_ion_concentrations)
% 22: Ca_i (millimolar) (in intracellular_ion_concentrations)
% 23: K_i (millimolar) (in intracellular_ion_concentrations)
% 24: Na_i (millimolar) (in intracellular_ion_concentrations)
% 25: V (millivolt) (in membrane)
% 26: h1 (dimensionless) (in sodium_current_h1_gate)
% 27: h2 (dimensionless) (in sodium_current_h2_gate)
% 28: m (dimensionless) (in sodium_current_m_gate)
% 29: a_ur (dimensionless) (in ultra_rapid_K_current_aur_gate)
% 30: i_ur (dimensionless) (in ultra_rapid_K_current_iur_gate)

%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------



ACh = 1.0e-24;   % millimolar (in ACh_dependent_K_current)
I_up_max = 2800.0;   % picoA (in Ca_handling_by_the_SR)
Vol_rel = 0.0000441;   % nanolitre (in Ca_handling_by_the_SR)
Vol_up = 0.0003969;   % nanolitre (in Ca_handling_by_the_SR)
alpha_rel = 200000.0;   % picoA_per_millimolar (in Ca_handling_by_the_SR)
k_cyca = 0.0003;   % millimolar (in Ca_handling_by_the_SR)
k_rel_d = 0.003;   % millimolar (in Ca_handling_by_the_SR)
k_rel_i = 0.0003;   % millimolar (in Ca_handling_by_the_SR)
k_srca = 0.5;   % millimolar (in Ca_handling_by_the_SR)
k_xcs = 0.4;   % dimensionless (in Ca_handling_by_the_SR)
r_recov = 0.815;   % per_millisecond (in Ca_handling_by_the_SR)
tau_tr = 0.01;   % millisecond (in Ca_handling_by_the_SR)
g_t = 8.25;   % nanoS (in Ca_independent_transient_outward_K_current)
E_Ca_app = 60.0;   % millivolt (in L_type_Ca_channel)
% segnaposto
g_Ca_L = 6.75*Const;   % nanoS (in L_type_Ca_channel)
k_Ca = 0.025;   % millimolar (in L_type_Ca_channel)
K_NaCa = 0.0374842;   % picoA_per_millimolar_4 (in Na_Ca_ion_exchanger_current)
d_NaCa = 0.0003;   % per_millimolar_4 (in Na_Ca_ion_exchanger_current)
gamma_Na = 0.45;   % dimensionless (in Na_Ca_ion_exchanger_current)
g_B_Ca = 0.078681;   % nanoS (in background_currents)
g_B_Na = 0.060599;   % nanoS (in background_currents)
Ca_b = settings.Ca_o;   % millimolar (in cleft_space_ion_concentrations)
K_b = settings.K_o;   % millimolar (in cleft_space_ion_concentrations)
Na_b = settings.Na_o;   % millimolar (in cleft_space_ion_concentrations)
Vol_c = 0.000800224;   % nanolitre (in cleft_space_ion_concentrations)
tau_Ca = 24.7;   % millisecond (in cleft_space_ion_concentrations)
tau_K = 10.0;   % millisecond (in cleft_space_ion_concentrations)
tau_Na = 14.3;   % millisecond (in cleft_space_ion_concentrations)
g_Kr = 0.5;   % nanoS (in delayed_rectifier_K_currents)
g_Ks = 1.0;   % nanoS (in delayed_rectifier_K_currents)
Mg_i = 2.5;   % millimolar (in intracellular_Ca_buffering)
Vol_d = 0.00011768;   % nanolitre (in intracellular_ion_concentrations)
Vol_i = 0.005884;   % nanolitre (in intracellular_ion_concentrations)
phi_Na_en = 0.0;   % picoA_per_pF (in intracellular_ion_concentrations)
tau_di = 0.01;   % millisecond (in intracellular_ion_concentrations)
g_K1 = 3.1;   % nanoS (in inward_rectifier)
F = 96487.0;   % coulomb_per_mole (in membrane)
R = 8314.0;   % millijoule_per_mole_kelvin (in membrane)
T = 306.15;   % kelvin (in membrane)
i_CaP_max = 4.0;   % picoA (in sarcolemmal_calcium_pump_current)
k_CaP = 0.0002;   % millimolar (in sarcolemmal_calcium_pump_current)
P_Na = 0.0018;   % nanolitre_per_second (in sodium_current)
K_NaK_K = 1.0;   % millimolar (in sodium_potassium_pump)
i_NaK_max = 68.55;   % picoA (in sodium_potassium_pump)
pow_K_NaK_Na_15 = 36.4829;   % millimolar15 (in sodium_potassium_pump)
g_kur = 2.25;   % nanoS (in ultra_rapid_K_current)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% i_KACh (picoA_per_pF) (in ACh_dependent_K_current)
% J_O_Calse (per_millisecond) (in Ca_handling_by_the_SR)
% i_rel (picoA_per_pF) (in Ca_handling_by_the_SR)
% i_rel_f2 (dimensionless) (in Ca_handling_by_the_SR)
% i_rel_factor (dimensionless) (in Ca_handling_by_the_SR)
% i_tr (picoA_per_pF) (in Ca_handling_by_the_SR)
% i_up (picoA_per_pF) (in Ca_handling_by_the_SR)
% r_Ca_d_factor (dimensionless) (in Ca_handling_by_the_SR)
% r_Ca_d_term (dimensionless) (in Ca_handling_by_the_SR)
% r_Ca_i_factor (dimensionless) (in Ca_handling_by_the_SR)
% r_Ca_i_term (dimensionless) (in Ca_handling_by_the_SR)
% r_act (per_millisecond) (in Ca_handling_by_the_SR)
% r_inact (per_millisecond) (in Ca_handling_by_the_SR)
% r_infinity (dimensionless) (in Ca_independent_transient_outward_K_current_r_gate)
% tau_r (millisecond) (in Ca_independent_transient_outward_K_current_r_gate)
% s_factor (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
% s_infinity (dimensionless) (in Ca_independent_transient_outward_K_current_s_gate)
% tau_s (millisecond) (in Ca_independent_transient_outward_K_current_s_gate)
% E_K (millivolt) (in Ca_independent_transient_outward_K_current)
% i_t (picoA_per_pF) (in Ca_independent_transient_outward_K_current)
% d_L_factor (dimensionless) (in L_type_Ca_channel_d_L_gate)
% d_L_infinity (dimensionless) (in L_type_Ca_channel_d_L_gate)
% tau_d_L (millisecond) (in L_type_Ca_channel_d_L_gate)
% f_L_factor (millivolt) (in L_type_Ca_channel_f_L1_gate)
% f_L_infinity (dimensionless) (in L_type_Ca_channel_f_L1_gate)
% tau_f_L1 (millisecond) (in L_type_Ca_channel_f_L1_gate)
% tau_f_L2 (millisecond) (in L_type_Ca_channel_f_L2_gate)
% f_Ca (dimensionless) (in L_type_Ca_channel)
% i_Ca_L (picoA_per_pF) (in L_type_Ca_channel)
% i_NaCa (picoA_per_pF) (in Na_Ca_ion_exchanger_current)
% E_Ca (millivolt) (in background_currents)
% i_B_Ca (picoA_per_pF) (in background_currents)
% i_B_Na (picoA_per_pF) (in background_currents)
% n_factor (dimensionless) (in delayed_rectifier_K_currents_n_gate)
% n_infinity (dimensionless) (in delayed_rectifier_K_currents_n_gate)
% tau_n (millisecond) (in delayed_rectifier_K_currents_n_gate)
% p_a_infinity (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
% pa_factor (dimensionless) (in delayed_rectifier_K_currents_pa_gate)
% tau_pa (millisecond) (in delayed_rectifier_K_currents_pa_gate)
% pip (dimensionless) (in delayed_rectifier_K_currents_pi_gate)
% factorAZD (dimensionless) (in delayed_rectifier_K_currents)
% i_Kr (picoA_per_pF) (in delayed_rectifier_K_currents)
% i_Ks (picoA_per_pF) (in delayed_rectifier_K_currents)
% J_O (per_millisecond) (in intracellular_Ca_buffering)
% J_O_C (per_millisecond) (in intracellular_Ca_buffering)
% J_O_TC (per_millisecond) (in intracellular_Ca_buffering)
% J_O_TMgC (per_millisecond) (in intracellular_Ca_buffering)
% J_O_TMgMg (per_millisecond) (in intracellular_Ca_buffering)
% i_di (picoA_per_pF) (in intracellular_ion_concentrations)
% i_K1 (picoA_per_pF) (in inward_rectifier)
% I (picoA_per_pF) (in membrane)
% Q_tot (millivolt) (in membrane)
% i_Stim (picoA_per_pF) (in membrane)
% past (millisecond) (in membrane)
% i_CaP (picoA_per_pF) (in sarcolemmal_calcium_pump_current)
% h_factor (dimensionless) (in sodium_current_h1_gate)
% h_infinity (dimensionless) (in sodium_current_h1_gate)
% tau_h1 (millisecond) (in sodium_current_h1_gate)
% tau_h2 (millisecond) (in sodium_current_h2_gate)
% m_factor (dimensionless) (in sodium_current_m_gate)
% m_infinity (dimensionless) (in sodium_current_m_gate)
% tau_m (millisecond) (in sodium_current_m_gate)
% E_Na (millivolt) (in sodium_current)
% i_Na (picoA_per_pF) (in sodium_current)
% i_NaK (picoA_per_pF) (in sodium_potassium_pump)
% pow_Na_i_15 (millimolar15) (in sodium_potassium_pump)
% a_ur_infinity (dimensionless) (in ultra_rapid_K_current_aur_gate)
% tau_a_ur (millisecond) (in ultra_rapid_K_current_aur_gate)
% i_ur_infinity (dimensionless) (in ultra_rapid_K_current_iur_gate)
% tau_i_ur (millisecond) (in ultra_rapid_K_current_iur_gate)
% i_Kur (picoA_per_pF) (in ultra_rapid_K_current)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

E_K = R*T/F*log(Y(12)/Y(23));
i_KACh = 1.0/50.0*10000.0/(1.0+9.13652*1.0^0.477811/ACh^0.477811)*(0.0517+0.4516/(1.0+exp((Y(25)+59.53)/17.18)))*(Y(25)-E_K)*0.05;
r_Ca_d_term = Y(21)/(Y(21)+k_rel_d);
r_Ca_i_term = Y(22)/(Y(22)+k_rel_i);
r_Ca_d_factor = r_Ca_d_term*r_Ca_d_term*r_Ca_d_term*r_Ca_d_term;
r_Ca_i_factor = r_Ca_i_term*r_Ca_i_term*r_Ca_i_term*r_Ca_i_term;
r_act = 203.8*(r_Ca_i_factor+r_Ca_d_factor);
r_inact = 33.96+339.6*r_Ca_i_factor;
dY(3, 1) = (r_recov*(1.0-Y(3)-Y(4))-r_act*Y(3))/1000.0;
dY(4, 1) = (r_act*Y(3)-r_inact*Y(4))/1000.0;
i_rel_f2 = Y(4)/(Y(4)+0.25);
i_rel_factor = i_rel_f2*i_rel_f2;
i_rel = alpha_rel*i_rel_factor*(Y(1)-Y(22))/50.0;
i_up = I_up_max/50.0*(Y(22)/k_cyca-k_xcs*k_xcs*Y(2)/k_srca)/((Y(22)+k_cyca)/k_cyca+k_xcs*(Y(2)+k_srca)/k_srca);
i_tr = (Y(2)-Y(1))*2.0*Vol_rel*F/(tau_tr*50.0);
J_O_Calse = 480.0*Y(1)*(1.0-Y(5))-400.0*Y(5);
dY(5, 1) = J_O_Calse/1000.0;
dY(2, 1) = 50.0*(i_up-i_tr)/(2.0*Vol_up*F)/1000.0;
dY(1, 1) = (50.0*(i_tr-i_rel)/(2.0*Vol_rel*F)-31.0*J_O_Calse)/1000.0;
i_t = g_t*Y(6)*Y(7)*(Y(25)-E_K)/50.0;
r_infinity = 1.0/(1.0+exp((Y(25)-1.0)/-11.0));
tau_r = 0.0035*exp(-Y(25)*Y(25)/30.0/30.0)+0.0015;
dY(6, 1) = (r_infinity-Y(6))/(tau_r*1000.0);
s_infinity = 1.0/(1.0+exp((Y(25)+40.5)/11.5));
s_factor = (Y(25)+52.45)/15.8827;
tau_s = 0.025635*exp(-s_factor*s_factor)+0.01414;
dY(7, 1) = (s_infinity-Y(7))/(tau_s*1000.0);
f_Ca = Y(21)/(Y(21)+k_Ca);
i_Ca_L = g_Ca_L*Y(8)*(f_Ca*Y(9)+(1.0-f_Ca)*Y(10))*(Y(25)-E_Ca_app)/50.0;    %%%%%
d_L_infinity = 1.0/(1.0+exp((Y(25)+9.0)/-5.8));
d_L_factor = (Y(25)+35.0)/30.0;
tau_d_L = 0.0027*exp(-d_L_factor*d_L_factor)+0.002;
dY(8, 1) = (d_L_infinity-Y(8))/(tau_d_L*1000.0);
f_L_infinity = 1.0/(1.0+exp((Y(25)+27.4)/7.1));
f_L_factor = Y(25)+40.0;
tau_f_L1 = 0.161*exp(-f_L_factor*f_L_factor/14.4/14.4)+0.01;
dY(9, 1) = (f_L_infinity-Y(9))/(tau_f_L1*1000.0);
tau_f_L2 = 1.3323*exp(-f_L_factor*f_L_factor/14.2/14.2)+0.0626;
dY(10, 1) = (f_L_infinity-Y(10))/(tau_f_L2*1000.0);
i_NaCa = 1.0/50.0*K_NaCa*(Y(24)*Y(24)*Y(24)*Y(11)*exp(F*Y(25)*gamma_Na/(R*T))-Y(13)*Y(13)*Y(13)*Y(22)*exp((gamma_Na-1.0)*Y(25)*F/(R*T)))/(1.0+d_NaCa*(Y(13)*Y(13)*Y(13)*Y(22)+Y(24)*Y(24)*Y(24)*Y(11)));   %%
E_Na = R*T/F*log(Y(13)/Y(24));
i_B_Na = g_B_Na*(Y(25)-E_Na)/50.0;
E_Ca = R*T/(2.0*F)*log(Y(11)/Y(22));
i_B_Ca = g_B_Ca*(Y(25)-E_Ca)/50.0;
i_CaP = i_CaP_max*Y(22)/(50.0*(Y(22)+k_CaP));   %%%%
dY(11, 1) = ((Ca_b-Y(11))/tau_Ca+50.0*(i_Ca_L+i_B_Ca+i_CaP-2.0*i_NaCa)/(2.0*Vol_c*F))/1000.0;
i_Kur = g_kur*Y(29)*Y(30)*(Y(25)-E_K)/50.0;    %%%
i_K1 = g_K1*(Y(12)/1.0)^0.4457*(Y(25)-E_K)/(50.0*(1.0+exp(1.5*(Y(25)-E_K+3.6)*F/(R*T))));    %%%%
i_Ks = g_Ks*Y(14)*(Y(25)-E_K)/50.0;    %%%
%factorAZD = (1.0-1.0/(1.0+exp(-(AZD-0.6)/0.12)))/0.9933;
pip = 1.0/(1.0+exp((Y(25)+55.0)/24.0));
i_Kr = g_Kr*Y(15)*pip*(Y(25)-E_K)/50.0;   %%%
pow_Na_i_15 = Y(24)^1.5;
i_NaK = i_NaK_max*Y(12)/(Y(12)+K_NaK_K)*pow_Na_i_15/(pow_Na_i_15+pow_K_NaK_Na_15)*(Y(25)+150.0)/(50.0*(Y(25)+200.0));
dY(12, 1) = ((K_b-Y(12))/tau_K+50.0*(i_t+i_Kur+i_K1+i_Ks+i_Kr+i_KACh-2.0*i_NaK)/(Vol_c*F))/1000.0;
i_Na = P_Na*Y(28)*Y(28)*Y(28)*(0.9*Y(26)+0.1*Y(27))*Y(13)*Y(25)*F*F/(R*T)*(exp((Y(25)-E_Na)*F/(R*T))-1.0)/(50.0*(exp(Y(25)*F/(R*T))-1.0));    %%%
dY(13, 1) = ((Na_b-Y(13))/tau_Na+50.0*(i_Na+i_B_Na+3.0*i_NaCa+3.0*i_NaK+phi_Na_en)/(Vol_c*F))/1000.0;
n_infinity = 1.0/(1.0+exp((Y(25)-19.9)/-12.7));
n_factor = (Y(25)-20.0)/20.0;
tau_n = 0.7+0.4*exp(-n_factor*n_factor);
dY(14, 1) = (n_infinity-Y(14))/(tau_n*1000.0);
p_a_infinity = 1.0/(1.0+exp((Y(25)+15.0)/-6.0));
pa_factor = (Y(25)+20.1376)/22.1996;
tau_pa = 0.03118+0.21718*exp(-pa_factor*pa_factor);
dY(15, 1) = (p_a_infinity-Y(15))/(tau_pa*1000.0);
J_O_C = 200000.0*Y(22)*(1.0-Y(17))-476.0*Y(17);
J_O_TC = 78400.0*Y(22)*(1.0-Y(18))-392.0*Y(18);
J_O_TMgC = 200000.0*Y(22)*(1.0-Y(19)-Y(20))-6.6*Y(19);
J_O_TMgMg = 2000.0*Mg_i*(1.0-Y(19)-Y(20))-666.0*Y(20);
dY(17, 1) = J_O_C/1000.0;
dY(18, 1) = J_O_TC/1000.0;
dY(19, 1) = J_O_TMgC/1000.0;
dY(20, 1) = J_O_TMgMg/1000.0;
J_O = 0.08*J_O_TC+0.16*J_O_TMgC+0.045*J_O_C;
dY(16, 1) = J_O/1000.0;
past = floor(time/settings.BCL)*settings.BCL;

if ((time-past >= settings.StimOffset) && (time-past <= settings.StimOffset+settings.Dur_stim))
  i_Stim = settings.Amp_stim;
else
  i_Stim = 0.0;
end

dY(23, 1) = -50.0*(i_t+i_Kur+i_K1+i_Ks+i_Kr+i_KACh-2.0*i_NaK+i_Stim)/(Vol_i*F)/1000.0;
dY(24, 1) = -50.0*(i_Na+i_B_Na+3.0*i_NaCa+3.0*i_NaK+phi_Na_en)/(Vol_i*F)/1000.0;
i_di = (Y(21)-Y(22))*2.0*Vol_d*F/(tau_di*50.0);
dY(22, 1) = (-50.0*(i_B_Ca+i_CaP+i_up-(i_di+i_rel+2.0*i_NaCa))/(2.0*Vol_i*F)-1.0*J_O)/1000.0;
dY(21, 1) = -50.0*(i_Ca_L+i_di)/(2.0*Vol_d*F)/1000.0;
I = i_Na+i_Ca_L+i_t+i_Kur+i_K1+i_Kr+i_Ks+i_B_Na+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_KACh+i_Stim;
dY(25, 1) = -I;
%Q_tot = 0.05*Y(25);
h_infinity = 1.0/(1.0+exp((Y(25)+63.6)/5.3));
h_factor = 1.0/(1.0+exp((Y(25)+35.1)/3.2));
tau_h1 = 0.03*h_factor+0.0003;
dY(26, 1) = (h_infinity-Y(26))/(tau_h1*1000.0);
tau_h2 = 0.12*h_factor+0.003;
dY(27, 1) = (h_infinity-Y(27))/(tau_h2*1000.0);
m_infinity = 1.0/(1.0+exp((Y(25)+27.12)/-8.21));
m_factor = (Y(25)+25.57)/28.8;
tau_m = 0.000042*exp(-m_factor*m_factor)+0.000024;
dY(28, 1) = (m_infinity-Y(28))/(tau_m*1000.0);
a_ur_infinity = 1.0/(1.0+exp(-(Y(25)+6.0)/8.6));
tau_a_ur = 0.009/(1.0+exp((Y(25)+5.0)/12.0))+0.0005;
dY(29, 1) = (a_ur_infinity-Y(29))/(tau_a_ur*1000.0);
i_ur_infinity = 1.0/(1.0+exp((Y(25)+7.5)/10.0));
tau_i_ur = 0.59/(1.0+exp((Y(25)+60.0)/10.0))+3.05;
dY(30, 1) = (i_ur_infinity-Y(30))/(tau_i_ur*1000.0);


% added to see currents values
% dY(31,1) = i_K1;       % K1 current
% dY(32,1) = i_Ca_L;     % late Ca current
% dY(33,1) = i_rel;      % Ca current in sarcoplasmatic reticulum
% dY(34,1) = i_t;        % transient outward K current (Ca dependent)
% dY(35,1) = i_Kr;       % rapid delay rectifier K current
% dY(36,1) = i_Ks;       % slow delay rectifier K current
% dY(37,1) = i_Na;       % fast Na current
% dY(38,1) = i_Kur;      % ultra rapid K current


if flag==1
    dY = [i_Ca_L i_K1 i_rel i_t i_Kr i_Ks i_Na i_Kur i_CaP];
end


%===============================================================================
% End of file
%===============================================================================
end
