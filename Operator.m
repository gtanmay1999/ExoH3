% clc;
% clear;
% close all;

%% Exo-H3
userWeight = 160; % lbs
userHeight = 68; % in

amplitude = 0.8; % Arbitrary value of step length
frequency = 1.2; % Step frequency
hipOffset = -17; % degrees
ramp = 2; % seconds

%% PARAMETERS
maxTorque = 35; % N*m *** Likely won't need
pressureDeadZone = 0.2*userWeight;
strategySaturation = 10; % Degrees
P_t_ZMP_gain_Ankle = 0; % Ankle strategy gain (0.1 - 0.9 for dependent, 1 for independent)
P_t_ZMP_gain_Hip = 0; % Hip strategy gain (0.1 - 0.9 for dependent, 1 for independent)

T = 2*pi/frequency;
eqn = @(gamma) P_t_ZMP_gain_Ankle/gamma*(exp(gamma*T)-1) - T;
gamma = fzero(eqn, 0.1)

Kp = 10; % Only for torque controller
KD = 0.39; % Only for torque controller
ankleStrategyInterval = 1.5; % 0 to 2
dt = 0.01;
K_high = 0.5; %1 % Gain when 1 foot is on ground
K_low = 0.5; %0.5 % Gain when 2 feet are on ground

% Lowpass filter
fc = 3*frequency/pi*1.1; % [Hz]
RC = 1/(2*pi*fc);
DCM_filter = 1 - dt/(dt+RC);
Strategy_filter = DCM_filter;

% Phase offsets
Phase_Right=2;
Phase_Left=2+pi;

% Hip trajectory for walking
A_0_h = 40.69; A_1_h = 23.22; A_2_h = -4.487; A_3_h = 0.3995;
A_4_h = 0.7047; A_5_h = 1.078; A_6_h = -0.2732;
B_1_h = -8.65; B_2_h = 3.338; B_3_h = 1.389; B_4_h = 0.7989;
B_5_h = 0.342; B_6_h = 0.06696;

% Knee trajectory for walking
A_0_k = 25.7; A_1_k = -3.827; A_2_k = -8.543; A_3_k = 1.913;
A_4_k = 1.087; A_5_k = 2.054; A_6_k = -0.3061;
B_1_k = -19.28; B_2_k = 17.93; B_3_k = 3.772; B_4_k = 1.504;
B_5_k = 0.5838; B_6_k = -0.8958;

% Ankle trajectory for walking
A_0_a = 4.011; A_1_a = 5.292; A_2_a = -8.582; A_3_a = -0.4755;
A_4_a = 1.688; A_5_a = -0.03504; A_6_a = 1.298;
B_1_a = 3.613; B_2_a = -5.954; B_3_a = 5.918; B_4_a = -2.118;
B_5_a = 1.023; B_6_a = -0.7245;

% Segment mass [lbs] (average male)
m1 = userWeight*(0.1072+0.0459+0.0169);
m2 = userWeight*(0.0783+0.4712+2*0.0299+2*0.0176);
m3 = userWeight*0.10;
m4 = userWeight*(0.0459+0.0169);

% Including Exo-H3 weights
m1 = m1 + 2.2048*(0.909+1.647+1.715);
m2 = m2 + 2.2048*7.591;
m3 = m3 + 2.2048*1.715;
m_Exo = 2.2048*1.647;
m_f = userWeight*0.0169;

% Segment length [in] (average male)
l1 = userHeight*(0.232+0.247);
l2 = userHeight*(0.19+0.52)/(0.19+0.52+0.40+0.43);
l3 = userHeight*0.232;
l4 = userHeight*0.247;
l_foot = 11.25;

% Segment center of mass location measured from proximal end [in]
% (average male) (assuming consistent with Exo-H3)
l_c1 = ((0.17*7.77+0.59*3.52+(0.40+0.43)*1.05)/(7.77+3.52+1.05))/(0.40+0.43)*l1;
l_c2 = ((0.29*34.66+0.63*6.07)/(6.07+34.66))/(0.52+0.19)*l2;
l_c3 = 0.17/0.40*l3;
l_c4 = ((0.19*3.52+0.43*1.05)/(3.52+1.05))/0.43*l4;


% % REINFORCEMENT LEARNING PARAMETERS
% 
% Reward Scaling Parameters (TUNE)
% DCMErrorPenalty    = 1;
% totalPowerPenalty  = 1;
% 
% LLM = 0.9; % Lower Limit Multiplier (Assumed)
% ULM = 1.1; % Upper Limit Multiplier (Assumed)
% 
% obsInfo = rlNumericSpec([52 1], ...
%     LowerLimit = [0, 0, 0, 0, ...
%         LLM*m1, LLM*m2, LLM*m3, LLM*m4, ...
%         LLM*l_c1, LLM*l_c2, LLM*l_c3, LLM*l_c4, ...
%         LLM*I1, LLM*I2, LLM*I3, LLM*I4, ...
%         -25, 0, -25, -25, 0, -25, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, ...
%         -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -35, -35, -35, -35, -35, -35, ...
%         0, 0, 0, 0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, 0]', ...
%     UpperLimit = [Inf, Inf, Inf, Inf, ...
%         ULM*m1, ULM*m2, ULM*m3, ULM*m4, ...
%         ULM*l_c1, ULM*l_c2, ULM*l_c3, ULM*l_c4, ...
%         ULM*I1, ULM*I2, ULM*I3, ULM*I4, ...
%         100, 100, 25, 100, 100, 25, Inf, Inf, Inf, Inf, Inf, Inf, ...
%         Inf, Inf, Inf, Inf, Inf, Inf, 35, 35, 35, 35, 35, 35, ...
%         Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf]');
% obsInfo.Name = "Observations";
% obsInfo.Description = "P_tZMP_gain, Kp, KD, Ankle Strategy Interval" + ...
%     "W_1, W_2, W_3, W_4, CoM_1, CoM_2, CoM_3, CoM_4, MOI_1, MOI_2, MOI_3, MOI_4" + ...
%     "Desired RHA, Desired RKA, Desired RAA, Desired LHA, Desired LKA, Desired LAA" + ...
%     "Actual RHA, Actual RKA, Actual RAA, Actual LHA, Actual LKA, Actual LAA" + ...
%     "RHT, RKT, RAT, LHT, LKT, LAT, RHMT, RKMT, RAMT, LHMT, LKMT, LAMT" + ...
%     "RHP, RTP, LHP, LTP, RHAV, RKAV, RAAV, LHAV, LKAV, LAAV" + ...
%     "DCM Error, Total Power";
% 
% Reinforcement Learning
% P_t_ZMP_gain_Limit = [0, 0.01]; % Upper limit assumed (Modify w/o RL first)
% gamma_Limit = [0.01, 1]; % Lower limit assumed (must be >0)
% 
% actInfo = rlNumericSpec([2 1], ...
%     'LowerLimit', [P_t_ZMP_gain_Limit(1); gamma_Limit(1)], ...
%     'UpperLimit', [P_t_ZMP_gain_Limit(2); gamma_Limit(2)]);
% actInfo.Name = "Actions";
% actInfo.Description = "P_tZMP_gain, gamma";
% 
% env = rlSimulinkEnv('ankleStrategy', 'ankleStrategy/RL Agent', obsInfo, actInfo);
% 
% initOpts = rlAgentInitializationOptions('NumHiddenUnit', 64);
% 
% criticOpts = rlOptimizerOptions('LearnRate',1e-03,'GradientThreshold',1);
% actorOpts = rlOptimizerOptions('LearnRate',1e-04,'GradientThreshold',1);
% 
% agentOpts = rlTD3AgentOptions(...
%     'SampleTime',-1,...
%     'TargetSmoothFactor',1e-3,...
%     'DiscountFactor',0.99,...
%     'MiniBatchSize',64,...
%     'ExperienceBufferLength',1e6,...
%     'TargetUpdateFrequency',4, ...
%     'CriticOptimizerOptions',criticOpts, ...
%     'ActorOptimizerOptions',actorOpts);
%     agentOpts.ExplorationModel.StandardDeviation = [0.1; 0.1];
% 
% agent = rlTD3Agent(obsInfo, actInfo, agentOpts);
% 
