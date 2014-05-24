close all; clear all; clc;

num_trials = 600;
num_ticks = 3000;

% load sensory_output_0
% load msn_output_0
% load pf_output_0
% load tan_output_0
% load motor_output_0
load response_0
load response_time_0
load w_ctx_msn_0
load w_pf_tan_0
load reward_0
load predicted_reward_0
load rpe_0
load DA_0

% load sensory_output_1
% load msn_output_1
% load pf_output_1
% load tan_output_1
% load motor_output_1
load response_1
load response_time_1
load w_ctx_msn_1
load w_pf_tan_1
load reward_1
load predicted_reward_1
load rpe_1
load DA_1

% pf_output_0 = reshape(pf_output_0, num_ticks, num_trials);
% tan_output_0 = reshape(tan_output_0, num_ticks, num_trials);
% sensory_output_0 = reshape(sensory_output_0, num_ticks, num_trials);
% msn_output_0 = reshape(msn_output_0, num_ticks, num_trials);
% motor_output_0 = reshape(motor_output_0, num_ticks, num_trials);
% 
% pf_output_1 = reshape(pf_output_1, num_ticks, num_trials);
% tan_output_1 = reshape(tan_output_1, num_ticks, num_trials);
% sensory_output_1 = reshape(sensory_output_1, num_ticks, num_trials);
% msn_output_1 = reshape(msn_output_1, num_ticks, num_trials);
% motor_output_1 = reshape(motor_output_1, num_ticks, num_trials);

%% plot learning
color = lines(3);

figure
subplot(2,3,1), hold all
% plot(reward_0)
plot(predicted_reward_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(DA_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
legend({'predicted reward','DA release'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

subplot(2,3,2), hold all
plot(w_ctx_msn_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(w_pf_tan_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
legend({'ctx-msn weight','pf-tan weight'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

subplot(2,3,3), hold all
plot(response_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', [0 0 0])
axis([0 num_trials+1 0 1])

subplot(2,3,4), hold all
% plot(reward_1)
plot(predicted_reward_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(DA_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
legend({'predicted reward','DA release'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

subplot(2,3,5), hold all
plot(w_ctx_msn_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(w_pf_tan_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
legend({'ctx-msn weight','pf-tan weight'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

subplot(2,3,6), hold all
plot(response_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', [0 0 0])
axis([0 num_trials+1 0 1])

%% plot circuit
% trial = 1;
% figure
% 
% subplot(2,2,1), hold all
% plot(sensory_output_0(:,trial))
% plot(pf_output(:,trial))
% legend({'sensory output', 'pf output'})
% 
% subplot(2,2,2), hold all
% plot(tan_output_0(:,trial))
% legend({'tan output'})
% 
% subplot(2,2,3), hold all
% plot(msn_output_0(:,trial))
% legend({'msn output'})
% 
% subplot(2,2,4), hold all
% plot(motor_output_0(:,trial))
% legend({'motor output'})

%%

%% plot learning
color = lines(3);

figure, hold
plot(predicted_reward_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(DA_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
axis square
set(gca,'XTick',100:100:num_trials, 'fontsize', 10, 'fontweight', 'b')
xlabel('Trial', 'fontsize', 18, 'fontweight', 'b')
ylabel('AU', 'fontsize', 18, 'fontweight', 'b')
legend({'Predicted Reward','DA Release'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

figure, hold
plot(w_ctx_msn_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(w_pf_tan_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
axis square
set(gca,'XTick',100:100:num_trials, 'fontsize', 10, 'fontweight', 'b')
xlabel('Trial', 'fontsize', 18, 'fontweight', 'b')
ylabel('AU', 'fontsize', 18, 'fontweight', 'b')
legend({'CTX-MSN Weight','Pf-TAN Weight'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

figure, hold
plot(predicted_reward_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(DA_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
axis square
set(gca,'XTick',100:100:num_trials, 'fontsize', 10, 'fontweight', 'b')
xlabel('Trial', 'fontsize', 18, 'fontweight', 'b')
ylabel('AU', 'fontsize', 18, 'fontweight', 'b')
legend({'Predicted Reward','DA Release'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

figure, hold
plot(w_ctx_msn_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(w_pf_tan_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
axis square
set(gca,'XTick',100:100:num_trials, 'fontsize', 10, 'fontweight', 'b')
xlabel('Trial', 'fontsize', 18, 'fontweight', 'b')
ylabel('AU', 'fontsize', 18, 'fontweight', 'b')
legend({'CTX-MSN Weight','Pf-TAN Weight'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff

figure, hold
plot(response_0, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,1))
plot(response_1, 'linewidth', 2, 'linesmoothing', 'on', 'color', color(:,3))
axis([0 num_trials+1 0 1])
axis square
set(gca,'XTick',100:100:num_trials, 'fontsize', 10, 'fontweight', 'b')
xlabel('Trial', 'fontsize', 18, 'fontweight', 'b')
ylabel('Response Probability', 'fontsize', 18, 'fontweight', 'b')
legend({'CRF','PRF'}, 'location', 'NorthEast', 'fontsize', 18)
legend boxoff