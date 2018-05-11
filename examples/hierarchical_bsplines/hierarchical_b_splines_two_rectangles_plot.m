%% hierarchical B-Splines mesh information, (c) Hoang Giang Bui, 2018
clc
clear
close all
hold on

hierarchical_b_splines_two_rectangles

param1.color = 'red';
param2.color = 'blue';

% plot_ctrl_points_hbsplines_2d(P1_P,P1_W,P1_Id,param1);
% plot_ctrl_points_hbsplines_2d(P2_P,P2_W,P2_Id,param2);

figure
plot_sampling_hbsplines_cells(P1_Xi,P1_Eta,P1_P,P1_W,P1_Id,P1_S,P1_C,P1_N,P1_params);
plot_sampling_hbsplines_cells(P2_Xi,P2_Eta,P2_P,P2_W,P2_Id,P2_S,P2_C,P2_N,P2_params);
