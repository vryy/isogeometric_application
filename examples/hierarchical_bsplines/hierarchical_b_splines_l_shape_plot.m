%% hierarchical B-Splines mesh information, (c) Hoang Giang Bui, 2018
clc
clear
close all
hold on

hierarchical_b_splines_l_shape

%param1.color = 'red';
%param2.color = 'blue';

% plot_ctrl_points_hbsplines_2d(P1_P,P1_W,P1_Id,param1);
% plot_ctrl_points_hbsplines_2d(P2_P,P2_W,P2_Id,param2);

P1_params.tol = 1e-6;
P1_params.max_xi = 1;
P1_params.max_eta = 1;
P1_params.method = 'bezier';

P2_params.tol = 1e-6;
P2_params.max_xi = 1;
P2_params.max_eta = 1;
P2_params.method = 'bezier';

P3_params.tol = 1e-6;
P3_params.max_xi = 1;
P3_params.max_eta = 1;
P3_params.method = 'bezier';

plot_sampling_hbsplines_cells(P1_Xi,P1_Eta,P1_P,P1_W,P1_EqId,P1_S,P1_C,P1_N,P1_params);
plot_sampling_hbsplines_cells(P2_Xi,P2_Eta,P2_P,P2_W,P2_EqId,P2_S,P2_C,P2_N,P2_params);
plot_sampling_hbsplines_cells(P3_Xi,P3_Eta,P3_P,P3_W,P3_EqId,P3_S,P3_C,P3_N,P3_params);

