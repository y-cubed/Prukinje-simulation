%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Katri Hituri 2013
% katri.hituri@tut.fi / katri.hituri@gmail.com
%
% A script to simulate IP3R models and 
% 
% Hituri K, Linne M-L. Comparison of Models for IP3 
% Receptor Kinetics using Stochastic Simulations.
% PLOS ONE. 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script is cell mode based. You can run each cell individually by 
% pressing the 'Evaluate cell' button or the whole script by typing 
% run_IP3R_P0 on MATLAB command line. 

%% Run the simulations
% You can run the simulations here or on your computers command line. 
% Please make sure that you have STEPS and NumPy installed.
% WARNING! These simulations will take several hours in total. 

!python models/sim_doi.py
!python models/sim_fd.py
!python models/sim_dli.py
!python models/sim_ot.py 


%% Import the simulations results

op_doi = ... 
    importdata('ip3r_doi_op_res.dat');
op_doi_ca_concs = ... 
    importdata('ip3r_doi_op_ca_concs.dat');

op_fd = ... 
    importdata('ip3r_fd_op_res.dat');

op_fd_ca_concs = ... 
    importdata('ip3r_fd_op_ca_concs.dat');

op_ot = ... 
    importdata('ip3r_ot_op_res.dat');

op_ot_ca_concs = ... 
    importdata('ip3r_ot_op_ca_concs.dat');

op_dli = ... 
    importdata('ip3r_dli_op_res.dat');

op_dli_ca_concs = ... 
    importdata('ip3r_dli_op_ca_concs.dat');


%% Plot the figure 2A 

figure();

% Ca concentrations are multiplied by 1e6 to have them in microM
fig = semilogx(op_ot_ca_concs*1e6, op_ot(:,1), 'LineWidth',3, 'Color', 'g');
hold on;
semilogx(op_dli_ca_concs*1e6, op_dli(:,1), 'LineWidth',3, 'Color', 'b');
semilogx(op_fd_ca_concs*1e6, op_fd(:,1), 'LineWidth',3, 'Color', 'r');
semilogx(op_doi_ca_concs*1e6, op_doi(:,1), 'LineWidth',3, 'Color', 'm');

axis([1e-3,1e2,0,0.8]);
xlabel('[Ca^{2+}] \muM', 'fontsize', 16);
ylabel('P_o', 'fontsize', 16);
set(findobj('type','axes'),'fontsize', 16);
title('Open probability of IP_3R as a function of cytosolic [Ca^{2+}] in different models');

legend('Othmer & Tang', 'Dawson, Lea & Irvine','Fraiman & Dawson', 'Doi et al.', 'Location', 'Northwest');

hold off;

%% Save the figure

saveas(fig, 'P0_Ca.fig');
saveas(fig, 'P0_Ca.jpg');
