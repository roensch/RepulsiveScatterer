clear all;

parameter_filename = 'reconstruction_parameters';

addpath(genpath('./Forward_Problem'));
addpath('./');
addpath(genpath('./Iteration_Scheme'));
addpath(genpath('./submodules'));
addpath(genpath('./Meshes'));
close all;

%% read parameter file
[p,reg_par] = feval(parameter_filename, 'spot','spot_test'); % specify experiment and filename for reconstruction

%% create operator
F = feval(p.op_name,p);

%% create synthetic data
[data,F] = F.create_synthetic_data(0); % function parameter specifies GPU-device. Defaults to 0 if empty.

%% create regularization method
Reg = feval(reg_par.method,reg_par);
    
%% perform inversion
[result, statistics] = Reg.solve(data,F);
