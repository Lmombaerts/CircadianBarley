% This script presents the code used to infer the circadian regulatory
% network of barley. 
% Code Author: Laurent Mombaerts @Octobre 2018

%% Network reconstruction Section: Performance analysis on simulated data
clear; close all;

% Load simulated data (Pokhilko et al. 2010)
load millar10.mat 
mRNA_idx = [1 4 7 10 12 14 16]; 
mRNA_names = {'LHY mRNA';'TOC1 mRNA';'Y mRNA';'PRR9 mRNA';'PRR7 mRNA';'NI mRNA';'GI mRNA'}; 
samplingTime = 4;

for simuNumber = 1:50
    mRNA_data = LL{1,simuNumber}(mRNA_idx,1:samplingTime:end);

    %% Apply 1st Order Linear Model ID
    % Resample to realistic data (48h of transients data, every 4 hours)
    tL = 0:samplingTime:(size(mRNA_data,2)-1)*samplingTime; 
    
    % Interpolation to 1h (piecewise cubic spline)
    tLI = 0:48;
    dataLI = pcs(mRNA_data, tL, tLI); 
    
    %Network Topology Inference (as a fitness Matrix)
    [fitness, models, ~, ~, ~] = just_tfest(1, 1, dataLI);
    [TP{simuNumber},FP{simuNumber},TN{simuNumber},FN{simuNumber},TPR{simuNumber},FPR{simuNumber},SPC{simuNumber},PPV{simuNumber},AUROC_DT{simuNumber},AUPREC_DT{simuNumber}] = ROC_Millar10(fitness); 
    fitness_ATA{simuNumber} = fitness;
end
save results

%% Model Comparison Section. 
%The following function is intended to be used when comparing two models
%generated from the "just_tfest" function under different experimental
%conditions. For circadian data, the frequency range is adapted between 16
%hours period and 32 hours. Change this to fit your data.

model_example_untreated = models{1,2};
model_example_treated = models{2,3};
[nugaps,freqs,distnu,w] = richgap(model_example_untreated, model_example_treated, 2*pi/16, 2*pi/32); 

