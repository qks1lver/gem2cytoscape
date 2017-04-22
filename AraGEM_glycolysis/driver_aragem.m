%% driver genCytoscapeNetwork
clear
clc

%% load model variables
load('20170414_ler_aragem_cytoscape.mat')

% create param
param.m = m;
param.pathways = pathways;
param.metFreq = metFreq;
param.rids = rids;
param.pName = 'Glycolysis / Gluconeogenesis';
param.fName = 'aragem_glycolysis';
param.maxMetFreq = 10;
param.excludeMets = {};
param.includeMets = [{'Pyruvate_c'} {'Acetate_c'} {'CO2_c'} {'alpha-D-Glucose_c'} {'beta-D-Fructose 6-phosphate'}];

% calculate line weights
load('20170415_all_sim_sols.mat')
d = 0.000000001;
timeframe = 1:timeParam.eodi-1;
v = ler_post;
v(abs(v) < d) = 0;
v = mean(v(:,timeframe),2);

%% generate network
param = genCytoscapeNetwork(param, v);

clearvars -except param
