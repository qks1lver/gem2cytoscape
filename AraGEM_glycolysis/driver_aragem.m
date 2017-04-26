%% driver genCytoscapeNetwork
clear
clc

%% load model variables
load('aragem_workspace.mat')

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

%% generate network
param = genCytoscapeNetwork(param, flux);

clearvars -except param
