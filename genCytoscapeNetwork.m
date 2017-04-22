function param = genCytoscapeNetwork(param, v)

%% generate a .sif file for Cytoscape for pathway name(s)
%
% Author: Jiun Yen
% Date: 2017.4.14
% Version: 2017.4.22
%
% param structure
%   .m COBRA-format GEM model
%   .pathways structure
%       .rids
%       .names
%       .values
%       .counts
%   .metFreq frequency of metabolite used in reactions
%   .rids indices of reactions with fluxes
%   .pName path name
%   .fName file name
%   .maxMetFreq upper limit for metabolite frequency in S matrix
%   .excludeMets metabolites excluded from network
%   .includeMets metabolites to include if they are in the pathways

%% find pathway rids
p_rids = param.pathways.rids(ismember(param.pathways.values,find(ismember(param.pathways.names, param.pName))));

%% remove inactive rids
try
    rids = intersect(param.rids, p_rids);
catch
    error('No rids');
end

%% Update param
S = full(param.m.S);
param.networkRids = rids;
param.networkRxns = cell(length(rids), 1);
param.networkInteractions = cell(2*length(rids), 1);
allmets = logical(sum(abs(S(:,rids)),2));
excludemets = ismember(param.m.metNames, param.excludeMets);
includemets = ismember(param.m.metNames, param.includeMets);
hifreqmets = param.metFreq > param.maxMetFreq;
excludemets = (excludemets | (allmets & hifreqmets)) & (allmets & ~includemets);
param.allMets = param.m.metNames(allmets);
param.excludeMets = param.m.metNames(excludemets);
param.networkMets = param.m.metNames(allmets & ~excludemets);
param.rids = rids;

%% construct network

% remove mets that are exclude from S-matrix
S(excludemets, :) = [];
param.m.metNames(excludemets) = [];

% LHS S-matrix (reactants)
SL = S < 0;

% RHS S-matrix (products)
SR = S > 0;

% write to file
f = fopen([param.fName '_network.sif'], 'w');
for i = 1:length(rids)
    rxnName = editName(param.m.rxnNames{rids(i)});
    metL = param.m.metNames(SL(:,rids(i)));
    metR = param.m.metNames(SR(:,rids(i)));
    tmp = sprintf('r%u',rids(i));
    param.networkRxns(i) = {[tmp ':' rxnName]};
    param.networkInteractions(2*i-1) = {[tmp 'L']};
    param.networkInteractions(2*i) = {[tmp 'R']};
    if ~isempty(metL) && ~isempty(metR)
        for a = 1:length(metL)
            fprintf(f, [param.networkRxns{i} '\t' param.networkInteractions{2*i-1} '\t' metL{a} '\n']);
        end
        for b = 1:length(metR)
            fprintf(f, [param.networkRxns{i} '\t' param.networkInteractions{2*i} '\t' metR{b} '\n']);
        end
    end
end
fclose(f);

%% construct node properties
f = fopen([param.fName '_node.txt'], 'w');
fprintf(f, 'name\ttype\n');
for i = 1:length(param.networkMets)
    fprintf(f, [param.networkMets{i} '\tmet\n']);
end
for i = 1:length(param.networkRxns)
    fprintf(f, [param.networkRxns{i} '\trxn\n']);
end
fclose(f);

%% construct edge properties

% check if necessary
if nargin < 2 || isempty(v)
    return;
end

% build directionality vector
dL = (v(rids) < 0) + 0;
dR = (v(rids) > 0) + 0;
v = abs(v(rids));
toDraw = (v ~= 0) + 0;

% build edge properties
f = fopen([param.fName '_edge.txt'], 'w');
fprintf(f, 'interaction\tdirection\tdraw\tweight\n');
for i = 1:length(rids)
    fprintf(f, [param.networkInteractions{2*i-1} '\t%u\t%u\t%4.4f\n'], [dL(i) toDraw(i) v(i)]);
    fprintf(f, [param.networkInteractions{2*i} '\t%u\t%u\t%4.4f\n'], [dR(i) toDraw(i) v(i)]);
end
fclose(f);

function name = editName(name)
%% trim rxnName to critical component

tmp = strsplit(name, ' / ');
name = strrep(tmp{1}, ', putative', '');
