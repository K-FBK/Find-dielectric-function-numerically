function [wl, n, k] = ImportRefractiveIndex(filename, dataLinesN, dataLinesK)
%IMPORTFILE Import data from a text file
%  [WL, N, K] = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as column vectors.
%
%  [WL, N, K] = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  [wl, n, k] = importfile(Ninomiya-o_CdSe_bulk.txt", [2, 501], [504, 1003]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 2023-01-26 09:26:44

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 3
    dataLinesN = [2, 501];
    dataLinesK = [504, 1003];
end

%% Set up the Import Options and import the N data
optsN = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
optsN.DataLines = dataLinesN;
optsN.Delimiter = "\t";

% Specify column names and types
optsN.VariableNames = ["wl", "n"];
optsN.VariableTypes = ["double", "double"];

% Specify file level properties
optsN.ExtraColumnsRule = "ignore";
optsN.EmptyLineRule = "read";

% Import the data
tblN = readtable(filename, optsN);

%% Set up the Import Options and import the K data
optsK = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
optsK.DataLines = dataLinesK;
optsK.Delimiter = "\t";

% Specify column names and types
optsK.VariableNames = ["wl", "k"];
optsK.VariableTypes = ["double", "double"];

% Specify file level properties
optsK.ExtraColumnsRule = "ignore";
optsK.EmptyLineRule = "read";

% Import the data
tblK = readtable(filename, optsK);

%% Convert to output type
wl = tblN.wl*1000; %Convert to nm
n = tblN.n;
k = tblK.k;
end