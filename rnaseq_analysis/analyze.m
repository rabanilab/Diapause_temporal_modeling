addpath('../MATLAB');
addpath('../rnaseq_analysis');

% -----------------------------------------------------------------------
% TSO order
% -----------------------------------------------------------------------

% test traveling salesmen solving
analyze_solve;

% run tso on RNA-seq data
analyze_tso;

% -----------------------------------------------------------------------
% fit sigmoid models
% -----------------------------------------------------------------------
analyze_model;
