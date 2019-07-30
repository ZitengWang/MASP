clear all; clc;

addpath audio

name = 'test_filenames.txt';
Names= textread(name,'%q');

no_files = size(Names,1);
SRMRstar = zeros(no_files,1);
ModSpecEnergy = cell(no_files,1);

% base_values = xlsread('Test_results.xls');


for j=1:no_files
    clc; disp(['Processing file ', num2str(j), ' of ', num2str(no_files)]);

   SRMRstar(j) = SRMR_main(char(Names(j)));
%    assert(abs(SRMRstar(j) - base_values(j,2)) < 1e100);
end