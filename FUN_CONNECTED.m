fprintf('\n=================================================================')
fprintf('\nDESCRIPTION: Finds connected sets -- one regular connected set and one leave-one-out connected set -- of workers and employers.')
fprintf('\n')
fprintf('\nINPUTS:      Nx3 matrix of worker IDs, firm IDs, and year IDs.')
fprintf('\n')
fprintf('\nOUTPUTS:     - M_{1}x1 matrix of worker IDs in regular connected set for firm FE model, where M_{1} <= N.')
fprintf('\n             - M_{2}x1 matrix of worker IDs in leave-one-out connected set for firm FE model, where M_{2} <= M_{1} <= N.')
fprintf('\n             - M_{3}x1 matrix of worker IDs in regular connected set for firm-year FE model, where M_{3} <= M_{1} <= N.')
fprintf('\n             - M_{4}x1 matrix of worker IDs in leave-one-out connected set for firm-year FE model, where M_{4} <= M_{3} <= M_{1} <= N.')
fprintf('\n')
fprintf('\nCREDITS:     Please cite the following paper:')
fprintf('\n             - Engbom, Niklas & Christian Moser & Jan Sauermann, 2021. "Firm Pay Dynamics," Working Paper.')
fprintf('\n')
fprintf('\nTIME STAMP:  October 27, 2021.')
fprintf('\n=================================================================')
fprintf('\n')


fprintf('\n==============================================================')
fprintf('\n OPENING HOUSEKEEPING')
fprintf('\n==============================================================')
fprintf('\n  --> start timer\n')
tic

fprintf('\n  --> clear memory\n')
clc
clear

fprintf('\n  --> stop and go into debugging model when encountering error\n')
dbstop if error

fprintf('\n  --> verify MATLAB version\n')
v = version;
v_old = (str2double(v(1:3)) < 9.1); % If MATLAB is older than version R2016b
if v_old
    fprintf('\nUSER WARNING: Implicit expansions require MATLAB R2016b or later version. Run slower routine instead.\n')
end
clear v v_old


fprintf('\n==============================================================')
fprintf('\n SET DIRECTORIES AND PARAMETERS')
fprintf('\n==============================================================')
fprintf('\n  --> set paths and macros not passed from Stata\n')
try_path_1 = 'P:/2019/107'; % IFAU/Uppsala server
try_path_2 = '/Users/cm3594/Dropbox (CBS)/AKM Conference'; % Chris' work iMac Pro
try_path_3 = '/Users/economoser/Dropbox (CBS)/AKM Conference/5_results/temp'; % Chris' private iMac Pro
try_path_4 = '/Users/niklasengbom/Dropbox/AKM'; % Nik's work computer
try_path_5 = '/shared/share_cmoser/20_AKM_Conference'; % Columbia GSB server
if exist(try_path_1, 'dir') % IFAU/Uppsala server
    DIR_MAIN = try_path_1;
elseif exist(try_path_2, 'dir') % Chris' work iMac Pro
    DIR_MAIN = try_path_2;
elseif exist(try_path_3, 'dir') % Chris' private iMac Pro
    DIR_MAIN = try_path_3;
elseif exist(try_path_4, 'dir') % Nik's work computer
    DIR_MAIN = try_path_4;
elseif exist(try_path_5, 'dir') % Columbia GSB server
    DIR_MAIN = try_path_5;
else
    fprintf('\nUSER ERROR: Directory not found.\n')
    exit
end
if strcmp(DIR_MAIN, try_path_1) % IFAU/Uppsala server
    DIR_INPUT = [DIR_MAIN, '/2 Data'];
    DIR_OUTPUT = [DIR_MAIN, '/2 Data'];
    DIR_LOG = [DIR_MAIN, '/3 Logs/'];
    DIR_LEAVEOUTTWOWAY = [DIR_MAIN, '/1 Code/_LeaveOutTwoWay_3_02_Chris'];
else % Chris' work iMac Pro or Chris' private iMac Pro or Nik's work computer or Columbia GSB server
    DIR_INPUT = [DIR_MAIN, '/5_results/temp'];
    DIR_OUTPUT = [DIR_MAIN, '/5_results/temp'];
    DIR_LOG = [DIR_MAIN, '/5_results/log'];
    DIR_LEAVEOUTTWOWAY = [DIR_MAIN, '/3_code/_work_in_progress/_LeaveOutTwoWay_3_02_Chris'];
end
cd(DIR_LEAVEOUTTWOWAY)
addpath(genpath(DIR_LEAVEOUTTWOWAY))
clear try_path_1 try_path_2 try_path_3 try_path_4 try_path_5 DIR_MAIN DIR_LEAVEOUTTWOWAY

fprintf('\n  --> read parameters\n')
params_dir = [DIR_INPUT '/parameters_connected.csv']; % Path to file containing AKM parameters.
file_params = fopen(params_dir);
data_params = fscanf(file_params, '%f %f %f %f', [4 inf])';
fclose(file_params);
year_min = data_params(:, 1); % Minimum year.
fprintf(['\nyear_min = ' num2str(year_min), '\n'])
year_max = data_params(:, 2); % Maximum year.
fprintf(['\nyear_max = ' num2str(year_max), '\n'])
inc_concept = data_params(:, 3); % Income threshold.
fprintf(['\ninc_concept = ' num2str(inc_concept), '\n'])
thresh = data_params(:, 4); % Minimum firm size threshold.
fprintf(['\nthresh = ' num2str(thresh), '\n'])
if inc_concept == 1
    inc_concept_str = 'earn';
elseif inc_concept == 2
    inc_concept_str = 'wage';
end
ext = ['_', num2str(year_min), '_', num2str(year_max)]; % ['_', inc_concept_str, '_', num2str(thresh), '_', num2str(year_min), '_', num2str(year_max)]
clear params_dir file_params data_params inc_concept inc_concept_str thresh

fprintf('\n  --> set read and write paths\n')
input_dir = [DIR_INPUT '/tomatlab_connected.csv']; % Path to file with Stata input to find connected sets: person_ID, employer_ID.
output_dir_connected_f_regular = [DIR_OUTPUT '/connected_f_regular.txt']; % Path to file where regular connected for firm FE model set is exported to.
output_dir_connected_f_leave_one_out = [DIR_OUTPUT '/connected_f_leave_one_out.txt']; % Path to file where leave-one-out connected set for firm FE model is exported to.
output_dir_connected_fy_regular = [DIR_OUTPUT '/connected_fy_regular.txt']; % Path to file where regular connected for firm-year FE model set is exported to.
output_dir_connected_fy_leave_one_out = [DIR_OUTPUT '/connected_fy_leave_one_out.txt']; % Path to file where leave-one-out connected set for firm-year FE model is exported to.
done_dir_connected = [DIR_OUTPUT '/done_connected.txt']; % Path to file that confirms being done with execution execution of FUN_CONNECTED.m.
success_dir_connected = [DIR_OUTPUT '/success_connected.txt']; % Path to file that confirms successful execution of FUN_CONNECTED.m.
clear DIR_INPUT

fprintf('\n  --> start diary (log) file\n')
FILE_LOG = [DIR_LOG, '/log_MATLAB_CONNECTED', ext, '.log'];
eval(['delete ''', FILE_LOG, ''''])
echo off all
eval(['diary ''', FILE_LOG, ''''])
fprintf('\n==============================================================')
fprintf('\n START LOG FILE FOR FUN_CONNECTED.m')
fprintf('\n==============================================================')
clear DIR_LOG FILE_LOG ext


fprintf('\n==============================================================')
fprintf('\n LOAD DATA')
fprintf('\n==============================================================')
fprintf('\n  --> prepare to read input data\n')
input_file = fopen(input_dir, 'r');
input_format = '%f %f %u'; % always load exactly 3 variables: id_pers_original, id_firm_original, year_akm
input_n = 3;
clear input_dir

fprintf('\n  --> read input data\n')
input_data = fscanf(input_file, input_format, [input_n inf])';
fclose(input_file);
id_pers_original = input_data(:, 1); % Input data column 1: worker ID.
id_firm_original = input_data(:, 2); % Input data column 2: employer ID.
year_akm = input_data(:, 3); % Input data column 3: year.
clear input_file input_format input_n input_data

fprintf('\n  --> generate firm-year ID\n')
id_firm_year = findgroups(id_firm_original, year_akm);
clear year_akm


fprintf('\n==============================================================')
fprintf('\n FIND CONNECTED SETS')
fprintf('\n==============================================================')
fprintf('\n  --> find regular connected set for firm FE model\n')
[id_pers_temp,id_firm_temp,id_pers_connected_f_regular,id_firm_connected_f_regular] = connected_set_Chris(id_pers_original,id_firm_original);
results = struct;
[~,~,id_pers_temp_unique] = unique(id_pers_temp);
results.N_pers_connected_regular_f = max(id_pers_temp_unique);
[~,~,id_firm_temp_unique] = unique(id_firm_temp);
results.N_firm_connected_regular_f = max(id_firm_temp_unique);
clear id_firm_original id_pers_temp_unique id_firm_temp_unique

fprintf('\n  --> find leave-one-out connected set for firm FE model\n')
[id_pers_temp,id_firm_temp,id_pers_connected_f_leave_one_out,~] = pruning_unbal_v3_Chris(id_pers_temp,id_firm_temp,id_pers_connected_f_regular,id_firm_connected_f_regular); % last argument not stored: id_firm_connected_f_leave_one_out
[~,~,id_pers_temp_unique] = unique(id_pers_temp);
results.N_pers_connected_leave_one_out_f = max(id_pers_temp_unique);
[~,~,id_firm_temp_unique] = unique(id_firm_temp);
results.N_firm_connected_leave_one_out_f = max(id_firm_temp_unique);
clear id_firm_temp id_firm_connected_f_regular id_pers_temp_unique id_firm_temp_unique

fprintf('\n  --> find regular connected set for firm-year FE model\n')
[id_pers_temp,id_firm_year_temp,id_pers_connected_fy_regular,id_firm_year_connected_fy_regular] = connected_set_Chris(id_pers_original,id_firm_year);
[~,~,id_pers_temp_unique] = unique(id_pers_temp);
results.N_pers_connected_regular_fy = max(id_pers_temp_unique);
[~,~,id_firm_year_temp_unique] = unique(id_firm_year_temp);
results.N_firm_year_connected_regular_fy = max(id_firm_year_temp_unique);
clear id_pers_original id_firm_year id_pers_temp_unique id_firm_year_temp_unique

fprintf('\n  --> find leave-one-out connected set for firm-year FE model\n')
[id_pers_temp,id_firm_year_temp,id_pers_connected_fy_leave_one_out,~] = pruning_unbal_v3_Chris(id_pers_temp,id_firm_year_temp,id_pers_connected_fy_regular,id_firm_year_connected_fy_regular); % last argument not stored: id_firm_connected_fy_leave_one_out
[~,~,id_pers_temp_unique] = unique(id_pers_temp);
results.N_pers_connected_leave_one_out_fy = max(id_pers_temp_unique);
[~,~,id_firm_year_temp_unique] = unique(id_firm_year_temp);
results.N_firm_year_connected_leave_one_out_fy = max(id_firm_year_temp_unique);
clear id_pers_temp id_firm_year_temp id_firm_year_connected_fy_regular id_pers_temp_unique id_firm_year_temp_unique

fprintf('\n  --> write summary of connected set to output file in MATLAB format\n')
save([DIR_OUTPUT, '/results_MATLAB_CONNECTED_', num2str(year_min), '_', num2str(year_max), '.mat'], 'results')
clear DIR_OUTPUT results year_min year_max


fprintf('\n==============================================================')
fprintf('\n EXPORT RESULTS')
fprintf('\n==============================================================')
fprintf('\n  --> export regular connected set for firm FE model to tab-delimited output file\n')
output_file_connected_f_regular = fopen(output_dir_connected_f_regular, 'w');
out_header_format = '%10s';
out_header = {'id_pers'};
out_data_format = '%10.0f';
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(output_file_connected_f_regular, out_header_format, out_header{:});
fprintf(output_file_connected_f_regular, out_data_format, unique(id_pers_connected_f_regular'));
fclose(output_file_connected_f_regular);
clear out_header_format out_header out_data_format id_pers_connected_f_regular output_file_connected_f_regular output_dir_connected_f_regular

fprintf('\n  --> export leave-one-out connected set for firm FE model to tab-delimited output file\n')
output_file_connected_f_leave_one_out = fopen(output_dir_connected_f_leave_one_out, 'w');
out_header_format = '%10s';
out_header = {'id_pers'};
out_data_format = '%10.0f';
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(output_file_connected_f_leave_one_out, out_header_format, out_header{:});
fprintf(output_file_connected_f_leave_one_out, out_data_format, unique(id_pers_connected_f_leave_one_out'));
fclose(output_file_connected_f_leave_one_out);
clear out_header_format out_header out_data_format id_pers_connected_f_leave_one_out output_file_connected_f_leave_one_out output_dir_connected_f_leave_one_out

fprintf('\n  --> export regular connected set for firm-year FE model to tab-delimited output file\n')
output_file_connected_fy_regular = fopen(output_dir_connected_fy_regular, 'w');
out_header_format = '%10s';
out_header = {'id_pers'};
out_data_format = '%10.0f';
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(output_file_connected_fy_regular, out_header_format, out_header{:});
fprintf(output_file_connected_fy_regular, out_data_format, unique(id_pers_connected_fy_regular'));
fclose(output_file_connected_fy_regular);
clear out_header_format out_header out_data_format id_pers_connected_fy_regular output_file_connected_fy_regular output_dir_connected_fy_regular

fprintf('\n  --> export leave-one-out connected set for firm-year FE model to tab-delimited output file\n')
output_file_connected_fy_leave_one_out = fopen(output_dir_connected_fy_leave_one_out, 'w');
out_header_format = '%10s';
out_header = {'id_pers'};
out_data_format = '%10.0f';
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(output_file_connected_fy_leave_one_out, out_header_format, out_header{:});
fprintf(output_file_connected_fy_leave_one_out, out_data_format, unique(id_pers_connected_fy_leave_one_out'));
fclose(output_file_connected_fy_leave_one_out);
clear out_header_format out_header out_data_format output_file_connected_fy_leave_one_out output_dir_connected_fy_leave_one_out

fprintf('\n  --> confirm being done with execution of FUN_CONNECTED.m\n')
done_file_connected = fopen(done_dir_connected, 'w');
out_header_format = '%10s';
out_header = {'done'};
out_data_format = '%10.0f';
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(done_file_connected, out_header_format, out_header{:});
fprintf(done_file_connected, out_data_format, 1);
fclose(done_file_connected);
clear out_header_format out_header out_data_format done_file_connected done_dir_connected

fprintf('\n  --> confirm successful execution of FUN_CONNECTED.m\n')
if exist('id_pers_connected_fy_leave_one_out','var')
    success_file_connected = fopen(success_dir_connected, 'w');
    out_header_format = '%10s';
    out_header = {'success'};
    out_data_format = '%10.0f';
    out_header_format = strcat(out_header_format, '\n');
    out_data_format = strcat(out_data_format, '\n');
    fprintf(success_file_connected, out_header_format, out_header{:});
    fprintf(success_file_connected, out_data_format, id_pers_connected_fy_leave_one_out(1));
    fclose(success_file_connected);
    clear id_pers_connected_fy_leave_one_out success_file_connected out_header_format out_header out_data_format success_dir_connected
end


fprintf('\n==============================================================')
fprintf('\n CLOSING HOUSEKEEPING')
fprintf('\n==============================================================')
fprintf('\n  --> summarize objects stored in memory\n')
whos

fprintf('\n  --> clear memory\n')
clear

fprintf('\n  --> end timer\n')
toc

fprintf('\n  --> close diary (log) file\n')
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n')
diary off

fprintf('\n  --> exit\n')
exit