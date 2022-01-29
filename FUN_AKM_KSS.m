fprintf('\n==============================================================')
fprintf('\nDESCRIPTION: FUN_AKM_KSS.m estimates KSS correction to AKM')
fprintf('\n             variance decomposition.')
fprintf('\n')
fprintf('\nINPUTS:      A file in .csv format prepared in Stata:')
fprintf('\n             - Matrix "input_dir" with at least three columns')
fprintf('\n               ordered as: inc, id_pers, id_firm.')
fprintf('\n')
fprintf('\nOUTPUTS:     KSS correction to AKM variance components.')
fprintf('\n')
fprintf('\nNOTES:       - Requires MATLAB R2011a or later due to changes')
fprintf('\n               in -ichol()- command and MATLAB R2016b or later')
fprintf('\n               for implicit expansion.')
fprintf('\n')
fprintf('\nCREDITS:     Please cite: Engbom, Niklas & Christian Moser &')
fprintf('\n             Jan Sauermann, 2021. "Firm Pay Dynamics," Working')
fprintf('\n             Paper.')
fprintf('\n')
fprintf('\nTIME STAMP:  October 28, 2021.')
fprintf('\n==============================================================')
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
clear v


fprintf('\n==============================================================')
fprintf('\n SET DIRECTORIES AND PARAMETERS')
fprintf('\n==============================================================')
fprintf('\n  --> Options for KSS correction (see description in leave_out_KSS.m or leave_out_KSS_Chris.m)\n')
leave_out_level = 'matches'; % 'obs' or 'matches'
type_algorithm = 'JLA'; % 'exact' or 'JLA;
simulations_JLA = 50; % default = 200
lincom_do = 0; % 0 or 1
Z_lincom = []; % matrix of observables to be used when projecting firm effects onto observables
labels_lincom = []; % vector of dimension rx1 that provides a label for each of the columns in Z_lincom

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'codes'); %this contains the main LeaveOut Routines.
path(path,'CMG'); % CMG package http://www.cs.cmu.edu/~jkoutis/cmg.html
[~,~] = evalc('installCMG(1)'); %installs CMG routine (silently) % CHRIS: replaced [result,output] with [~,~] on LHS so as to not store output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n  --> read parameters\n')
params_dir = [DIR_INPUT '/parameters_akm_kss.csv']; % Path to file containing AKM parameters.
file_params = fopen(params_dir);
data_params = fscanf(file_params, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [16 inf])';
fclose(file_params);
n_params = 1;
year_min = data_params(:, n_params); % Minimum year.
fprintf(['\nyear_min = ' num2str(year_min), '\n'])
n_params = n_params + 1;
year_max = data_params(:, n_params); % Maximum year.
fprintf(['\nyear_max = ' num2str(year_max), '\n'])
n_params = n_params + 1;
akm_age_poly_order = data_params(:, n_params); % Age polynomial order (0 = no age controls; 1 = restricted profile of age indicators; 2 / 3 / etc. = include 2nd order term / 2nd and 3rd order terms / etc.).
fprintf(['\nakm_age_poly_order = ' num2str(akm_age_poly_order), '\n'])
n_params = n_params + 1;
akm_age_flat_min = data_params(:, n_params); % Minimum age for which income-age profile is  restricted to be flat (only relevant if akm_age_poly_order == 1).
fprintf(['\nage_flat_min = ' num2str(akm_age_flat_min), '\n'])
n_params = n_params + 1;
akm_age_flat_max = data_params(:, n_params); % Maximum age for which income-age profile is restricted to be flat (only relevant if akm_age_poly_order == 1).
fprintf(['\nage_flat_max = ' num2str(akm_age_flat_max), '\n'])
n_params = n_params + 1;
akm_age_norm = data_params(:, n_params); % Normalization age.
fprintf(['\nage_norm = ' num2str(akm_age_norm), '\n'])
n_params = n_params + 1;
akm_age_min = data_params(:, n_params); % Minimum age.
fprintf(['\nage_min = ' num2str(akm_age_min), '\n'])
n_params = n_params + 1;
akm_age_max = data_params(:, n_params); % Maximum age.
fprintf(['\nage_max = ' num2str(akm_age_max), '\n'])
n_params = n_params + 1;
akm_year_dummies = data_params(:, n_params); % Year dummies (0 = no; 1 = yes).
fprintf(['\nakm_year_dummies = ' num2str(akm_year_dummies), '\n'])
n_params = n_params + 1;
akm_dem_inter = data_params(:, n_params); % Demographics (gender x education) interactions (0 = no; 1 = yes).
fprintf(['\ndem_inter = ' num2str(akm_dem_inter), '\n'])
n_params = n_params + 1;
akm_hours = data_params(:, n_params); % Hours controls (0 = no; 1 = yes).
fprintf(['\nakm_hours = ' num2str(akm_hours), '\n'])
n_params = n_params + 1;
akm_occ = data_params(:, n_params); % Occupation controls (0 = no; 1 = yes).
fprintf(['\nakm_occ = ' num2str(akm_occ), '\n'])
n_params = n_params + 1;
akm_tenure = data_params(:, n_params); % Tenure controls (0 = no; 1 = yes).
fprintf(['\nakm_tenure = ' num2str(akm_tenure), '\n'])
n_params = n_params + 1;
akm_model = data_params(:, n_params); % Estimate model with firm FEs (akm_model = 1) or model with firm-year FEs (akm_model = 2).
fprintf(['\nakm_model = ' num2str(akm_model), '\n'])
n_params = n_params + 1;
save_space = data_params(:, n_params); % Save space by first saving and later loading temporarily redundant variables (0 = do not save space; 1 = save space).
fprintf(['\nsave_space = ' num2str(save_space), '\n'])
n_params = n_params + 1;
parallel_max = data_params(:, n_params); % Maximum number of parallel workers.
fprintf(['\nparallel_max = ' num2str(parallel_max), '\n'])
if akm_age_poly_order == 1 && (akm_age_flat_min >= akm_age_flat_max || akm_age_flat_min >= akm_age_max || akm_age_flat_max <= akm_age_min)
    fprintf('\nUSER ERROR: Requested to estimate age effects (akm_age_poly_order = 1) but restricted-to-be-flat age region is too short or falls outside of selected age range.\n')
    exit
elseif akm_age_poly_order >= 2 && (akm_age_norm < akm_age_min || akm_age_norm > akm_age_max)
    fprintf('\nUSER ERROR: Requested to estimate higher-order age terms (akm_age_poly_order >= 2) but normalization age falls outside of selected age range.\n')
    exit
end
akm_model_f = (akm_model == 1);
akm_model_fy = (akm_model == 2);
clear params_dir file_params n_params akm_age_min akm_age_max data_params akm_model

fprintf('\n  --> set parallel environment\n')
delete(gcp('nocreate')); %clear parallel environment
c = parcluster('local');  %tell me # of available cores
ext = randi(999999, 1);
job_storage_loc = strcat(DIR_OUTPUT,'/ext_',num2str(ext));
eval(['mkdir ''', job_storage_loc, ''''])
c.JobStorageLocation = job_storage_loc; % recommended to avoid storing conflicting job information (see section "Running Multiple PCT Matlab Jobs" at https://rcc.uchicago.edu/docs/software/environments/matlab/)
nw = c.NumWorkers; %tell me # of available cores
if parallel_max == 0
    nw_max = round(3 + 34/(year_max - year_min + 1)); % OLD: ceil(120/(year_max - year_min + 1))
else
    nw_max = parallel_max;
end
nw = min(nw, nw_max); % restrict number of workers to be not too large -- inserted for SWE server, but maybe better to adjust/remove for other servers!?
parpool(c,nw,'IdleTimeout',Inf); % a total of nw cores will be assigned to MATLAB for parallelization
clear job_storage_loc ext c nw nw_max parallel_max

fprintf('\n  --> set read and write paths\n')
% namesrc = [DIR_INPUT, '/temp_AKM_KSS.csv']; %where original data is
filename=[DIR_OUTPUT '/results_AKM_KSS']; %output file name without extension (implicitly adding .csv)
input_dir = [DIR_INPUT '/tomatlab_akm_kss.csv']; % Path to file with Stata input to estimate AKM equation: income, person ID, firm ID, year, etc.
% output_dir = [DIR_OUTPUT '/tostata_akm_kss.txt']; % Path to file where MATLAB output is exported to.
done_dir_akm_kss = [DIR_OUTPUT '/done_akm_kss.txt']; % Path to file that confirms being done with execution execution of FUN_AKM_KSS.m.
success_dir_akm_kss = [DIR_OUTPUT '/success_akm_kss.txt']; % Path to file that confirms successful execution of FUN_AKM_KSS.m.
clear DIR_INPUT

fprintf('\n  --> start diary (log) file\n')
if akm_model_f
    FILE_LOG = [DIR_LOG, '/log_MATLAB_AKM_KSS_f_', num2str(year_min), '_', num2str(year_max), '.log'];
elseif akm_model_fy
    FILE_LOG = [DIR_LOG, '/log_MATLAB_AKM_KSS_fy_', num2str(year_min), '_', num2str(year_max), '.log'];
end
eval(['delete ''', FILE_LOG, ''''])
echo off all
eval(['diary ''', FILE_LOG, ''''])
fprintf('\n==============================================================')
fprintf('\n START LOG FILE FOR FUN_AKM_KSS.m')
fprintf('\n==============================================================')
clear DIR_LOG FILE_LOG


fprintf('\n==============================================================')
fprintf('\n LOAD DATA')
fprintf('\n==============================================================')
fprintf('\n  --> read input data\n')
input_file = fopen(input_dir, 'r');
input_format = '%f %f %f %u'; % always load at least 4 variables: inc, id_pers, id_firm, year
input_n = 4;
if akm_dem_inter
    input_format = [input_format, ' %u'];
    input_n = input_n + 1;
end
if akm_age_poly_order
    input_format = [input_format, ' %u'];
    input_n = input_n + 1;
end
if akm_hours
    input_format = [input_format, ' %u'];
    input_n = input_n + 1;
end
if akm_occ
    input_format = [input_format, ' %u'];
    input_n = input_n + 1;
end
if akm_tenure
    input_format = [input_format, ' %u'];
    input_n = input_n + 1;
end
fprintf('\n')
disp(['number of variables in input data = ' int2str(input_n) ' (format = ' input_format ')']);
input_data = fscanf(input_file, input_format, [input_n inf])';
fclose(input_file);
inc = input_data(:, 1); % Input data column 1: log income.
id_pers = input_data(:, 2); % Input data column 2: worker ID.
id_firm = input_data(:, 3); % Input data column 3: firm ID.
year = input_data(:, 4); % Input data column 4: year.
any_controls = (akm_year_dummies | akm_age_poly_order | akm_dem_inter | akm_hours | akm_occ | akm_tenure);
if any_controls
    col_counter = 4;
    if akm_dem_inter
        col_counter = col_counter + 1;
        dem = input_data(:, col_counter); % Additional input data: demographics (gender x education).
    end
    if akm_age_poly_order
        col_counter = col_counter + 1;
        age = input_data(:, col_counter); % Additional input data: age.
    end
    if akm_hours
        col_counter = col_counter + 1;
        hours = input_data(:, col_counter); % Additional input data: hours.
    end
    if akm_occ
        col_counter = col_counter + 1;
        occ = input_data(:, col_counter); % Additional input data: occupation codes.
    end
    if akm_tenure
        col_counter = col_counter + 1;
        tenure = input_data(:, col_counter); % Additional input data: tenure.
    end
end
clear input_dir input_file input_format input_n input_data col_counter
    

fprintf('\n==============================================================')
fprintf('\n ESTIMATE AKM EQUATION')
fprintf('\n==============================================================')
fprintf('\n  --> create unique indicators for workers, firms, years, etc.\n')
NT = length(inc); % Total number of worker-years.
[~, ~, id_pers_unique_index] = unique(id_pers);
results.N_workers_pop = max(id_pers_unique_index);
clear id_pers
% N = results.N_workers_pop; % Number of unique elements in id_pers = number of unique workers.
[~, ~, id_firm_unique_index] = unique(id_firm);
results.N_firms_pop = max(id_firm_unique_index);
id_firm_year = findgroups(id_firm, year);
[~, ~, id_firm_year_unique_index] = unique(id_firm_year);
results.N_firm_years_pop = max(id_firm_year_unique_index);
clear id_firm id_firm_year id_firm_year_unique_index
% J = max(id_firm_unique_index); % Number of unique elements in id_firm = number of unique firms.
if any_controls
    if akm_year_dummies
        [~, ~, year_unique_index] = unique(year);
        results.N_years = max(year_unique_index);
        T = results.N_years; % Number of unique elements in year = number of unique years.
    end
    if akm_age_poly_order == 1
        [~, ~, age_unique_index] = unique(age);
        results.N_ages_unrestricted = max(age_unique_index);
        clear age_unique_index
        age(age >= akm_age_flat_min & age <= akm_age_flat_max) = akm_age_flat_min; % Restrict income-age profile to be flat between ages akm_age_flat_min and akm_age_flat_max.
        clear akm_age_flat_min akm_age_flat_max
        [~, ~, age_unique_index] = unique(age); % Note: Command -unique()- returns unique values of age. Note: Column index (3rd assignment object) maps vector of unique entries (1st assignment object) into the original vector (argument of -unique()-).
        results.N_ages_restricted = max(age_unique_index);
        clear age
        N_Age = results.N_ages_restricted; % Number of unique elements in age = number of unique ages.
    end
    if akm_dem_inter
        [~, ~, dem_unique_index] = unique(dem);
        results.N_dems = max(dem_unique_index);
        clear dem
        N_Dem = results.N_dems; % Number of unique elements in dem = number of unique demographic (gender x education) groups.
    end
    if akm_hours
        [~, ~, hours_unique_index] = unique(hours);
        results.N_hours = max(hours_unique_index);
        clear hours
    %     N_H = results.N_hours; % Number of unique elements in hours = number of unique contractual hours levels.
    end
    if akm_occ
        [~, ~, occ_unique_index] = unique(occ);
        results.N_occs = max(occ_unique_index);
        clear occ
    %     N_O = results.N_occs; % Number of unique elements in occ = number of unique occupation codes.
    end
    if akm_tenure
        [~, ~, tenure_unique_index] = unique(tenure);
        results.N_tenure = max(tenure_unique_index);
        clear tenure
    %     N_Ten = results.N_tenure; % Number of unique elements in tenure = number of unique tenure levels.
    end
end

% fprintf('\n  --> create worker indicators\n')
% W = sparse(1:NT, id_pers_unique_index', 1); % Dimension: NT x N. Note: This is the only independent variable for which no category is dropped.
% W_len = size(W, 2);
% clear id_pers_unique_index
% 
% fprintf('\n  --> create firm indicators\n')
% F = sparse(1:NT, id_firm_unique_index', 1);
% F = F(:, 2:end); % Drop indicator for first firm, or else collinear with worker effects. Dimension: NT x (J - 1).
% F_len = size(F, 2);
% clear id_firm_unique_index

if any_controls
    fprintf('\n  --> create (demographics-specific) year indicators\n')
    if akm_year_dummies
        if akm_dem_inter
            Y = sparse(1:NT, T*(dem_unique_index' - 1) + year_unique_index', 1);
            Y(:, T*(dem_unique_index' - 1) + 1) = []; % Drop indicator for first year of every demographic group, or else collinear with worker effects. Dimension: NT x (T - 1)*N_Dem.
        else
            Y = sparse(1:NT, year_unique_index', 1);
            Y = Y(:, 2:end); % Drop indicator for first year, or else collinear with worker effects. Dimension: NT x (T - 1).
            % Y = Y(:, 2:end) - sparse(Y(:, end)*ones(1, size(Y, 2) - 1)); % Work in progress: use Deaton Method as an alternative to previous line: Subtract replication of last column from entire matrix, so year effects sum to zero.  Dimension: NT x (T - 1).
        end
        Y_len = size(Y, 2);
        clear year_unique_index T
    end

    fprintf('\n  --> create age indicators or higher-order (2nd and up) terms\n')
    if akm_age_poly_order == 1
        if akm_dem_inter
            A = sparse(1:NT, N_Age*(dem_unique_index' - 1) + age_unique_index', 1);
            A(:, N_Age*(dem_unique_index' - 1) + 1) = []; % Drop indicator for first age of every demographic group, or else collinear with worker effects. Dimension: NT x (N_Age - 1)*N_Dem.
        else
            A = sparse(1:NT, age_unique_index', 1);
            A = A(:, 2:end); % Drop indicator for first age, or else collinear with worker effects. Dimension: NT x (N_Age - 1).
        end
        clear age_unique_index N_Age dem_unique_index N_Dem
        A_len = size(A, 2);
    elseif akm_age_poly_order >= 2
        age = (age - akm_age_norm)/akm_age_norm;  % Rescale to avoid big numbers.
        if akm_dem_inter
            E = sparse(1:NT, dem_unique_index', 1);
            if v_old
                A = zeros(NT, N_Dem*(akm_age_poly_order - 1));            
                for n = 2:akm_age_poly_order
                    A(:, N_Dem*(n - 2) + 1:N_Dem*(n - 1)) = bsxfun(@times, E, age.^n);
                end
                A = sparse(A); % Dimension: NT x (N_Age - 1)*N_Dem.
            else
                A = repmat(E, 1, akm_age_poly_order - 1).*repmat(age, 1, N_Dem*(akm_age_poly_order - 1)).^kron((2:akm_age_poly_order), ones(1, N_Dem)); % Dimension: NT x (N_Age - 1)*N_Dem.
            end
            clear E age N_Dem dem_unique_index
        else
            if v_old
                A = zeros(NT, akm_age_poly_order - 1);
                for n = 2:akm_age_poly_order
                    A(:, n - 1) = age.^n;
                end
                A = sparse(A); % Dimension: NT x (N_Age - 1).
            else
                A = sparse(age.^(2:akm_age_poly_order)); % Dimension: NT x (N_Age - 1).
            end
            clear age
        end
        A_len = size(A, 2);
    end
    clear akm_age_norm v_old akm_dem_inter

    fprintf('\n  --> create hours indicators\n')
    if akm_hours
        H = sparse(1:NT, hours_unique_index', 1);
        H = H(:, 2:end); % Drop indicator for first hours level, or else collinear with worker fixed effects. Dimension: NT x (N_H - 1).
        H_len = size(H, 2);
        clear hours_unique_index
    end

    fprintf('\n  --> create occupation indicators\n')
    if akm_occ
        O = sparse(1:NT, occ_unique_index', 1);
        O = O(:, 2:end); % Drop indicator for first occupation code, or else collinear with worker fixed effects. Dimension: NT x (N_O - 1).
        O_len = size(O, 2);
        clear occ_unique_index
    end
    
    fprintf('\n  --> create tenure indicators\n')
    if akm_tenure
        Ten = sparse(1:NT, tenure_unique_index', 1);
        Ten = Ten(:, 2:end); % Drop indicator for first tenure level, or else collinear with worker fixed effects. Dimension: NT x (N_O - 1).
        Ten_len = size(Ten, 2);
        clear tenure_unique_index
    end
end

fprintf('\n  --> generate design matrix (X) and Gramian matrix (X_prime*X)\n')
% X = [W, F];
% X_len = W_len + F_len; % = N + (J - 1).
controls = [];
controls_len = 0; % = N + (J - 1).
if any_controls
    if akm_year_dummies
        controls = [controls, Y];
        controls_len = controls_len + Y_len;
        clear Y Y_len
    end
    if akm_age_poly_order >= 1
        controls = [controls, A];
        controls_len = controls_len + A_len;
        clear A A_len
    end
    if akm_hours
        controls = [controls, H];
        controls_len = controls_len + H_len;
        clear H H_len
    end
    if akm_occ
        controls = [controls, O];
        controls_len = controls_len + O_len;
        clear O O_len
    end
    if akm_tenure
        controls = [controls, Ten];
        controls_len = controls_len + Ten_len;
        clear Ten Ten_len
    end
    fprintf('\n')
    disp(['dimensions of the controls matrix (controls) = ' num2str(NT) 'x' num2str(controls_len)])
    clear akm_age_poly_order akm_hours akm_occ akm_tenure NT controls_len
end
clear any_controls

fprintf('\n  --> Run KSS code\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %CALL
                        
%Note: This is where you should specify where is the input .csv file
%      and how you would like to call the log-file created by the
%      leave out functions. 
%      
%      The user should also specify where to save and name:
%      1. Log File
%      2. Saved Results (will be in .csv)
        

%Make sure that the input .csv file is sorted by worker ID and year ID
%(xtset id_pers year in Stata). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data=importdata(namesrc);
% id_pers=data(:,1); % = worker ID
% id_firm=data(:,2); % = firm ID
% year=data(:,3); % = year ID
% y=data(:,4); % = log income
% clear data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %RUN LEAVE-OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Run
rng('default') % set seed to default (=0)
% [sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_COMPLETE(inc,id_pers,id_firm,leave_out_level,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE,type_of_algorithm,epsilon,filename);
% [sigma2_psi,sigma_psi_alpha,sigma2_alpha] = leave_out_KSS(inc,id_pers,id_firm,controls,leave_out_level,type_algorithm,simulations_JLA,lincom_do,Z_lincom,labels_lincom,filename);
% [~,~,~] = leave_out_KSS(inc,id_pers_unique_index,id_firm_unique_index,controls,leave_out_level,type_algorithm,simulations_JLA,lincom_do,Z_lincom,labels_lincom,filename);
[sigma2_psi,~,~,results] = leave_out_KSS_Chris(inc,id_pers_unique_index,id_firm_unique_index,controls,leave_out_level,type_algorithm,simulations_JLA,lincom_do,Z_lincom,labels_lincom,filename,akm_model_f,akm_model_fy,akm_year_dummies,year,results,save_space,DIR_OUTPUT);
clear year inc id_pers_unique_index id_firm_unique_index controls leave_out_level type_algorithm simulations_JLA lincom_do Z_lincom labels_lincom filename akm_year_dummies save_space

fprintf('\n  --> write summary of estimation results to output file in MATLAB format\n')
if akm_model_f
    save([DIR_OUTPUT, '/results_MATLAB_AKM_KSS_f_', num2str(year_min), '_', num2str(year_max), '.mat'], 'results')
elseif akm_model_fy
    save([DIR_OUTPUT, '/results_MATLAB_AKM_KSS_fy_', num2str(year_min), '_', num2str(year_max), '.mat'], 'results')
end
clear DIR_OUTPUT results year_min year_max akm_model_f akm_model_fy

fprintf('\n  --> confirm being done with execution of FUN_AKM_KSS.m\n')
done_file_akm_kss = fopen(done_dir_akm_kss, 'w');
out_header_format = '%10s';
out_header = {'done'};
out_data_format = '%10.0f';
clear done_dir_akm_kss
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(done_file_akm_kss, out_header_format, out_header{:});
fprintf(done_file_akm_kss, out_data_format, 1);
fclose(done_file_akm_kss);
clear done_file_akm_kss done_dir_akm_kss out_data_format out_header

fprintf('\n  --> confirm successful execution of FUN_AKM_KSS.m\n')
if exist('sigma2_psi','var')
    success_file_akm_kss = fopen(success_dir_akm_kss, 'w');
    out_header_format = '%10s';
    out_header = {'success'};
    out_data_format = '%10.0f';
    clear success_dir_akm_kss
    out_header_format = strcat(out_header_format, '\n');
    out_data_format = strcat(out_data_format, '\n');
    fprintf(success_file_akm_kss, out_header_format, out_header{:});
    fprintf(success_file_akm_kss, out_data_format, sigma2_psi);
    fclose(success_file_akm_kss);
    clear sigma2_psi success_file_akm_kss success_dir_akm_kss out_header out_header_format out_data_format
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

fprintf('\n  --> close parallel pool\n')
delete(gcp('nocreate')); % clear parallel environment

fprintf('\n  --> close diary (log) file\n')
fprintf('\n\n\n\n\n\n\n\n\n\n\n\n')
diary off

fprintf('\n  --> exit\n')
exit