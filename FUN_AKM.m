fprintf('\n==============================================================')
fprintf('\nDESCRIPTION: Estimates worker fixed effects, firm FEs, time')
fprintf('\n             FEs, returns to demographics, and additional')
fprintf('\n             controls based on Abowd, Kramarz, and Margolis')
fprintf('\n             (Econometrica, 1999) using the algorithm by Card,')
fprintf('\n             Heining, and Kline (QJE, 2013).')
fprintf('\n')
fprintf('\nINPUTS:      Two files in .csv format prepared in Stata:')
fprintf('\n             - A matrix "input_dir" with seven (7) columns or')
fprintf('\n               variables ordered as: income, person ID, firm ID,')
fprintf('\n               year, education, age, lagged firm ID.')
fprintf('\n')
fprintf('\nOUTPUTS:     File containing worker fixed effects, firm fixed')
fprintf('\n             effects, time fixed effects, and returns to')
fprintf('\n             demographics.')
fprintf('\n')
fprintf('\nNOTES:       - Requires MATLAB R2011a or later due to changes')
fprintf('\n               in -ichol()- command and MATLAB R2016b or later')
fprintf('\n               for implicit expansion.')
fprintf('\n')
fprintf('\nCREDITS:     Please cite the following paper:')
fprintf('\n             - Engbom, Niklas & Christian Moser & Jan')
fprintf('\n               Sauermann, 2021. "Firm Pay Dynamics," Working')
fprintf('\n               Paper.')
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
params_dir = [DIR_INPUT '/parameters_akm_kss.csv']; % Path to file containing AKM parameters.
file_params = fopen(params_dir);
data_params = fscanf(file_params, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', [17 inf])';
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
if akm_age_poly_order == 1 && (akm_age_flat_min >= akm_age_flat_max || akm_age_flat_min >= akm_age_max || akm_age_flat_max <= akm_age_min)
    fprintf('\nUSER ERROR: Requested to estimate age effects (akm_age_poly_order = 1) but restricted-to-be-flat age region is too short or falls outside of selected age range.\n')
    exit
elseif akm_age_poly_order >= 2 && (akm_age_norm < akm_age_min || akm_age_norm > akm_age_max)
    fprintf('\nUSER ERROR: Requested to estimate higher-order age terms (akm_age_poly_order >= 2) but normalization age falls outside of selected age range.\n')
    exit
end
akm_model_f = (akm_model == 1);
akm_model_fy = (akm_model == 2);
clear params_dir file_params akm_age_min akm_age_max data_params n_params akm_model

fprintf('\n  --> set read and write paths\n')
input_dir = [DIR_INPUT '/tomatlab_akm_kss.csv']; % Path to file with Stata input to estimate AKM equation: income, person ID, firm ID, year, etc.
output_dir = [DIR_OUTPUT '/tostata_akm.txt']; % Path to file where MATLAB output is exported to.
done_dir_akm = [DIR_OUTPUT '/done_akm.txt']; % Path to file that confirms being done with execution execution of FUN_AKM.m.
success_dir_akm = [DIR_OUTPUT '/success_akm.txt']; % Path to file that confirms successful execution of FUN_AKM.m.
clear DIR_INPUT

fprintf('\n  --> start diary (log) file\n')
if akm_model_f
    FILE_LOG = [DIR_LOG, '/log_MATLAB_AKM_f_', num2str(year_min), '_', num2str(year_max), '.log'];
    emp_str = 'firm';
elseif akm_model_fy
    FILE_LOG = [DIR_LOG, '/log_MATLAB_AKM_fy_', num2str(year_min), '_', num2str(year_max), '.log'];
    emp_str = 'firm-year';
end
eval(['delete ''', FILE_LOG, ''''])
echo off all
eval(['diary ''', FILE_LOG, ''''])
fprintf('\n==============================================================')
fprintf('\n START LOG FILE FOR FUN_AKM.m')
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
clear input_dir input_file input_format input_n input_data col_counter

fprintf('\n  --> summarize income distribution\n')
results = struct;
results.mean_y_pop = mean(inc);
results.var_y_pop = var(inc);
results.N_pop = length(inc);
disp(['Mean(ln(income)) = ' num2str(results.mean_y_pop)])
disp(['Var(ln(income)) = ' num2str(results.var_y_pop)])
disp(['Number of observations = ' num2str(results.N_pop)])


fprintf('\n==============================================================')
fprintf('\n ESTIMATE AKM EQUATION')
fprintf('\n==============================================================')
fprintf('\n  --> create unique indicators for workers, firms, years, etc.\n')
NT = length(inc); % Total number of worker-years.
id_pers_old = id_pers;
if save_space
    save([DIR_OUTPUT, '/id_pers_old.mat'], 'id_pers_old', '-v7.3')
    clear id_pers_old
end
[~, ~, id_pers_unique_index] = unique(id_pers);
results.N_workers_pop = max(id_pers_unique_index);
clear id_pers
% N = results.N_workers_pop; % Number of unique elements in id_pers = number of unique workers.
[~, ~, id_firm_unique_index] = unique(id_firm);
results.N_firms_pop = max(id_firm_unique_index);
if ~akm_model_f % if estimating model with firm FEs
    clear id_firm_unique_index
end
% J = results.N_firms_pop; % Number of unique elements in id_firm = number of unique firms.

id_firm_year = findgroups(id_firm, year);
[~, ~, id_firm_year_unique_index] = unique(id_firm_year);
results.N_firm_years_pop = max(id_firm_year_unique_index);
if ~akm_model_fy
    clear id_firm_year id_firm_year_unique_index
end
if ~akm_year_dummies
    clear id_firm_year
end
if akm_model_fy && akm_year_dummies
    clear id_firm_year_unique_index
end
% JT = results.N_firm_years_pop; % Number of unique elements in id_firm_year = number of unique firm-years.
clear id_firm
[~, ~, year_unique_index] = unique(year);
results.N_years = max(year_unique_index);
if (~akm_model_fy || ~akm_year_dummies) && save_space
    save([DIR_OUTPUT, '/year.mat'], 'year', '-v7.3')
    clear year
end
if akm_year_dummies
    T = results.N_years; % Number of unique elements in year = number of unique years.
else
    clear year_unique_index
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

fprintf('\n  --> create worker indicators\n')
W = sparse(1:NT, id_pers_unique_index', 1); % Dimension: NT x N. Note: This is the only independent variable for which no category is dropped.
W_len = size(W, 2);
clear id_pers_unique_index

fprintf(['\n  --> create ', emp_str, ' indicators\n'])
if akm_model_f % if estimating model with firm FEs
    F = sparse(1:NT, id_firm_unique_index', 1);
    F = F(:, 2:end); % Drop indicator for first firm, or else collinear with worker effects. Dimension: NT x (J - 1).
    clear id_firm_unique_index
elseif akm_model_fy && ~akm_year_dummies % if estimating model with non-normalized year firm-year FEs
    F = sparse(1:NT, id_firm_year_unique_index', 1);
    F = F(:, 2:end); % Drop indicator for first firm-year, or else collinear with worker effects. Dimension: NT x (JT - 1).
    clear id_firm_year_unique_index
elseif akm_model_fy && akm_year_dummies % if estimating model with normalized-to-be-mean-zero-each year firm-year FEs and year FEs
    loop_year_min = min(year);
    loop_year_max = max(year);
    for yy = loop_year_min:loop_year_max
        [~, ~, id_firm_year_unique_index_temp] = unique(id_firm_year.*(year == yy));
        F_temp = sparse(1:NT, id_firm_year_unique_index_temp', 1);
%         F_temp = F_temp(:, 2:end - 1) - F_temp(:, end); % Normalization that sums unweighted firm-year FEs to zero each year. Note: First column is redundant because it corresponds to id_firm_year_unique_index_temp for year ~= yy, and last column is the firm-year in a given year that is normalized.
        F_temp = F_temp(:, 2:end - 1) - sum(F_temp(:, 2:end - 1), 1).*(F_temp(:, end)/sum(F_temp(:, end), 1)); % Normalization that sums worker-weighted firm-year FEs to zero each year. Note: First column is redundant because it corresponds to id_firm_year_unique_index_temp for year ~= yy, and last column is the weighted firm-year in a given year that is normalized.
        if yy == loop_year_min
            F = F_temp;
        else
            F = cat(2, F, F_temp);
        end
    end
    clear id_firm_year_unique_index_temp F_temp id_firm_year loop_year_min loop_year_max yy
    if save_space
        save([DIR_OUTPUT, '/year.mat'], 'year', '-v7.3')
        clear year
    end
end
F_len = size(F, 2);

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
    clear age_unique_index N_Age N_Dem
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
        clear E age N_Dem
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
clear akm_age_norm v_old

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
    Ten = Ten(:, 2:end); % Drop indicator for first tenure level, or else collinear with worker fixed effects. Dimension: NT x (N_Ten - 1).
    Ten_len = size(Ten, 2);
    clear tenure_unique_index
end

fprintf('\n  --> clear demographic variable that is no longer needed\n')
if akm_dem_inter
    clear dem_unique_index
end

fprintf('\n  --> generate design matrix (X) and Gramian matrix (X_prime*X)\n')
X = [W, F];
X_len = W_len + F_len; % = N + (J - 1).
if save_space
    save([DIR_OUTPUT, '/W.mat'], 'W', '-v7.3')
    save([DIR_OUTPUT, '/F.mat'], 'F', '-v7.3')
    clear W F
end
if akm_year_dummies
    X = [X, Y];
    X_len = X_len + Y_len;
    if save_space
        save([DIR_OUTPUT, '/Y.mat'], 'Y', '-v7.3')
        clear Y
    end
end
if akm_age_poly_order >= 1
    X = [X, A];
    X_len = X_len + A_len;
    if save_space
        save([DIR_OUTPUT, '/A.mat'], 'A', '-v7.3')
        clear A
    end
end
if akm_hours
    X = [X, H];
    X_len = X_len + H_len;
    if save_space
        save([DIR_OUTPUT, '/H.mat'], 'H', '-v7.3')
        clear H
    end
end
if akm_occ
    X = [X, O];
    X_len = X_len + O_len;
    if save_space
        save([DIR_OUTPUT, '/O.mat'], 'O', '-v7.3')
        clear O
    end
end
if akm_tenure
    X = [X, Ten];
    X_len = X_len + Ten_len;
    if save_space
        save([DIR_OUTPUT, '/Ten.mat'], 'Ten', '-v7.3')
        clear Ten
    end
end
fprintf('\n')
disp(['dimensions of the design matrix (X) = ' num2str(NT) 'x' num2str(X_len)])
fprintf('\n')
disp(['dimensions of the Gramian matrix (X_prime*X) = ' num2str(X_len) 'x' num2str(X_len)])
fprintf('\n')
disp(['dimensions of the moment matrix (X_prime*y) = ' num2str(X_len) 'x' num2str(1)])
clear NT

fprintf('\n  --> perform incomplete Cholesky factorization of Gramian matrix (X_prime*X) as preconditioner for matrix inversion\n')
whos
X_prime_X = X'*X;
diagcomp_linear_n = 1;
diagcomp_linear_step = .1;
diagcomp_linear_list = diagcomp_linear_step.*(1:diagcomp_linear_n);
diagcomp_max = max(sum(abs(X_prime_X), 2)./diag(X_prime_X)) - 2; % value recommended by MATLAB to guarantee successful execution of -ichol()-, see https://www.mathworks.com/help/matlab/ref/ichol.html
diagcomp_candidate_n = 20;
diagcomp_candidate_base = 1.5;
diagcomp_overshoot = 5;
diagcomp_factor = 1./(diagcomp_candidate_base.^(diagcomp_candidate_n:-1:-diagcomp_overshoot));
diagcomp_candidates = diagcomp_max.*diagcomp_factor;
diagcomp_exponential_list = diagcomp_candidates(diagcomp_candidates > diagcomp_linear_step*diagcomp_linear_n);
diagcomp_list = [diagcomp_linear_list diagcomp_exponential_list];
for diagcomp = diagcomp_list
    try
        L = ichol(X_prime_X, struct('type', 'ict', 'droptol', 1e-2, 'diagcomp', diagcomp)); % Note: Command -ichol()- uses iterative approximation algorithm.
        fprintf('\n')
        disp(['NOTE: function -ichol()- with diagcomp = ' num2str(diagcomp) ' succeeded!'])
        break % exit for loop after successful evaluation of -ichol()-
    catch
        disp(['USER WARNING: Function -ichol()- with diagcomp = ' num2str(diagcomp) ' failed!'])
        if diagcomp == diagcomp_list(end)
            fprintf('\n')
            error(['USER ERROR: Function -ichol()- did not execute successfully for any value of diagcomp <= ' num2str(diagcomp)])
        end
    end
end
clear diagcomp diagcomp_candidate_base diagcomp_candidate_n diagcomp_candidates diagcomp_exponential_list diagcomp_factor diagcomp_linear_list diagcomp_linear_n diagcomp_linear_step diagcomp_list diagcomp_max diagcomp_overshoot

fprintf('\n  --> estimate coefficient vector\n')
b = pcg(X_prime_X, X'*inc, 1e-10, 2.5e3, L, L'); % Note: Command -pcg()- uses preconditioned conjugate gradients method to compute b = (X'*X)^-1*X'*inc. Dimension: N + (J - 1) + (T - 1)*N_Dem + ... if akm_dem_inter, or else N + (J - 1) + (T - 1) + ....
clear X_prime_X
fprintf('\n')
disp(['dimensions of the coefficient vector (b) = ' num2str(X_len) 'x' num2str(1)])
clear L X X_len

fprintf('\n  --> extract estimated paramters\n')
index_start = 1; % Beginning of index for worker fixed effects.
index_end = W_len; % End of index for worker fixed effects.
if save_space
    load([DIR_OUTPUT, '/W.mat'], 'W')
    eval(['delete ''', [DIR_OUTPUT, '/W.mat'], ''''])
end
pe = full(W*b(index_start:index_end)); % Worker fixed effects = N worker indicators * estimated worker returns vector. Dimension: NT x 1.
clear W W_len
index_start = index_end + 1; % Beginning of index for firm fixed effects.
index_end = index_end + F_len; % End of index for firm fixed effects.
if save_space
    load([DIR_OUTPUT, '/F.mat'], 'F')
    eval(['delete ''', [DIR_OUTPUT, '/F.mat'], ''''])
end
fe = full(F*b(index_start:index_end)); % Firm fixed effects (first one dropped) = J - 1 firm indicators (first one dropped) * estimated firm returns vector. Dimension: NT x 1.
clear F F_len
if akm_year_dummies
    index_start = index_end + 1; % Beginning of index for (demographics-specific) year fixed effects.
    index_end = index_end + Y_len; % End of index for (demographics-specific) year fixed effects.
    if save_space
        load([DIR_OUTPUT, '/Y.mat'], 'Y')
        eval(['delete ''', [DIR_OUTPUT, '/Y.mat'], ''''])
    end
    xb_y = full(Y*b(index_start:index_end)); % (Demographics-specific) Year fixed effects = T - 1 year indicators (first one dropped) * estimated year returns vector. Dimension: NT x N_Dem x 1.
    clear Y Y_len
end
if akm_age_poly_order >= 1
    index_start = index_end + 1; % Beginning of index for (demographics-specific) age fixed effects.
    index_end = index_end + A_len; % End of index for (demographics-specific) age fixed effects.
    if save_space
        load([DIR_OUTPUT, '/A.mat'], 'A')
        eval(['delete ''', [DIR_OUTPUT, '/A.mat'], ''''])
    end
    xb_a = full(A*b(index_start:index_end)); % Age fixed effects or higher-order age terms = N_Age - 1 age indicators or higher-order age terms * estimated age returns vector. Dimension: NT x 1.
    clear A A_len
end
if akm_hours
    index_start = index_end + 1; % Beginning of index for hours fixed effects.
    index_end = index_end + H_len; % End of index for hours fixed effects.
    if save_space
        load([DIR_OUTPUT, '/H.mat'], 'H')
        eval(['delete ''', [DIR_OUTPUT, '/H.mat'], ''''])
    end
    xb_h = full(H*b(index_start:index_end)); % Hours fixed effects.
    clear H H_len
end
if akm_occ
    index_start = index_end + 1; % Beginning of index for occupation fixed effects.
    index_end = index_end + O_len; % End of index for occupation fixed effects.
    if save_space
        load([DIR_OUTPUT, '/O.mat'], 'O')
        eval(['delete ''', [DIR_OUTPUT, '/O.mat'], ''''])
    end
    xb_o = full(O*b(index_start:index_end)); % Occupation fixed effects.
    clear O O_len
end
if akm_tenure
    index_start = index_end + 1; % Beginning of index for tenure fixed effects.
    index_end = index_end + Ten_len; % End of index for tenure fixed effects.
    if save_space
        load([DIR_OUTPUT, '/Ten.mat'], 'Ten')
        eval(['delete ''', [DIR_OUTPUT, '/Ten.mat'], ''''])
    end
    xb_t = full(Ten*b(index_start:index_end)); % Tenure fixed effects.
    clear Ten Ten_len
end
clear b index_start index_end


fprintf('\n==============================================================')
fprintf('\n ANALYZE ESTIMATION RESULTS')
fprintf('\n==============================================================')
fprintf(['\n  --> normalize worker, ', emp_str, ', and year fixed effect distributions such that worker fixed effects are mean zero\n'])
m_pe = mean(pe);
pe = pe - m_pe;
if akm_model_f || (akm_model_fy && ~akm_year_dummies)
    fe = fe + m_pe;
elseif akm_model_fy && akm_year_dummies
    xb_y = xb_y + m_pe;
end
clear m_pe

fprintf(['\n  --> means of worker and ', emp_str, ' fixed effects\n'])
m = mean([pe, fe]);
disp(m)
clear m

fprintf('\n  --> create AKM residual\n')
resid = inc - pe - fe;
if akm_year_dummies
    resid = resid - xb_y;
end
if akm_age_poly_order >= 1
    resid = resid - xb_a;
end
if akm_hours
    resid = resid - xb_h;
end
if akm_occ
    resid = resid - xb_o;
end
if akm_tenure
    resid = resid - xb_t;
end

fprintf('\n  --> full covariance matrix of estimated AKM wage components\n')
cov_mat = [inc, pe, fe];
clear inc
if akm_year_dummies
    cov_mat = [cov_mat, xb_y]; 
end
if akm_age_poly_order >= 1
    cov_mat = [cov_mat, xb_a];
end
if akm_hours
    cov_mat = [cov_mat, xb_h];
end
if akm_occ
    cov_mat = [cov_mat, xb_o];
end
if akm_tenure
    cov_mat = [cov_mat, xb_t];
end
cov_mat = [cov_mat, resid];
C = full(cov(cov_mat));
clear cov_mat resid
cov_str = [' income    person    ', emp_str, '  '];
if akm_year_dummies
    if akm_dem_inter
        cov_str = strcat(cov_str, '  dem-year');
    else
        cov_str = strcat(cov_str, '      year');
    end
end
if akm_age_poly_order >= 1
    if akm_dem_inter
        cov_str = strcat(cov_str, '  dem-age');
    else
        cov_str = strcat(cov_str, '      age');
    end
end
if akm_hours
    cov_str = strcat(cov_str, '   hours');
end
if akm_occ
    cov_str = strcat(cov_str, '     occup');
end
if akm_occ
    cov_str = strcat(cov_str, '     tenure');
end
cov_str = strcat(cov_str, '   residual');
fprintf('\n')
disp(cov_str)
disp(C)
var_y = C(1, 1);
var_y_share = 100*var_y/var_y;
var_pe = C(2, 2);
results.var_pe = var_pe;
var_pe_share = 100*var_pe/var_y;
var_fe = C(3, 3);
if akm_model_f
    results.var_fe = var_fe;
elseif akm_model_fy
    results.var_fye = var_fe;
end
var_fe_share = 100*var_fe/var_y;
if akm_model_f
    results.cov_pe_fe = C(2, 3);
elseif akm_model_fy
    results.cov_pe_fye = C(2, 3);
end
cov_sum_2 = var_y - var_pe - var_fe;
cov_index = 3;
if akm_year_dummies
    cov_index = cov_index + 1;
    var_year = C(cov_index, cov_index);
    results.var_year = var_year;
    var_year_share = 100*var_year/var_y;
    cov_sum_2 = cov_sum_2 - var_year;
end
if akm_age_poly_order >= 1
    cov_index = cov_index + 1;
    var_age = C(cov_index, cov_index);
    results.var_age = var_age;
    var_age_share = 100*var_age/var_y;
    cov_sum_2 = cov_sum_2 - var_age;
end
if akm_hours
    cov_index = cov_index + 1;
    var_hours = C(cov_index, cov_index);
    results.var_hours = var_hours;
    var_hours_share = 100*var_hours/var_y;
    cov_sum_2 = cov_sum_2 - var_hours;
end
if akm_occ
    cov_index = cov_index + 1;
    var_occ = C(cov_index, cov_index);
    results.var_occ = var_occ;
    var_occ_share = 100*var_occ/var_y;
    cov_sum_2 = cov_sum_2 - var_occ;
end
if akm_tenure
    cov_index = cov_index + 1;
    var_tenure = C(cov_index, cov_index);
    results.var_tenure = var_tenure;
    var_tenure_share = 100*var_tenure/var_y;
    cov_sum_2 = cov_sum_2 - var_tenure;
end
cov_index = cov_index + 1;
var_resid = C(cov_index, cov_index);
results.var_resid = var_resid;
var_resid_share = 100*var_resid/var_y;
cov_sum_2 = cov_sum_2 - var_resid;
results.cov_sum_2 = cov_sum_2;
cov_sum_2_share = 100*cov_sum_2/var_y;
clear cov_str C cov_index

fprintf('\n  --> display variance decomposition\n')
disp(['Var(income) = ' num2str(var_y) ' (' num2str(var_y_share) '%)'])
disp(['Var(worker FE) = ' num2str(var_pe) ' (' num2str(var_pe_share) '%)'])
disp(['Var(', emp_str, ' FE) = ' num2str(var_fe) ' (' num2str(var_fe_share) '%)'])
if akm_year_dummies
    if akm_dem_inter
        disp(['Var(dem-year FE) = ' num2str(var_year) ' (' num2str(var_year_share) '%)'])
    else
        disp(['Var(year FE) = ' num2str(var_year) ' (' num2str(var_year_share) '%)'])
    end
    clear var_year var_year_share
end
if akm_age_poly_order >= 1
    if akm_dem_inter
        disp(['Var(dem-age FE) = ' num2str(var_age) ' (' num2str(var_age_share) '%)'])
    else
        disp(['Var(age FE) = ' num2str(var_age) ' (' num2str(var_age_share) '%)'])
    end
    clear var_age var_age_share
end
if akm_hours
    disp(['Var(hours FE) = ' num2str(var_hours) ' (' num2str(var_hours_share) '%)'])
    clear var_hours var_hours_share
end
if akm_occ
    disp(['Var(occ FE) = ' num2str(var_occ) ' (' num2str(var_occ_share) '%)'])
    clear var_occ var_occ_share
end
if akm_tenure
    disp(['Var(tenure FE) = ' num2str(var_tenure) ' (' num2str(var_tenure_share) '%)'])
    clear var_tenure var_tenure_share
end
disp(['2*sum(Cov(.)) = ' num2str(cov_sum_2) ' (' num2str(cov_sum_2_share) '%)'])
disp(['Var(resid) = ' num2str(var_resid) ' (' num2str(var_resid_share) '%)'])
clear var_y var_y_share var_pe var_pe_share var_fe var_fe_share akm_dem_inter cov_sum_2 cov_sum_2_share var_resid var_resid_share

fprintf(['\n  --> correlation coefficient between AKM worker and ', emp_str, ' fixed effects\n']);
rho = corr(pe, fe);
if akm_model_f
    results.rho_pe_fe = rho;
elseif akm_model_fy
    results.rho_pe_fye = rho;
end
disp(['Corr(worker FE, ', emp_str, ' FE) = ' num2str(rho)])
fprintf('\n')
clear emp_str rho

fprintf('\n  --> write estimation results to output file in tab-delimited format\n')
output_file = fopen(output_dir, 'w');
out_header_format = '%10s\t %2s\t %6s\t %6s';
out_header = {'id_pers', 'y', 'pe', 'fe'};
out_data_format = '%10.0f\t %2.0f\t %6.4f\t %6.4f';
if save_space
    load([DIR_OUTPUT, '/id_pers_old.mat'], 'id_pers_old')
    eval(['delete ''', [DIR_OUTPUT, '/id_pers_old.mat'], ''''])
    load([DIR_OUTPUT, '/year.mat'], 'year')
    eval(['delete ''', [DIR_OUTPUT, '/year.mat'], ''''])
end
out_data = [id_pers_old'; year'; pe'; fe'];
clear save_space output_dir id_pers_old year pe
if akm_year_dummies
    out_header_format = strcat(out_header_format, '\t %6s');
    out_header{end + 1} = 'xb_y';
    out_data_format = strcat(out_data_format, '\t %6.4f');
    out_data = [out_data; xb_y'];
    clear xb_y
end
if akm_age_poly_order >= 1
    out_header_format = strcat(out_header_format, '\t %6s');
    out_header{end + 1} = 'xb_a';
    out_data_format = strcat(out_data_format, '\t %6.4f');
    out_data = [out_data; xb_a'];
    clear xb_a
end
if akm_hours
    out_header_format = strcat(out_header_format, '\t %6s');
    out_header{end + 1} = 'xb_h';
    out_data_format = strcat(out_data_format, '\t %6.4f');
    out_data = [out_data; xb_h'];
    clear xb_h
end
if akm_occ
    out_header_format = strcat(out_header_format, '\t %6s');
    out_header{end + 1} = 'xb_o';
    out_data_format = strcat(out_data_format, '\t %6.4f');
    out_data = [out_data; xb_o'];
    clear xb_o
end
if akm_tenure
    out_header_format = strcat(out_header_format, '\t %6s');
    out_header{end + 1} = 'xb_t';
    out_data_format = strcat(out_data_format, '\t %6.4f');
    out_data = [out_data; xb_t'];
    clear xb_t
end
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(output_file, out_header_format, out_header{:});
fprintf(output_file, out_data_format, out_data);
fclose(output_file);
clear out_header_format out_header out_data_format out_data output_file akm_age_poly_order akm_hours akm_occ akm_tenure akm_year_dummies

fprintf('\n  --> write summary of estimation results to output file in MATLAB format\n')
if akm_model_f
    save([DIR_OUTPUT, '/results_MATLAB_AKM_f_', num2str(year_min), '_', num2str(year_max), '.mat'], 'results')
elseif akm_model_fy
    save([DIR_OUTPUT, '/results_MATLAB_AKM_fy_', num2str(year_min), '_', num2str(year_max), '.mat'], 'results')
end
clear DIR_OUTPUT results year_min year_max akm_model_f akm_model_fy

fprintf('\n  --> confirm being done with execution of FUN_AKM.m\n')
done_file_akm = fopen(done_dir_akm, 'w');
out_header_format = '%10s';
out_header = {'done'};
out_data_format = '%10.0f';
clear done_dir_akm
out_header_format = strcat(out_header_format, '\n');
out_data_format = strcat(out_data_format, '\n');
fprintf(done_file_akm, out_header_format, out_header{:});
fprintf(done_file_akm, out_data_format, 1);
fclose(done_file_akm);
clear done_file_akm out_data_format out_header out_header_format

fprintf('\n  --> confirm successful execution of FUN_AKM.m\n')
if exist('fe','var')
    success_file_akm = fopen(success_dir_akm, 'w');
    out_header_format = '%10s';
    out_header = {'success'};
    out_data_format = '%10.0f';
    clear success_dir_akm
    out_header_format = strcat(out_header_format, '\n');
    out_data_format = strcat(out_data_format, '\n');
    fprintf(success_file_akm, out_header_format, out_header{:});
    fprintf(success_file_akm, out_data_format, fe);
    fclose(success_file_akm);
    clear success_file_akm fe out_data_format out_header out_header_format
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