function [sigma2_psi,sigma_psi_alpha,sigma2_alpha,results] = leave_out_KSS_Chris(y,id,firmid,controls,leave_out_level,type_algorithm,simulations_JLA,lincom_do,Z_lincom,labels_lincom,filename,akm_model_f,akm_model_fy,akm_year_dummies,year,results,save_space,DIR_OUTPUT)
    %% Author: Raffaele Saggio
    %Email: rsaggio@mail.ubc.ca
    
    %% EDITED BY CHRIS MOSER (October 12, 2021)

    %% DESCRIPTION
    %This function computes the bias-corrected components in a two-way model as
    %described in Kline, Saggio and Soelvsten (2020, ECTA -- KSS Henceforth).

    %This function can be applied to any two-way fixed effects model
    %(student-teachers, patient-doctors, etc). We use AKM jargon (workers, firms) 
    %when describing the code.

    %% WARNING!
    %This function only works properly if the user has properly sorted the data
    %by id-year (e.g. xtset id year in Stata).

    %% MANDATORY INPUTS
    %y: outcome. Dimensions: N* x 1; N*= # of person-year observations.
    %--
    %id: worker indicators. Dimensions: N* x 1
    %--
    %firmid: firm indicators. Dimensions: N* x 1

    %% NON-MANDATORY INPUTS    
    %leave_out_level: string variable that takes two values:

        %'obs': perform leave-out by leaving a person-year observation out (default)

        %'matches': perform leave-out by leaving an entire person-firm match out.

        %Default: `matches'.

    %controls:
        %Matrix of controls with dimensions: N* x P. This matrix of controls must
        %be appropriately defined by the user ex-ante. For instance, if the user 
        %wants to include time effects then the user should include in the matrix 
        %'controls' the set of dummy variables associated to a particular year
        %effect, making sure to avoid potential collinearity issues. 

        %These controls are going to be partialled out as follows. Let the
        %model be

        %y=D*alpha + F*psi + X*b + e

        %The code is going to estimate the above model, get an estimate of b
        %compute then ynew=y-Xb and run a variance decomposition on 

        %ynew=D*alpha + F*psi

    %type_algorithm: This takes two values: "exact" or "JLA".                    

        %   "exact": performs exact computation of (Bii,Pii). 

        %   "JLA": perform random projection methods to approximate (Bii,Pii) as 
        %   detailed in the Computational Appendix of KSS. 

        %In larger datasets, the user should always set type_algorithm='JLA'.

        %Default: JLA if # of observations in the original data is >10,000

    %simulations_JLA: a natural number.    

        %This governs the # of simulations in the JLA algorithm to approximate 
        %(Bii,Pii).

        %Default: 200.

    %lincom_do: binary. 

        %If =1, the code regresses the firm effects on the set of covariates
        %specied by Z_lincom and report the correct t-statistic
        %associated with this regression.

        %Default: 0

     %Z_lincom: matrix of regressors with dimension N* x r.

        %Matrix of observables to be used when projecting the firm effects into
        %observables.

     %labels_lincom: string

        %vector of dimension rx1 that provides a label for each of 
        %the columns in Z_lincom.

    %filename: string. 
        %Where saved results should be stored and named. Use name like 
        %"leave_out_results" and not "leave_out_results.csv" %Default is 'leave_out_estimates';

    %% OUTPUTS

    %sigma2_psi:      Variance of the firm effects.

    %sigma_psi_alpha: Covariance of the firm, person effects.

    %sigma2_alpha:    Variance of the person effects. 
    %                 When leaving a match-out, this parameter is only
    %                 identified among movers, i.e. individuals that moved b/w
    %                 firms across periods. The log file is going to report bounds
    %                 for this parameter + provide its estimate when leaving a
    %                 single observation out in the KSS procedure.
    %

    %The function  is also going to save on disk one .csv file and one .mat file. 
    %The .csv contains information on the leave-out connected set. 
    %First column reports the outcome variable, second and third columns the 
    %worker and the firm identifier. 
    %The fourth column reports the stastistical leverages. 
    %If the code is reporting a leave-out correction at the match-level, 
    %the .csv will be collapsed at the match level. 


    %%% MANUALLY-SET PARAMETERS BY CHRIS MOSER
    droptol = 1e-2;
    % diagcomp = 1e-1; % loop through values of this later on


    %% READ
    no_controls=0;
    no_algo=0;
    no_scale=0;
    no_labels=0;
    if nargin < 3
        error('More arguments needed');
    end

    if nargin == 3
        no_controls=1;
        controls=ones(size(y,1),1);
        leave_out_level='matches';
        no_algo=1;
        no_scale=1;
        lincom_do=0;
        no_labels=1;
        filename='leave_out_estimates';
    end

    if nargin == 4   
        leave_out_level='matches';
        no_algo=1;
        no_scale=1; 
        lincom_do=0;
        no_labels=1;
        filename='leave_out_estimates';
    end

    if nargin == 5
        no_algo=1;
        no_scale=1;   
        filename='leave_out_estimates';
        lincom_do=0;
        no_labels=1;
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
    end

    if nargin == 6    
        no_scale=1;   
        filename='leave_out_estimates'; 
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        lincom_do=0;
        no_labels=1;
    end

    if nargin == 7      
        filename='leave_out_estimates';
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        if size(simulations_JLA,2)==0
        no_scale=1;
        end
        lincom_do=0;
        no_labels=1;
    end

    if nargin == 8      
        filename='leave_out_estimates';
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        if size(simulations_JLA,2)==0
        no_scale=1;
        end
        if size(lincom_do,2)==0
        lincom_do=0;
        end
        if lincom_do==1 
        disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates')    
        lincom_do=0;
        end
        no_labels=1;
    end

    if nargin == 9      
        filename='leave_out_estimates';
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        if size(simulations_JLA,2)==0
        no_scale=1;
        end
        if size(lincom_do,2)==0
        lincom_do=0;
        end
        if lincom_do==1 && size(Z_lincom,2)==0 
        disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates,no lincom results will be shown')    
        lincom_do=0;
        end
        no_labels=1;

    end

    if nargin == 10      
        filename='leave_out_estimates';
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        if size(simulations_JLA,2)==0
        no_scale=1;
        end
        if size(lincom_do,2)==0
        lincom_do=0;
        end
        if lincom_do==1 && size(Z_lincom,2)==0 
        disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates,no lincom results will be shown')    
        lincom_do=0;
        end
        if size(labels_lincom,2)==0
        no_labels=1;
        end

    end

    if nargin == 11      
        if size(leave_out_level,2)==0
        leave_out_level='matches';
        end
        if size(type_algorithm,2)==0
        type_algorithm='JLA';
        end
        if size(simulations_JLA,2)==0
        no_scale=1;
        end
        if size(lincom_do,2)==0
        lincom_do=0;
        end
        if lincom_do==1 && size(Z_lincom,2)==0 
        disp('Warning: user wants to project the firm effects on some covariates Z but the user did not specify the set of covariates, no lincom results will be shown')    
        lincom_do=0;
        end
        if size(labels_lincom,2)==0
        no_labels=1;
        end

    end


    %Read empty controls
    if size(controls,2)==0
        no_controls=1;
        controls=ones(size(y,1),1);
    end

    %Read number of outputs
    if  nargout==1
        n_of_parameters=1;
    end

    if nargout==2
        n_of_parameters=2;
    end

    if nargout>=3
        n_of_parameters=3;
    end



    % automatically set number of simulations (Read # of FE)
    if no_scale == 1
        simulations_JLA=200; %default rule
    end
    clear no_scale
    
    % automatically set type of algorithm
    if no_algo == 1
        if size(y,1)>10000
            type_algorithm='JLA';
        end
        if size(y,1)<=10000
            type_algorithm='exact';
        end
    end
    clear no_algo

    if 1 == 1
    %Listing options
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp('Running KSS Correction with the following options')
    if strcmp(leave_out_level,'matches') 
    s='Leave Out Strategy: Leave match out';
    disp(s);
    end
    if strcmp(leave_out_level,'obs') 
    s='Leave Out Strategy: Leave person-year observation out';
    disp(s);
    end
    if strcmp(type_algorithm,'exact')
    s='Algorithm for Computation of Statistical Leverages: Exact';
    disp(s);
    end
    if strcmp(type_algorithm,'JLA')
    s=['Algorithm for Computation of Statistical Leverages: JLA with ' num2str(simulations_JLA) ' simulations.'];
    disp(s);
    end
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    end


    %% STEP 1: FIND CONNECTED SET
    %As first step in our analysis, we run estimation of a standard AKM model
    %on the largest connected set

    % CHRIS: create firm-year indicators, if requested
    [~,~,firmid_original]=unique(firmid);
    if akm_model_fy
        firmid = findgroups(firmid, year);
    end

    %Lagfirmid
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    clear gcs

    %Find connected set
    if lincom_do==1
    K=size(controls,2);
    controls=[controls Z_lincom];
    end
    clear Z_lincom
    
    results.mean_y_before_connected = mean(y);
    results.var_y_before_connected = var(y);
    results.N_before_connected = length(y);
    [~,~,id_unique]=unique(id);
    results.N_workers_before_connected = max(id_unique);
    clear id_unique
    [~,~,firmid_unique]=unique(firmid_original);
    results.N_firms_before_connected = max(firmid_unique);
    clear firmid_unique
    firmyearid = findgroups(firmid_original, year);
    [~,~,firmyearid_unique]=unique(firmyearid);
    results.N_firm_years_before_connected = max(firmyearid_unique);
    clear firmyearid firmyearid_unique
    
    [y,id,firmid,id_old,firmid_old,firmid_original,controls,year] = connected_set(y,id,firmid,lagfirmid,firmid_original,controls,year,akm_model_f);
    clear lagfirmid
    
    results.mean_y_connected = mean(y);
    results.var_y_connected = var(y);
    results.N_connected = length(y);
    [~,~,id_unique]=unique(id);
    results.N_workers_connected = max(id_unique);
    clear id_unique
    [~,~,firmid_unique]=unique(firmid_original);
    results.N_firms_connected = max(firmid_unique);
    clear firmid_unique
    firmyearid = findgroups(firmid_original, year);
    [~,~,firmyearid_unique]=unique(firmyearid);
    results.N_firm_years_connected = max(firmyearid_unique);
    clear firmyearid firmyearid_unique


    %% STEP 2: LEAVE ONE OUT CONNECTED SET
    %Here we compute the leave out connected set as defined in the Computational
    %Appendix of KSS. 
    %The input data is represented by the largest connected set. After applying
    %the function 'pruning_unbal_v3', the output data will be a connected set
    %such that the associated bipartite graph between workers and firms remains
    %connected after removing any particular worker from the graph.


    %Leave One Out Largest Connected Set
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp('SECTION 1')
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    [~,y,firmid,id,id_old,firmid_old,firmid_original,controls,~,year] = evalc('pruning_unbal_v3(y,firmid,id,id_old,firmid_old,firmid_original,controls,year)'); % Note: First output variable 'results' contains the text output from the output captured (i.e., not printed to display) by -evalc()-.
%     [y,firmid,id,id_old,firmid_old,firmid_original,controls,~,year] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,firmid_original,controls,year);
    
    results.mean_y_connected_leave_one_out = mean(y);
    results.var_y_connected_leave_one_out = var(y);
    results.N_connected_leave_one_out = length(y);
    [~,~,id_unique]=unique(id);
    results.N_workers_connected_leave_one_out = max(id_unique);
    clear id_unique
    [~,~,firmid_unique]=unique(firmid_original);
    results.N_firms_connected_leave_one_out = max(firmid_unique);
    clear firmid_unique
    firmyearid = findgroups(firmid_original, year);
    [~,~,firmyearid_unique]=unique(firmyearid);
    results.N_firm_years_connected_leave_one_out = max(firmyearid_unique);
    clear firmyearid firmyearid_unique

    %%%Drop stayers with a single person year observation
    T=accumarray(id,1);
    T=T(id);
    sel=T>1;
    clear T
    y=y(sel,:);
    firmid=firmid(sel,:);
    id=id(sel,:);
    id_old=id_old(sel,:);
    firmid_old=firmid_old(sel,:);
    firmid_original=firmid_original(sel,:);
    controls=controls(sel,:);
    year=year(sel);
    clear sel
    
    results.mean_y_connected_leave_one_out_singletons = mean(y);
    results.var_y_connected_leave_one_out_singletons = var(y);
    results.N_connected_leave_one_out_singletons = length(y);
    [~,~,id_unique]=unique(id);
    results.N_workers_connected_leave_one_out_singletons = max(id_unique);
    clear id_unique
    [~,~,firmid_unique]=unique(firmid_original);
    results.N_firms_connected_leave_one_out_singletons = max(firmid_unique);
    clear firmid_unique
    firmyearid = findgroups(firmid_original, year);
    [~,~,firmyearid_unique]=unique(firmyearid);
    results.N_firm_years_connected_leave_one_out_singletons = max(firmyearid_unique);
    clear firmyearid firmyearid_unique

    %Resetting ids one last time.
    [~,~,firmid]=unique(firmid);
    [~,~,id]=unique(id);
    
    %Important Auxiliaries
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    stayer=(firmid==lagfirmid);
    clear lagfirmid
    stayer(gcs==1)=1;
    clear gcs
    stayer=accumarray(id,stayer);
    T=accumarray(id,1);
    stayer=T==stayer;
%     T=T(id);
    clear T
    movers=stayer~=1;
    clear stayer
    movers=movers(id);
    [~,~,n]=unique(id(movers));
    clear movers
    Nmovers=max(n);
    clear n
    NT=size(y,1);
    D=sparse(1:NT,id',1);
    N=size(D,2);
    % F=sparse(1:NT,firmid',1); % OLD BY KSS
    % J=size(F,2); % OLD BY KSS
    J=max(firmid); % NEW BY CHRIS
%     var_den=var(y);
    clear id_mover
    %Summarize
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s)
    s='Info on the leave one out connected set:';
    disp(s);
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s)
    s=['mean wage: ' num2str(mean(y))];
    disp(s)
    s=['variance of wage: ' num2str(var(y))];
    disp(s)
    s=['# of Movers: ' num2str(Nmovers)];
    disp(s);
    if akm_model_f
        s=['# of Firms: ' num2str(max(firmid))];
    elseif akm_model_fy
        s=['# of Firm-Years: ' num2str(max(firmid))];
    end
    disp(s);
    s=['# of Person-Year Observations: ' num2str(size(y,1))];
    disp(s);
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s)
    clear Nmovers


    %% STEP 3: Residualizing and Collapsing
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp('SECTION 2')
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp(s);
    %If the model includes controls, we're going to estimate the coefficients
    %of those controls and partial them out as follows
    if lincom_do==1
       Z=controls(:,K+1:end);
       controls=controls(:,1:K); 
    end

    %Residualize
    no_controls_save = no_controls;
    if no_controls == 0
        if akm_model_f % if estimating model with firm FEs
            F=sparse(1:NT,firmid',1);
            F=F*[speye(J-1);sparse(-zeros(1,J-1))]; % N+JxN+J-1 restriction matrix --> Drop indicator for last firm, or else collinear with worker effects.
        elseif akm_model_fy && (akm_year_dummies == 0 || akm_year_dummies == 2) % if estimating model with non-normalized year firm-year FEs
            F = sparse(1:NT, firmid', 1);
            F = F(:, 2:end); % Drop indicator for first firm-year, or else collinear with worker effects.
        elseif akm_model_fy && akm_year_dummies == 1 % if estimating model with normalized-to-be-mean-zero-each year firm-year FEs and year FEs
            loop_year_min = min(year);
            loop_year_max = max(year);
            for yy = loop_year_min:loop_year_max
                [~, ~, firmid_temp] = unique(firmid.*(year == yy));
                F_temp = sparse(1:NT, firmid_temp', 1);
                F_temp = F_temp(:, 2:end - 1) - sum(F_temp(:, 2:end - 1), 1).*(F_temp(:, end)/sum(F_temp(:, end), 1)); % Normalization that sums worker-weighted firm-year FEs to zero each year. Note: First column is redundant because it corresponds to firmid_temp for year ~= yy, and last column is the weighted firm-year in a given year that is normalized.
                if yy == loop_year_min
                    F = F_temp;
                else
                    F = cat(2, F, F_temp);
                end
            end
            clear F_temp
        end
        X=[D,F,controls];
        clear D F controls
        xx=X'*X;
        xy=X'*y;
        diagcomp_linear_n = 1;
        diagcomp_linear_step = .1;
        diagcomp_linear_list = diagcomp_linear_step.*(1:diagcomp_linear_n);
        diagcomp_max = max(sum(abs(xx), 2)./diag(xx)) - 2; % value recommended by MATLAB to guarantee successful execution of -ichol()-, see https://www.mathworks.com/help/matlab/ref/ichol.html
        diagcomp_candidate_n = 20;
        diagcomp_candidate_base = 1.5;
        diagcomp_overshoot = 5;
        diagcomp_factor = 1./(diagcomp_candidate_base.^(diagcomp_candidate_n:-1:-diagcomp_overshoot));
        diagcomp_candidates = diagcomp_max.*diagcomp_factor;
        diagcomp_exponential_list = diagcomp_candidates(diagcomp_candidates > diagcomp_linear_step*diagcomp_linear_n);
        diagcomp_list = [diagcomp_linear_list diagcomp_exponential_list];
        Lchol = [];
        for diagcomp = diagcomp_list
            try
                Lchol=ichol(xx,struct('type','ict','droptol',droptol,'diagcomp',diagcomp));
                fprintf('\n')
                disp(['NOTE: function -ichol()- with diagcomp = ' num2str(diagcomp) ' succeeded!'])
                break % exit for loop after successful evaluation of -ichol()-
            catch
                disp(['USER WARNING: Function -ichol()- with diagcomp = ' num2str(diagcomp) ' failed!'])
                if diagcomp == diagcomp_list(end)
                    fprintf('\n')
                    disp(['USER WARNING: Function -ichol()- did not execute successfully for any value of diagcomp <= ' num2str(diagcomp)])
                end
            end
        end
        clear diagcomp_linear_n diagcomp_linear_step diagcomp_linear_list diagcomp_max diagcomp_candidate_n diagcomp_candidate_base diagcomp_overshoot diagcomp_factor diagcomp_candidates diagcomp_exponential_list diagcomp diagcomp_list
        if size(Lchol,1) > 0 % if one of the -ichol()- evaluations succeeded, then use preconditioner
            b=pcg(xx,xy,1e-10,2.5e3,Lchol,Lchol');
        else % else, run brute-force conjugate-gradient method
            b=pcg(xx,xy,1e-10,2.5e3);
        end
        clear xx xy Lchol
        y=y-X(:,N+J:end)*b(N+J:end); %variance decomposition will be based on this residualized outcome.
        clear X b
%         no_controls=1;
    end
    
    results.mean_y_residualized = mean(y);
    results.var_y_residualized = var(y);
    results.N_residualized = length(y);
    [~,~,id_unique]=unique(id);
    results.N_workers_residualized = max(id_unique);
    clear id_unique
    [~,~,firmid_unique]=unique(firmid_original);
    results.N_firms_residualized = max(firmid_unique);
    clear firmid_unique
    firmyearid = findgroups(firmid_original, year);
    [~,~,firmyearid_unique]=unique(firmyearid);
    results.N_firm_years_residualized = max(firmyearid_unique);
    clear firmyearid firmyearid_unique
    var_y=results.var_y_residualized; % ADDED BY CHRIS, 01/15/2021
    
    clear no_controls firmid_original year


    %Collapsing
    %If the user wants to run the KSS correction leaving a match out, we're
    %going to collapse the data down to match-means and weight this collapsed
    %data by the length of a given match. 

    peso=ones(size(y,1),1);
    if strcmp(leave_out_level,'matches') 
        y_py=y;
        if save_space
            save([DIR_OUTPUT, '/y_py.mat'], 'y_py', '-v7.3')
            clear y_py
        end
        [~,~,match_id] 		= unique([id firmid],'rows','stable');
        peso				= accumarray(match_id,1);
        id					= accumarray(match_id,id,[],@(x)mean(x));
        firmid			    = accumarray(match_id,firmid,[],@(x)mean(x));
        id_old              = accumarray(match_id,id_old,[],@(x)mean(x));
        if save_space
            save([DIR_OUTPUT, '/id_old.mat'], 'id_old', '-v7.3')
            clear id_old
        end
        firmid_old          = accumarray(match_id,firmid_old,[],@(x)mean(x));
        if save_space
            save([DIR_OUTPUT, '/firmid_old.mat'], 'firmid_old', '-v7.3')
            clear firmid_old
        end
        y					= accumarray(match_id,y,[],@(x)mean(x)); 
        if lincom_do==1
            Z               = accumarray(match_id,Z,[],@(x)mean(x));
        end
        clear match_id
    end


    %% STEP 4: COMPUTATION OF (Pii,Bii)
    %This is the computationally expensive part of the code where we compute
    %the terms (Pii,Bii) as defined in KSS.

    %Build Design
        NT=size(y,1);
        D=sparse(1:NT,id',1);
        if save_space
            save([DIR_OUTPUT, '/id.mat'], 'id', '-v7.3')
            clear id
        end
        
        % OPTION 1: original code by KSS
        F=sparse(1:NT,firmid',1);

        % OPTION 2: alternative code by Chris -- need to discuss with Raffa?!
    %     if akm_model_f % if estimating model with firm FEs
    %         F=sparse(1:NT,firmid',1);
    %         F=F*[speye(J-1);sparse(-zeros(1,J-1))]; % N+JxN+J-1 restriction matrix --> Drop indicator for last firm, or else collinear with worker effects.
    %     elseif akm_model_fy && (akm_year_dummies == 0 || akm_year_dummies == 2) % if estimating model with non-normalized year firm-year FEs
    %         F = sparse(1:NT, firmid', 1);
    %         F = F(:, 2:end); % Drop indicator for first firm-year, or else collinear with worker effects.
    %     elseif akm_model_fy && akm_year_dummies == 1 % if estimating model with normalized-to-be-mean-zero-each year firm-year FEs and year FEs
    %         loop_year_min = min(year);
    %         loop_year_max = max(year);
    %         for yy = loop_year_min:loop_year_max
    %             [~, ~, firmid_temp] = unique(firmid.*(year == yy));
    %             F_temp = sparse(1:NT, firmid_temp', 1);
    %             F_temp = F_temp(:, 2:end - 1) - sum(F_temp(:, 2:end - 1), 1).*(F_temp(:, end)/sum(F_temp(:, end), 1)); % Normalization that sums worker-weighted firm-year FEs to zero each year. Note: First column is redundant because it corresponds to firmid_temp for year ~= yy, and last column is the weighted firm-year in a given year that is normalized.
    %             if yy == loop_year_min
    %                 F = F_temp;
    %             else
    %                 F = cat(2, F, F_temp);
    %             end
    %         end
    %     end
    %     clear F_temp
    
    clear akm_year_dummies
        
        if save_space
            save([DIR_OUTPUT, '/firmid.mat'], 'firmid', '-v7.3')
            clear firmid
        end
        
        N=size(D,2);
        J=size(F,2);
        X=[D,-F]; %shaped in a pure Laplacian format (NOTE BY CHRIS: the pure Laplacian form does not impose any restrictions on coefficients, which aids computationally speed but might not be desirable from a model selection perspective).
        if save_space
            save([DIR_OUTPUT, '/D.mat'], 'D', '-v7.3')
            clear D
            save([DIR_OUTPUT, '/F.mat'], 'F', '-v7.3')
            clear F
        end
        
    %Weighting Matrices
        X_fe=[sparse(NT,N) X(:,N+1:end)];
        X_fe=repelem(X_fe,peso,1); %weight by lenght of the spell
        X_pe=[X(:,1:N) sparse(NT,J)];
        X_pe=repelem(X_pe,peso,1); %weight by lenght of the spell
        PESO_MAT=sparse(1:NT,(1:NT)',peso.^0.5,NT,NT);
        if save_space
            save([DIR_OUTPUT, '/peso.mat'], 'peso', '-v7.3')
            clear peso
        end
        y_untransformed=y;
        if save_space
            save([DIR_OUTPUT, '/y_untransformed.mat'], 'y_untransformed', '-v7.3')
            clear y_untransformed
        end
        X=PESO_MAT*X;% TO ACCOUNT FOR WEIGHTING (FGLS)
        y=PESO_MAT*y;%FGLS transformation
        if save_space
            save([DIR_OUTPUT, '/PESO_MAT.mat'], 'PESO_MAT', '-v7.3')
            clear PESO_MAT
            save([DIR_OUTPUT, '/y.mat'], 'y', '-v7.3')
            clear y
        end
        xx=X'*X;
        disp('Calculating the statistical leverages...')
        [~,Lchol] = evalc('cmg_sdd(xx)'); %preconditioner for Laplacian matrices. % Note: First output variable 'results' contains the text output from the output captured (i.e., not printed to display) by -evalc()-.
%         Lchol = cmg_sdd(xx); %preconditioner for Laplacian matrices.
        
    whos % ADDED BY CHRIS

    tic
    if n_of_parameters==1
        [Pii, Mii, correction_JLA, Bii_fe]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
    elseif n_of_parameters==2
        [Pii, Mii, correction_JLA, Bii_fe, Bii_cov]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
    elseif n_of_parameters==3
        [Pii, Mii, correction_JLA, Bii_fe, Bii_cov, Bii_pe]=leverages(X_fe,X_pe,X,xx,Lchol,type_algorithm,simulations_JLA);
    end
    clear X_fe X_pe X xx Lchol type_algorithm simulations_JLA

    disp('Done!')
    toc


    %% STEP 5: ESTIMATION OF VARIANCE COMPONENTS
    %We use the statistical leverages, Pii, and the Bii associated with a given
    %variance component to bias correct these quantities using the KSS approach
    
    if save_space
        load([DIR_OUTPUT, '/F.mat'], 'F')
        eval(['delete ''', [DIR_OUTPUT, '/F.mat'], ''''])
    end
    X_fe                = [sparse(NT,N) F];
    if save_space
        load([DIR_OUTPUT, '/peso.mat'], 'peso')
        eval(['delete ''', [DIR_OUTPUT, '/peso.mat'], ''''])
    end
    X_fe                = repelem(X_fe,peso,1); %weight by lenght of the spell
    if save_space
        load([DIR_OUTPUT, '/D.mat'], 'D')
        eval(['delete ''', [DIR_OUTPUT, '/D.mat'], ''''])
    end
    X_pe                = [D sparse(NT,J)];
    X_pe                = repelem(X_pe,peso,1); %weight by lenght of the spell
    S                   = speye(J-1);
    S                   = [S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
    X                   = [D,F*S]; %back to grounded Laplacian. 
    clear D F S J
    if save_space
        load([DIR_OUTPUT, '/PESO_MAT.mat'], 'PESO_MAT')
        eval(['delete ''', [DIR_OUTPUT, '/PESO_MAT.mat'], ''''])
    end
    X                   = PESO_MAT*X;% TO ACCOUNT FOR WEIGHTING (FGLS)
    clear PESO_MAT
    xx                  = X'*X;
    clear NT N

    diagcomp_linear_n = 1;
    diagcomp_linear_step = .1;
    diagcomp_linear_list = diagcomp_linear_step.*(1:diagcomp_linear_n);
    diagcomp_max = max(sum(abs(xx), 2)./diag(xx)) - 2; % value recommended by MATLAB to guarantee successful execution of -ichol()-, see https://www.mathworks.com/help/matlab/ref/ichol.html
    diagcomp_candidate_n = 20;
    diagcomp_candidate_base = 1.5;
    diagcomp_overshoot = 5;
    diagcomp_factor = 1./(diagcomp_candidate_base.^(diagcomp_candidate_n:-1:-diagcomp_overshoot));
    diagcomp_candidates = diagcomp_max.*diagcomp_factor;
    diagcomp_exponential_list = diagcomp_candidates(diagcomp_candidates > diagcomp_linear_step*diagcomp_linear_n);
    diagcomp_list = [diagcomp_linear_list diagcomp_exponential_list];
    for diagcomp = diagcomp_list
        try
            Lchol=ichol(xx,struct('type','ict','droptol',droptol,'diagcomp',diagcomp));
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
    clear diagcomp_linear_n diagcomp_linear_step diagcomp_linear_list diagcomp_max diagcomp_candidate_n diagcomp_candidate_base diagcomp_overshoot diagcomp_factor diagcomp_candidates diagcomp_exponential_list diagcomp diagcomp_list droptol
    
    if save_space
        load([DIR_OUTPUT, '/y.mat'], 'y')
        eval(['delete ''', [DIR_OUTPUT, '/y.mat'], ''''])
    end
    xy                  = X'*y;
    [b, ~]              = pcg(xx,xy,1e-10,2.5e3,Lchol,Lchol'); % second output variable that was purged: flag
    clear xx xy Lchol
    eta                 = y-X*b;
%     var_eta             = var(eta); % ADDED BY CHRIS, 01/15/2021
    eta_h				= eta./Mii; %Leave one out residual
    clear eta Mii
    sigma_i				= (y-mean(y)).*eta_h; %KSS estimate of individual variance.
    clear eta_h
    sigma_i				= sigma_i.*correction_JLA; %to adjust for non-linear bias induced by JLA.
    clear correction_JLA

    X_fe                = X_fe(:,1:end-1);
    X_pe                = X_pe(:,1:end-1);

    if strcmp(leave_out_level,'matches')
        if save_space
            load([DIR_OUTPUT, '/id.mat'], 'id')
            eval(['delete ''', [DIR_OUTPUT, '/id.mat'], ''''])
        end
        T               = accumarray(id,1);
        stayers         = T==1;
        clear T
        stayers         = stayers(id);
        if save_space
            load([DIR_OUTPUT, '/y_py.mat'], 'y_py')
            eval(['delete ''', [DIR_OUTPUT, '/y_py.mat'], ''''])
            load([DIR_OUTPUT, '/firmid.mat'], 'firmid')
            eval(['delete ''', [DIR_OUTPUT, '/firmid.mat'], ''''])
        end
        sigma_stayers   = sigma_for_stayers(y_py,id,firmid,peso,b);
        clear y_py id firmid
        sigma_i(stayers)= sigma_stayers(stayers);
    end
    clear leave_out_level

    if n_of_parameters==1
        [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
    end

    if n_of_parameters==2
        [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
        [sigma_alpha_psi_AKM, sigma_psi_alpha]  = kss_quadratic_form(sigma_i,X_fe,X_pe,b,Bii_cov);
    end

    if n_of_parameters==3
        [sigma_2_psi_AKM, sigma2_psi]           = kss_quadratic_form(sigma_i,X_fe,X_fe,b,Bii_fe);
        [sigma_alpha_psi_AKM, sigma_psi_alpha]  = kss_quadratic_form(sigma_i,X_fe,X_pe,b,Bii_cov);
        [sigma_2_alpha_AKM, sigma2_alpha]       = kss_quadratic_form(sigma_i,X_pe,X_pe,b,Bii_pe);
        explained_var_AKM                       = (sigma_2_psi_AKM+2*sigma_alpha_psi_AKM+sigma_2_alpha_AKM)/var_y; % ADDED BY CHRIS, 01/15/2021
        explained_var_leave_out                 = (sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var_y; % ADDED BY CHRIS, 01/15/2021
    end
    
    clear X_pe b
    
    % AKM plug-in results:
    if akm_model_f
        results.var_fe = sigma_2_psi_AKM;
        results.cov_pe_fe = sigma_alpha_psi_AKM;
    elseif akm_model_fy
        results.var_fye = sigma_2_psi_AKM;
        results.cov_pe_fye = sigma_alpha_psi_AKM;
    end
    results.var_pe = sigma_2_alpha_AKM;
    if n_of_parameters==3
        rho = sigma_alpha_psi_AKM/(sqrt(sigma_2_psi_AKM)*sqrt(sigma_2_alpha_AKM));
        if akm_model_f
            results.rho_pe_fe = rho;
        elseif akm_model_fy
            results.rho_pe_fye = rho;
        end
        R2 = explained_var_AKM;
        results.R2 = R2;
    end
    
    % KSS bias-corrected results:
    if akm_model_f
        results.var_fe_KSS = sigma2_psi;
        results.cov_pe_fe_KSS = sigma_psi_alpha;
    elseif akm_model_fy
        results.var_fye_KSS = sigma2_psi;
        results.cov_pe_fye_KSS = sigma_psi_alpha;
    end
    results.var_pe_KSS = sigma2_alpha;
    if n_of_parameters==3
        rho_KSS = sigma_psi_alpha/(sqrt(sigma2_psi)*sqrt(sigma2_alpha));
        if akm_model_f
            results.rho_pe_fe_KSS = rho_KSS;
        elseif akm_model_fy
            results.rho_pe_fye_KSS = rho_KSS;
        end
        R2_KSS = explained_var_leave_out;
        results.R2_KSS = R2_KSS;
    end
    

    %% STEP 6: PRINTING RESULTS
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp(s);
    disp('SECTION 3')
    s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
    disp(s);
    disp(s);
    fprintf('\n\n')
    if no_controls_save==0 % ADDED BY CHRIS, 01/15/2021
        s=['Variance of Residualized Log Income (cond. on controls): ' num2str(var_y)]; % ADDED BY CHRIS, 01/15/2021
        disp(s); % ADDED BY CHRIS, 01/15/2021
    else % ADDED BY CHRIS, 01/15/2021
        s=['Variance of Log Income (no controls): ' num2str(var_y)]; % ADDED BY CHRIS, 01/15/2021
        disp(s); % ADDED BY CHRIS, 01/15/2021
    end % ADDED BY CHRIS, 01/15/2021
    clear no_controls_save var_y
%     s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
%     disp(s);
%     disp(s);
    fprintf('\n\n')
    s='BIASED PLUG-IN ESTIMATES';
    disp(s)
    if akm_model_f
        s=['Variance of Firm Effects: ' num2str(sigma_2_psi_AKM)];
    elseif akm_model_fy
        s=['Variance of Firm-Year Effects: ' num2str(sigma_2_psi_AKM)];
    end
    disp(s)
    if n_of_parameters>=2
        if akm_model_f
            s=['Covariance of Firm, Person Effects: ' num2str(sigma_alpha_psi_AKM)];
        elseif akm_model_fy
            s=['Covariance of Firm-Year, Person Effects: ' num2str(sigma_alpha_psi_AKM)];
        end
        disp(s);
    end
    if n_of_parameters==3
        s=['Variance of Person Effects: ' num2str(sigma_2_alpha_AKM)];
        disp(s);
        if akm_model_f
            s=['Correlation of Firm, Person Effects: ', num2str(rho)];
        elseif akm_model_fy
            s=['Correlation of Firm-Year, Person Effects: ', num2str(rho)];
        end
        disp(s);
    %     s=['Explained Variance Share (R^2) - Plugin: ' num2str((sigma_2_psi_AKM+2*sigma_alpha_psi_AKM+sigma_2_alpha_AKM)/var_den)];
        s=['Explained Variance Share (R^2) - Plugin: ', num2str(R2)]; % ADDED BY CHRIS, 01/15/2021
        disp(s);
        disp(' ')
    end
%     s='-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*';
%     disp(s);
%     disp(s);
    fprintf('\n\n')
    s='BIAS-CORRECTED KSS ESTIMATES';
    disp(s)
    if akm_model_f
        s=['Variance of Firm Effects: ' num2str(sigma2_psi)];
    elseif akm_model_fy
        s=['Variance of Firm-Year Effects: ' num2str(sigma2_psi)];
    end
    disp(s);
    if n_of_parameters>=2
        if akm_model_f
            s=['Covariance of Firm, Person Effects: ' num2str(sigma_psi_alpha)];
        elseif akm_model_fy
            s=['Covariance of Firm-Year, Person Effects: ' num2str(sigma_psi_alpha)];
        end
        disp(s);
    end
    if n_of_parameters==3
        s=['Variance of Person Effects: ' num2str(sigma2_alpha)];
        disp(s);
        if akm_model_f
            s=['Correlation of Firm, Person Effects: ' num2str(rho_KSS)];
        elseif akm_model_fy
            s=['Correlation of Firm-Year, Person Effects: ' num2str(rho_KSS)];
        end
        disp(s);
    %     s=['Explained Variance Share (R^2) - Leave-Out: ' num2str((sigma2_psi+2*sigma_psi_alpha+sigma2_alpha)/var_den)];
        s=['Explained Variance Share (R^2) - Leave-Out: ' num2str(R2_KSS)]; % ADDED BY CHRIS, 01/15/2021
        disp(s);
        fprintf('\n\n')
    end
    clear n_of_parameters


    %% STEP 7: SAVING OUTPUT
    %Export csv with data from leave out connected set
    if save_space
        load([DIR_OUTPUT, '/y_untransformed.mat'], 'y_untransformed')
        eval(['delete ''', [DIR_OUTPUT, '/y_untransformed.mat'], ''''])
        load([DIR_OUTPUT, '/id_old.mat'], 'id_old')
        eval(['delete ''', [DIR_OUTPUT, '/id_old.mat'], ''''])
        load([DIR_OUTPUT, '/firmid_old.mat'], 'firmid_old')
        eval(['delete ''', [DIR_OUTPUT, '/firmid_old.mat'], ''''])
    end
    out=[y_untransformed,id_old,firmid_old,full(Pii)];
    s=[filename '.csv'];
    dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
    clear y_untransformed id_old Pii filename DIR_OUTPUT


    %% STEP 8: LINCOM
    if lincom_do == 1 
        Z                   = repelem(Z,peso,1);
        Transform           = X_fe; %regress the firm effects (+ regression is person-year weighted)
        if akm_model_f
            disp('Regressing the firm effects on observables...')
        elseif akm_model_fy
            disp('Regressing the firm-year effects on observables...')
        end
        if no_labels == 0
            [~]                 = lincom_KSS(y,X,Z,Transform,sigma_i,labels_lincom);
        elseif no_labels == 1
            [~]                 = lincom_KSS(y,X,Z,Transform,sigma_i);
        end
    end
    clear lincom_do akm_model_f akm_model_fy peso X_fe y sigma_i no_labels labels_lincom

end

