% close all; clear all; clc
dbstop if error
try_path_1 = 'P:/2019/107'; % IFAU/Uppsala server input/output directory
try_path_2 = '/Users/cm3594/Dropbox (CBS)/AKM Conference/5_results/temp'; % local (iMac Pro at the office)
try_path_3 = '/Users/economoser/Dropbox (CBS)/AKM Conference/5_results/temp'; % local (iMac Pro at home)
try_path_4 = '/Users/niklasengbom/Dropbox/AKM'; % local (iMac Pro at home)
if exist(try_path_1, 'dir') % server
    directories.code = [ try_path_1 '/1 Code/' ] ;
    directories.data = [ try_path_1 '/2 Data/' ] ;
    directories.data2 = [ try_path_1 '/4 Export/out/' ] ;
    directories.data3 = [ try_path_1 '/4 Export/csv/' ] ;
    directories.tables = [ try_path_1  '/4 Export/tex/' ] ;
    addpath(try_path_1)
    addpath(directories.code)
    addpath(directories.data)
    addpath(directories.tables)
elseif exist(try_path_2, 'dir') % local (Columbia work iMac Pro)
    directories.tables = try_path_2 ;
elseif exist(try_path_3, 'dir') % local (Columbia work MacBook)
    directories.tables = try_path_3 ;
elseif exist(try_path_4, 'dir') % Nik personal
    directories.data1 = [ try_path_4 '/4_logs/2021.11.05/'] ;
    directories.data2 = [ try_path_4 '/5_results/10_IFAU_20211209/out/'] ;
    directories.data3 = [ try_path_4 '/5_results/10_IFAU_20211209/csv/'] ;
    directories.tables = [ try_path_4 '/5_results/export/tex/'] ;
else
    error('\nUSER ERROR: Directory not found.\n')
end
clear try_path_1 try_path_2 try_path_3 try_path_4




%% TABLE 1: SUMMARY STATS
% read data from stata 
data = csv2mat_numeric([directories.data2 'Table1_workers.out']) ;
% write table
fid = fopen([directories.tables 'Table1.tex'],'w');
fprintf(fid,'\\begin{tabular}{l c c c}\n');
fprintf(fid,'\\hline \\hline \n');
fprintf(fid,' && Mean & St.d. \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\\\[-.15in] \n');
fprintf(fid,'\\multicolumn{4}{l}{\\textit{Panel A. Worker-level variables}} \\\\ \n');
clear output
output.age = 'Worker age (years)' ;
output.college = 'Share with college degree' ;
output.female = 'Share female' ;
output.earnings = 'Monthly earnings (log SEK)' ;
for moment = fieldnames(data)'
    if ~contains(moment{1},'_sd')
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f \\\\ \n',output.(moment{1}), ...
                                                 data.(moment{1}), ...
                                                 data.([moment{1} '_sd'])) ;
    end
end
data = csv2mat_numeric([directories.data2 'Table1_firms.out']) ;
fprintf(fid,'\\\\[.0in] \n');
fprintf(fid,'\\multicolumn{4}{l}{\\textit{Panel B. Firm-level variables}} \\\\ \n');
clear output
output.assets = 'Log capital' ;
output.f_size = 'Log number of workers' ;
output.va_pw = 'Labor value added per worker' ;
% output.assets_pw = 'Capital per worker, $\log ( K/N )$' ;
for moment = fieldnames(output)'
    if ~contains(moment{1},'_sd')
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f \\\\ \n',output.(moment{1}), ...
                                                 data.(moment{1}), ...
                                                 data.([moment{1} '_sd'])) ;
    end
end
data = csv2mat_numeric([directories.data2 'Table1_observations.out']) ;
fprintf(fid,'\\\\[.0in] \n');
fprintf(fid,'\\multicolumn{4}{l}{\\textit{Panel C. Observations}} \\\\ \n');
clear output
output.workeryears = 'Number of worker-years' ;
output.workers = 'Number of unique workers' ;
output.firmyears = 'Number of firm-years' ;
output.firms = 'Number of unique firms' ;
for moment = fieldnames(output)'
        fprintf(fid,'\\hspace{.05in} %s && %s \\\\ \n',output.(moment{1}), ...
                                         comma_separator(data.(moment{1}))) ;
end
fprintf(fid,'\\\\[-.15in] \n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);















%% TABLE 2: AUTOCORRELATION OF FIRM MEAN PAY
% read data from stata 
data = csv2mat_numeric([directories.data2 'Table4.out']) ;
% write table
for var = { 'earnings_demeaned' }
    for spec = { '' , '_b' }
        fid = fopen([directories.tables 'Table2' var{1} spec{1} '.tex'],'w');
        fprintf(fid,'\\setlength{\\tabcolsep}{3.5pt} \n') ;
        fprintf(fid,'\\begin{tabular}{l ccccc ccccc ccccc ccccc ccccc ccccc cc}\n');
        fprintf(fid,'\\hline \\hline \n');
        fprintf(fid,'\\\\[-.1in] \n');
        for a=1985:2016
            fprintf(fid,' & %4.0f',a);
        end
        fprintf(fid,'\\\\ \n');
        fprintf(fid,'\\hline \n');
        fprintf(fid,'\\\\[-.08in] \n');        
        for version  = { '' , '_w' }
            if strcmp(version{1},'')
                fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel A. Unweighted}} \\\\ \n');
            else
                fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel B. Weighted}} \\\\ \n');
            end
            for a = 1985:2016
                fprintf(fid,'%4.0f',a) ;
                for f = 1985:2016
                    if f<a
                        fprintf(fid,' & ') ;
                    else
                        fprintf(fid,' & %4.3f',data.([ var{1} '_cor' spec{1} version{1} num2str(a)])(f-a+1)) ;
                    end
                end
                fprintf(fid,'\\\\ \n');
            end
            fprintf(fid,'\\\\[-.08in] \n');
        end
        fprintf(fid,'\\hline\n');
        fprintf(fid,'\\end{tabular}');
        fclose(fid);
        
    end
end



















%% TABLE 2: AKM RESULTS
% read data from stata 
clear akm
for s = { 'f' , 'fy' , 'KSS_f' 'KSS_fy' }
    load([directories.data 'results_MATLAB_AKM_' s{1} '_1985_2016.mat']) ;
    akm.(s{1}) = results ;
end
% need to add these fields manually to be able to loop over them
akm.KSS_fy.var_fe = akm.KSS_fy.var_fye ;
akm.KSS_fy.var_fe_KSS = akm.KSS_fy.var_fye_KSS ;
akm.KSS_fy.cov_pe_fe = akm.KSS_fy.cov_pe_fye ;
akm.KSS_fy.cov_pe_fe_KSS = akm.KSS_fy.cov_pe_fye_KSS ;
if ~isfield(akm.KSS_fy,'rho_pe_fe')
    akm.KSS_fy.rho_pe_fe = akm.KSS_fy.rho_pe_fye ;
    akm.KSS_fy.rho_pe_fe_KSS = akm.KSS_fy.rho_pe_fye_KSS ;
end
    
% write table
fid = fopen([directories.tables 'Table3.tex'],'w');
fprintf(fid,'\\begin{tabular}{l c cc c cc}\n');
fprintf(fid,'\\hline \\hline \n');
fprintf(fid,'\\\\[-.15in] \n');
fprintf(fid,'&& \\multicolumn{2}{c}{Plug-in estimates} && \\multicolumn{2}{c}{KSS bias-correction} \\\\ \n');
fprintf(fid,'\\cline{3-4} \\cline{6-7} \\\\[-.15in] \n');
fprintf(fid,'&& Firm & Firm-year && Firm & Firm-year \\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\\\[-.15in] \n');
clear output
output.var_y_residualized = 'Var($y_{ijt}$)' ;
output.var_pe = 'Var($\widehat{\alpha}_{i}$)    ' ;
output.var_fe = 'Var($\widehat{\psi}_{jt}$)' ;
output.cov_pe_fe = '2$*$Cov($\widehat{\alpha}_{i},\widehat{\psi}_{jt}$)' ;
fprintf(fid,'\\multicolumn{7}{l}{\\textit{Panel A. AKM variance decomposition}} \\\\ \n');
for moment = fieldnames(output)'
    if strcmp(moment{1},'var_y_residualized')
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f && %4.3f & %4.3f \\\\[.02in] \n',output.(moment{1}), ...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1}),...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1})) ;
        fprintf(fid,'\\hspace{.05in} \\textit{(\\%% of total)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} \\\\[.02in] \n', ...
                     100, ...
                     100, ...
                     100, ...
                     100) ;
    else
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f && %4.3f & %4.3f \\\\[.02in] \n',output.(moment{1}), ...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1}),...
                     akm.KSS_f.([moment{1} '_KSS']),akm.KSS_fy.([moment{1} '_KSS'])) ;
        fprintf(fid,'\\hspace{.05in} \\textit{(\\%% of total)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} \\\\[.02in] \n', ...
                     100*akm.KSS_f.(moment{1})./akm.KSS_f.var_y_residualized, ...
                     100*akm.KSS_fy.(moment{1})./akm.KSS_f.var_y_residualized,...
                     100*akm.KSS_f.([moment{1} '_KSS'])./akm.KSS_f.var_y_residualized, ...
                     100*akm.KSS_fy.([moment{1} '_KSS'])./akm.KSS_f.var_y_residualized) ;
    end
end
fprintf(fid,'\\\\[.0in] \n');
clear output
fprintf(fid,'\\multicolumn{7}{l}{\\textit{Panel B. Descriptive statistics}} \\\\ \n');
output.rho_pe_fe = 'Corr($\widehat{\alpha}_{i},\widehat{\psi}_{jt}$)' ;
for moment = fieldnames(output)'
    if strcmp(moment{1},'var_y_residualized')
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f && %4.3f & %4.3f \\\\[.02in] \n',output.(moment{1}), ...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1}),...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1})) ;
    else
        fprintf(fid,'\\hspace{.05in} %s && %4.3f & %4.3f && %4.3f & %4.3f \\\\[.02in] \n',output.(moment{1}), ...
                     akm.KSS_f.(moment{1}),akm.KSS_fy.(moment{1}),...
                     akm.KSS_f.([moment{1} '_KSS']),akm.KSS_fy.([moment{1} '_KSS'])) ;
    end
end
clear output
% output.N_before_connected = 'Worker-years' ;
% output.N_firm_years_before_connected = 'Firm-years' ;
% output.N_workers_before_connected = 'Unique workers' ;
% output.N_firms_before_connected  = 'Unique firms' ;
output.N = 'Worker-years' ;
output.N_firm_years = 'Firm-years' ;
output.N_workers = 'Unique workers' ;
output.N_firms = 'Unique firms' ;
div.N = 1.e00 ;
div.N_firm_years = 1.e00 ;
div.N_workers = 1.e00 ;
div.N_firms = 1.e00 ;
version = 'pop' ;
version_kss= 'connected_leave_one_out_singletons' ;
for moment = fieldnames(output)'
        str = 0.000 ;
%         share_f = 100*akm.KSS_f.([ moment{1} '_' version]) ./ akm.KSS_f.([ moment{1} '_' version_base]) ;
%         share_fy = 100*akm.KSS_fy.([ moment{1} '_' version]) ./ akm.KSS_fy.([ moment{1} '_' version_base]) ;
        fprintf(fid,'\\hspace{.05in} %s && %s & %s && %s & %s \\\\ \n',output.(moment{1}), ...
                                         comma_separator(round(akm.f.([ moment{1} '_' version])./div.(moment{1}))) , ...
                                         comma_separator(round(akm.fy.([ moment{1} '_' version])./div.(moment{1}))) , ...
                                         comma_separator(round(akm.KSS_f.([ moment{1} '_' version_kss])./div.(moment{1}))) , ...
                                         comma_separator(round(akm.KSS_fy.([ moment{1} '_' version_kss])./div.(moment{1}))) );
        fprintf(fid,'\\hspace{.05in} \\textit{(\\%% of population)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} && \\textit{(%4.1f\\%%)} & \\textit{(%4.1f\\%%)} \\\\[.02in] \n', ...
                                         str, str , str , str );
                                     
end
fprintf(fid,'\\\\[-.15in] \n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);


















%% TABLE 3: AKM MODEL TO MODEL COMPARISON
% read data from stata
yearlist = { '1985_1986' , '2015_2016' , '1985_1988' , '2013_2016' , '1985_1992' , '2009_2016' , '1985_2000' , '2001_2016' , '1985_2016' } ;
clear akm n
for s = { 'KSS_f' 'KSS_fy' }
    for yy = yearlist 
%         load([directories.data 'results_MATLAB_AKM_' s{1} '_' yy{1} '.mat']) ;
%         load([directories.data1 '../AKM_' yy{1} '/results_MATLAB_AKM_' s{1} '_' yy{1} '.mat']) ;
        load([directories.data '/results_MATLAB_AKM_' s{1} '_' yy{1} '.mat']) ;
        ss = extractBefore(s{1},'_') ;
        sss = extractAfter(s{1},'_') ; 
        ssss = extractAfter(sss,'f') ;
        akm.(['pe' ssss ]).(['p' yy{1}]) = results.var_pe ;
        akm.(['pe' ssss ss]).(['p' yy{1}]) = results.var_pe_KSS ;
        akm.(['fe' ssss ]).(['p' yy{1}]) = results.(['var_' sss 'e']) ;
        akm.(['fe' ssss ss]).(['p' yy{1}]) = results.(['var_' sss 'e_KSS']) ;        
        n.firmyears.(['p' yy{1}]) = results.N_firm_years_pop ;
        n.firms.(['p' yy{1}]) = results.N_firms_pop ;
    end
end

% write table
% write table
fid = fopen([directories.tables 'Table4.tex'],'w');
fprintf(fid,'\\begin{tabular}{l c ');
len = 2 ;
for yy = yearlist(1:2:end-1)
    fprintf(fid,' cc c');
    len = len+3 ;
end
fprintf(fid,' ccc');
len = len+1 ;
fprintf(fid,'} \\hline \\hline \n');
%fprintf(fid,'&') ;
for yy = yearlist(1:2:end-1)
    s1 = extractBefore(yy{1},'_');
    s2 = extractAfter(yy{1},'_');
    fprintf(fid,'&& \\multicolumn{2}{c}{%2.0f years}',str2num(s2)-str2num(s1)+1);
end
for yy = yearlist(end)
    s1 = extractBefore(yy{1},'_');
    s2 = extractAfter(yy{1},'_');
    fprintf(fid,'&& %2.0f years',str2num(s2)-str2num(s1)+1);
end
fprintf(fid,'\\\\ \n');
xx = 3 ;
for yy = yearlist(1:2:end-1)
    fprintf(fid,'\\cline{%s-%s}',num2str(xx),num2str(xx+1));
    xx = xx+3 ;
end
fprintf(fid,'\\cline{%s-%s}',num2str(xx),num2str(xx));
fprintf(fid,'\\\\[-.15in] \n');
i = 1 ;
for yy = yearlist
    s1 = extractBefore(yy{1},'_');
    s2 = extractAfter(yy{1},'_');
    if mod(i,2)==1
        fprintf(fid,'&') ;
    end
    fprintf(fid,'& %s--%s',s1,s2) ;    
    i=i+1 ;
end
fprintf(fid,'\\\\ \\hline \n');
clear name
name.fe = '$Var(\hat{\psi}_{j})$' ;
name.fey = '$Var(\hat{\psi}_{jt})$' ;
name.fediff = 'Difference (\%)' ;
% output.cov_pe_fe = '$2*Cov(\hat{\alpha}_{i},\hat{\psi}_{jt})$' ;
for version = { '' , 'KSS' }
    fprintf(fid,'\\\\[-.1in] \n');
    if strcmp( version{1} , '' )
        fprintf(fid,'\\multicolumn{ %s }{l}{\\textit{Panel A. Plug-in estimates}} \\\\ \n',num2str(len));
    else
        fprintf(fid,'\\multicolumn{ %s }{l}{\\textit{Panel B. KSS bias-corrected estimates}} \\\\ \n',num2str(len));
    end
    i = 1 ;
    for s = { '' , 'y' , 'diff' }
        ss = { 'fe' } ;
        fprintf(fid,'\\hspace{.1in} %s &',name.([ss{1} s{1}])) ;
        j = 2 ;
        for yy = yearlist
            if strcmp(s{1},'diff')
                fprintf(fid,'& %4.1f ',100*(akm.([ss{1} 'y' version{1}]).(['p' yy{1}])./akm.([ss{1} version{1}]).(['p' yy{1}])-1)) ;
            else
                fprintf(fid,'& %4.3f ',akm.([ss{1} s{1} version{1}]).(['p' yy{1}])) ;
            end
            if mod(j,2)==1
                fprintf(fid,'&') ;
            end
            j = j+1 ;
        end
        fprintf(fid,'\\\\ \n');
        i = i+1 ;
    end
end
fprintf(fid,'\\\\[-.1in] \n');
fprintf(fid,'\\multicolumn{ %s }{l}{\\textit{Panel C. Descriptive statistics}} \\\\ \n',num2str(len)) ;
clear name
name.firmyears = 'Firm-years' ;
name.firms = 'Unique firms' ;
for s = fieldnames(name)'
    fprintf(fid,'\\hspace{.05in} %s &',name.(s{1})) ;
    j = 2 ;
    for yy = yearlist 
        fprintf(fid,'& %s ',comma_separator(n.(s{1}).(['p' yy{1}]))) ;
        if mod(j,2)==1
            fprintf(fid,'&') ;
        end
        j = j+1 ;
    end
    fprintf(fid,'\\\\ \n');
end
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);



% write table
fid = fopen([directories.tables 'Table4expanded.tex'],'w');
fprintf(fid,'\\begin{tabular}{l c ');
len = 2 ;
for yy = yearlist 
    fprintf(fid,'c');
    len = len+1 ;
end
fprintf(fid,'} \\hline \\hline \n');
fprintf(fid,'&');
for yy = yearlist 
    s1 = extractBefore(yy{1},'_');
    s2 = extractAfter(yy{1},'_');
    fprintf(fid,'& %2.0f years ',str2num(s2)-str2num(s1)+1);
end
fprintf(fid,'\\\\ \\hline \n');
clear name
name.pe = '$Var(\hat{\alpha}_{i})$' ;
name.pey = '$Var(\hat{\alpha}_{i})$' ;
name.fe = '$Var(\hat{\psi}_{j})$' ;
name.fey = '$Var(\hat{\psi}_{jt})$' ;
name.pediff = '\textit{Person}' ;
name.fediff = '\textit{Firm}' ;
% output.cov_pe_fe = '$2*Cov(\hat{\alpha}_{i},\hat{\psi}_{jt})$' ;
for version = { '' , 'KSS' }
    fprintf(fid,'\\\\[-.1in] \n');
    if strcmp( version{1} , '' )
        panel = 'A' ;
        fprintf(fid,'\\multicolumn{ %s }{l}{\\textbf{Panel A. Plug-in estimates}} \\\\ \n',num2str(len));
    else
        panel = 'B' ;
        fprintf(fid,'\\multicolumn{ %s }{l}{\\textbf{Panel B. KSS bias-corrected estimates}} \\\\ \n',num2str(len));
    end
    i = 1 ;
    for s = { '' , 'y' , 'diff' }
        if strcmp(s{1},'') 
            fprintf(fid,'\\multicolumn{ %s }{l}{\\hspace{.05in} Panel %s%1.0i. Firm fixed effect model} \\\\ \n',num2str(len),panel,i);
        elseif strcmp(s{1},'y') 
            fprintf(fid,'\\multicolumn{ %s }{l}{\\hspace{.05in} Panel %s%1.0i. Firm-year fixed effect model} \\\\ \n',num2str(len),panel,i);
        else
            fprintf(fid,'\\multicolumn{ %s }{l}{\\hspace{.05in} Panel %s%1.0i. Difference between firm-year FE model and firm FE model (\\%%)} \\\\ \n',num2str(len),panel,i);
        end
        for ss = { 'pe' , 'fe' }
            fprintf(fid,'\\hspace{.1in} %s &',name.([ss{1} s{1}])) ;
            for yy = yearlist
                if strcmp(s{1},'diff')
                    fprintf(fid,'& %4.3f ',100*(akm.([ss{1} 'y' version{1}]).(['p' yy{1}])./akm.([ss{1} version{1}]).(['p' yy{1}])-1)) ;
                else
                    fprintf(fid,'& %4.3f ',akm.([ss{1} s{1} version{1}]).(['p' yy{1}])) ;
                end
            end
            fprintf(fid,'\\\\ \n');
        end
        fprintf(fid,'\\\\[-.1in] \n');
        i = i+1 ;
    end
end
% fprintf(fid,'\\\\[.0in] \n');
fprintf(fid,'\\multicolumn{ %s }{l}{\\textbf{Panel C. Descriptive statistics}} \\\\ \n',num2str(len)) ;
clear name
name.firmyears = 'Firm-years' ;
name.firms = 'Unique firms' ;
for s = fieldnames(name)'
    fprintf(fid,'\\hspace{.05in} %s &',name.(s{1})) ;
    for yy = yearlist 
        fprintf(fid,'& %s ',comma_separator(n.(s{1}).(['p' yy{1}]))) ;
    end
    fprintf(fid,'\\\\ \n');
end
fprintf(fid,'\\hspace{.05in} Years &');
for yy = yearlist 
    s1 = extractBefore(yy{1},'_');
    s2 = extractAfter(yy{1},'_');
    fprintf(fid,'& %s--%s',s1,s2);
end
fprintf(fid,'\\\\[.05in] \n');
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);














%% TABLE 4: AUTOCOVARIANCE BY YEAR
% read data from stata 
data = csv2mat_numeric([directories.data2 'Table4.out']) ;
% write table
for var = { 'fye' }
    for version = { '' , '_w' }
        for spec = { '' , '_b' }
            fid = fopen([directories.tables 'Table5' var{1} spec{1} version{1} '.tex'],'w');
            fprintf(fid,'\\setlength{\\tabcolsep}{3.5pt} \n') ;
            fprintf(fid,'\\begin{tabular}{l ccccc ccccc ccccc ccccc ccccc ccccc cc}\n');
            fprintf(fid,'\\hline \\hline \n');
            fprintf(fid,'\\\\[-.1in] \n');
            for a=1985:2016
                fprintf(fid,' & %4.0f',a);
            end
            fprintf(fid,'\\\\ \n');
            fprintf(fid,'\\hline \n');
            fprintf(fid,'\\\\[-.08in] \n');
            fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel A. Autocorrelation}} \\\\ \n');
            for a = 1985:2016
                fprintf(fid,'%4.0f',a) ;
                for f = 1985:2016
                    if f<a
                        fprintf(fid,' & ') ;
                    else
                        fprintf(fid,' & %4.3f',data.([ var{1} '_cor' spec{1} version{1} num2str(a)])(f-a+1)) ;
                    end
                end
                fprintf(fid,'\\\\ \n');
            end
            fprintf(fid,'\\\\[-.08in] \n');
            fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel B. Autocovariance}} \\\\ \n');
            for a=1985:2016
                fprintf(fid,'%4.0f',a) ;
                for f = 1985:2016
                    if f<a
                        fprintf(fid,' & ') ;
                    else
                        fprintf(fid,' & %4.3f',data.([ var{1} '_cov' spec{1} version{1} num2str(a)])(f-a+1)) ;
                    end
                end
                fprintf(fid,'\\\\ \n');
            end    
            fprintf(fid,'\\hline\n');
            fprintf(fid,'\\end{tabular}');
            fclose(fid);

        end
    end
end

















%% TABLE 5: SHARE OF FIRMS BY YEAR
data = csv2mat_numeric([directories.data2 'Table5.out']) ;
% write table
fid = fopen([directories.tables 'TableA1.tex'],'w');
fprintf(fid,'\\setlength{\\tabcolsep}{3.5pt} \n') ;
fprintf(fid,'\\begin{tabular}{l ccccc ccccc ccccc ccccc ccccc ccccc cc}\n');
fprintf(fid,'\\hline \\hline \n');
for a=1985:2016
    fprintf(fid,' & %4.0f',a);
end
fprintf(fid,'\\\\ \n');
fprintf(fid,'\\hline \n');
for version = { '' , '_w' }
    fprintf(fid,'\\\\[-.08in] \n');
    if strcmp(version{1},'_w')
        fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel B. Weighted}} \\\\ \n');
    else
        fprintf(fid,'\\multicolumn{33}{c}{\\textit{Panel A. Unweighted}} \\\\ \n');    
    end
    for a = 1985:2016
        fprintf(fid,'%4.0f',a) ;
        for f = 1985:2016
            if f<a
                fprintf(fid,' & ') ;
            else
                fprintf(fid,' & %4.3f',data.(['active' version{1} '_' num2str(a) '_' num2str(f)]));
            end
        end
        fprintf(fid,'\\\\ \n');
    end
end
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);


















%% TABLE 6: SECOND STAGE REGRESSIONS
data = csv2mat_numeric([directories.data2 'Table6.out']) ;

% reorganize the data
for moment = fieldnames(data)'
    if contains(moment{1},'f_')
        field = extractBefore(string(moment{1}),'f_')+extractAfter(string(moment{1}),"f_") ;
        data.(field) = data.(moment{1}) ;
        data = rmfield(data,moment{1}) ;
    end
    if contains(moment{1},'_pw')
        field = extractBefore(string(moment{1}),'_pw')+extractAfter(string(moment{1}),"_pw") ;
        data.(field) = data.(moment{1}) ;
        data = rmfield(data,moment{1}) ;
    end
end
clear fields
fields = { 'test' } ;
for moment = fieldnames(data)'
    bef = char(extractBefore(string(moment{1}),'_')) ;
    if ~strcmp(bef,'joint') & ~contains(fields,bef)
        fields = { fields{1:end} , bef } ;
    end
end
fields = fields(2:end) ;

% write table
lags = { '' , '_d1' , '_d3' , '_d5' , '_d10' } ;
for weight = { 'w' , 'uw' }
    for winsor = [ 0 , 1 , 5 , 10 ] 
        fid = fopen([directories.tables 'Table6' weight{1} '_winsor' num2str(winsor) '.tex'],'w');
        fprintf(fid,'\\begin{tabular}{l c ccccc c ccccc}\n');
        fprintf(fid,'\\hline \\hline \n');
        fprintf(fid,'\\\\[-.1in] \n');
        fprintf(fid,' && \\multicolumn{5}{c}{Univariate} && \\multicolumn{5}{c}{Multivariate} \\\\ \n');
        fprintf(fid,'\\cline{2-7} \\cline{9-13} \\\\[-.1in] \n');
        fprintf(fid,' && Level & 1-year & 3-year & 5-year & 10-year && Level & 1-year & 3-year & 5-year & 10-year \\\\ \n');
        fprintf(fid,'\\hline \n');
        fprintf(fid,'\\\\[-.1in] \n');
        clear output
        output.size = 'Log number of workers' ;
        output.employees = 'Log number of workers' ;
        output.va = 'Log value added per worker' ;
        output.assets = 'Log capital' ;
        for moment = fields
            fprintf(fid,'%s & ',output.(moment{1})) ;
            for field = lags
                if strcmp(field{1},'')
                    fprintf(fid,'& %4.3f',data.([moment{1} '_' weight{1} field{1} '_beta'])) ;
                else
                    fprintf(fid,'& %4.3f',data.([moment{1} '_' weight{1} field{1} '_w' num2str(winsor) '_beta'])) ;
                end
            end
            fprintf(fid,'&') ;
            for field = lags
                if strcmp(field{1},'')
                    fprintf(fid,'& %4.3f',data.(['joint_' moment{1} '_' weight{1} field{1} '_beta'])) ;
                else
                    fprintf(fid,'& %4.3f',data.(['joint_' moment{1} '_' weight{1} field{1} '_w' num2str(winsor) '_beta'])) ;
                end
            end    
            fprintf(fid,'\\\\ \n');
            fprintf(fid,'&') ;
            for field = lags 
                if strcmp(field{1},'')
                    fprintf(fid,'& (%4.3f)',data.([moment{1} '_' weight{1} field{1} '_se'])) ;
                else
                    fprintf(fid,'& (%4.3f)',data.([moment{1} '_' weight{1} field{1} '_w' num2str(winsor) '_se'])) ;
                end
            end
            fprintf(fid,'&') ;
            for field = lags 
                if strcmp(field{1},'')
                    fprintf(fid,'& (%4.3f)',data.(['joint_' moment{1} '_' weight{1} field{1} '_se'])) ;
                else
                    fprintf(fid,'& (%4.3f)',data.(['joint_' moment{1} '_' weight{1} field{1} '_w' num2str(winsor) '_se'])) ;
                end
            end
            fprintf(fid,'\\\\[.05in] \n');
        end
        fprintf(fid,'\\\\[0in] \n');

        fprintf(fid,'$R^2$ & ') ;
        for field = lags
            fprintf(fid,'&') ;
        end
        fprintf(fid,'&') ;
        for field = lags
            if strcmp(field{1},'')
                fprintf(fid,'& %4.3f',data.(['joint_' weight{1} field{1} '_r2'])) ;
            else
                fprintf(fid,'& %4.3f',data.(['joint_' weight{1} field{1} '_w' num2str(winsor) '_r2'])) ;
            end
        end
        fprintf(fid,'\\\\[.05in] \n');
        clear output
        output.va = 'Firm-years' ;
        moment = fields(1) ;
        fprintf(fid,'Firm-years & ') ;
        for field = lags
            if strcmp(field{1},'')
                fprintf(fid,'& %s',comma_separator(data.([moment{1} '_' weight{1} field{1} '_n']))) ;
            else
                fprintf(fid,'& %s',comma_separator(data.([moment{1} '_' weight{1} field{1} '_w' num2str(winsor) '_n']))) ;
            end
        end
        fprintf(fid,'&') ;
        for field = lags
            if strcmp(field{1},'')
                fprintf(fid,'& %s',comma_separator(data.(['joint_' weight{1} field{1} '_n']))) ;
            else
                fprintf(fid,'& %s',comma_separator(data.(['joint_' weight{1} field{1} '_w' num2str(winsor) '_n']))) ;
            end
        end
        fprintf(fid,'\\\\[.05in] \n');
        fprintf(fid,'\\hline\n');
        fprintf(fid,'\\end{tabular}');
        fclose(fid);
    end
end

% OLD VERSION
% data = csv2mat_numeric([directories.data2 'Table6.out']) ;
% 
% % write table
% fid = fopen([directories.tables 'Table6.tex'],'w');
% fprintf(fid,'\\begin{tabular}{l c cccc c cccc}\n');
% fprintf(fid,'\\hline \\hline \n');
% fprintf(fid,' && \\multicolumn{4}{c}{Unweighted} && \\multicolumn{4}{c}{Weighted} \\\\ \n');
% fprintf(fid,'\\cline{2-6} \\cline{8-11} \\\\[-.15in] \n');
% fprintf(fid,' && \\multicolumn{2}{c}{Univariate} & \\multicolumn{2}{c}{Multivariate} && \\multicolumn{2}{c}{Univariate} & \\multicolumn{2}{c}{Multivariate} \\\\ \n');
% fprintf(fid,'\\hline \n');
% clear output
% output.f_age = 'Firm age' ;
% output.f_size = 'Firm size' ;
% output.rev_pw = 'Sales per worker' ;
% output.va_pw = 'Va.p.w.' ;
% output.assets = 'Assets' ;
% output.assets_pw = 'Assets p.w.' ;
% output.debt_pw = 'Debt per worker' ;
% output.equity_pw = 'Equity per worker' ;
% output.inv_pw = 'Investment per worker' ;
% for str = fieldnames(data)'
%     if ~strcmp(str{1},'firms') && ~strcmp(str{1},'firmyears') && contains(str{1},'fe_uw_sd')
%         temp = extractBefore(string(str{1}),"_") ;
%         if strcmp(temp,'f')
%             moment = extractBefore(extractAfter(string(str{1}),"f_"),"_") ;
%             moment = "f_" + moment ;
%         elseif ~contains(str{1},"_pw")
%             moment = extractBefore(string(str{1}),"_");
%         else
%             moment = extractBefore(string(str{1}),"_")+"_pw" ;
%         end
%         fprintf(fid,'%s && %4.3f & %4.3f & %4.3f & %4.3f && %4.3f & %4.3f & %4.3f & %4.3f \\\\ \n', ...
%                                                  output.(moment{1}), ...
%                                                  data.([moment{1} '_uw']),data.([moment{1} '_fe_uw']), ...
%                                                  data.([moment{1} '_uw_joint']),data.([moment{1} '_fe_uw_joint']), ...
%                                                  data.([moment{1} '_w']),data.([moment{1} '_fe_w']), ...
%                                                  data.([moment{1} '_w_joint']),data.([moment{1} '_fe_w_joint']) ) ;
%         fprintf(fid,'&& (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f) && (%4.3f) & (%4.3f) & (%4.3f) & (%4.3f) \\\\ \n', ...
%                                                  data.([moment{1} '_uw_sd']),data.([moment{1} '_fe_uw_sd']), ...
%                                                  data.([moment{1} '_uw_joint_sd']),data.([moment{1} '_fe_uw_joint_sd']), ...
%                                                  data.([moment{1} '_w_sd']),data.([moment{1} '_fe_w_sd']), ...
%                                                  data.([moment{1} '_w_joint_sd']),data.([moment{1} '_fe_w_joint_sd']) ) ;
%         fprintf(fid,'\\\\[-.1in] \n');
%     end
% end
% 
% fprintf(fid,'\\\\[.0in] \n');
% clear output
% output.r2 = '$R^2$' ;
% output.r2within = '$R^2$ (within)' ;
% for moment = fieldnames(output)'
%     if contains(moment{1},'within')
%         fprintf(fid,'%s && & & & %4.3f && & & & %4.3f \\\\ \n',output.(moment{1}),data.([moment{1} '_joint_fe_uw']),data.([moment{1} '_joint_fe_w'])) ;
%     else
%         fprintf(fid,'%s && & & %4.3f & %4.3f && & & %4.3f & %4.3f \\\\ \n',output.(moment{1}),data.([moment{1} '_joint_uw']),data.([moment{1} '_joint_fe_uw']),data.([moment{1} '_joint_w']),data.([moment{1} '_joint_fe_w'])) ;
%     end 
% end
% fprintf(fid,'\\\\[.0in] \n');
% clear output
% output.firmyears = 'Firm-years' ;
% output.firms = 'Firms' ;
% for moment = fieldnames(output)'
%     temp = comma_separator(data.(moment{1})) ;
%     fprintf(fid,'%s && %s & %s & %s & %s && %s & %s & %s & %s \\\\ \n',output.(moment{1}),temp,temp,temp,temp,temp,temp,temp,temp) ;
% end
% fprintf(fid,'\\\\[.0in] \n');
% fprintf(fid,'Firm FEs && No & Yes & No & Yes && No & Yes & No & Yes \\\\ \n') ;
% fprintf(fid,'\\\\[-.15in] \n');
% fprintf(fid,'\\hline\n');
% fprintf(fid,'\\end{tabular}');
% fclose(fid);












%% TABLE 7: AUTOCOVARIANCE BY FIRM AGE
% read data from stata 
data = csv2mat_numeric([directories.data2 'Table7.out']) ;

% write table
for version = { '' , '_w' }
    fid = fopen([directories.tables 'Table7' version{1} '.tex'],'w');
    fprintf(fid,'\\begin{tabular}{l ccccc ccccc ccccc ccccc}\n');
    fprintf(fid,'\\hline \\hline \n');
    fprintf(fid,'\\\\[-.1in] \n');
    fprintf(fid,'\\multicolumn{21}{c}{\\textit{Panel A. Unbalanced panel}} \\\\ \n');
    for a=0:19
        fprintf(fid,' & %2.0f',a);
    end
    fprintf(fid,'\\\\ \n');
    fprintf(fid,'\\hline \n');
    for a = 0:19
        fprintf(fid,'%2.0f',a) ;
        for f = 1:20
            if f<=a
                fprintf(fid,' & ') ;
            else
                fprintf(fid,' & %4.3f',data.(['autocorr' version{1} num2str(a)])(f-a)) ;
            end
        end
        fprintf(fid,'\\\\ \n');
    end
    fprintf(fid,'\\\\[-.1in] \n');
    fprintf(fid,'\\multicolumn{21}{c}{\\textit{Panel B. Balanced panel}} \\\\ \n');
    for a = 0:19
        fprintf(fid,'%2.0f',a) ;
        for f = 1:20
            if f<=a
                fprintf(fid,' & ') ;
            else
                fprintf(fid,' & %4.3f',data.(['autocorr_balanced' version{1} num2str(a)])(f-a)) ;
            end
        end
        fprintf(fid,'\\\\ \n');
    end    
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\end{tabular}');
    fclose(fid);
    
end


























%% TABLE 5: STABILITY OF FIRM YEAR FES ACROSS MODELS
% read data from stata 
data = csv2mat_numeric([directories.data3 '7_comp_short_full.csv']) ;
fid = fopen([directories.tables 'Table8.tex'],'w');
fprintf(fid,'\\begin{tabular}{l c c}\n');
fprintf(fid,'\\hline \\hline \n');
fprintf(fid,'\\\\[-.1in] \n');
fprintf(fid,' & Levels & Differences');
fprintf(fid,'\\\\ \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\\\[-.08in] \n');
fprintf(fid,'$\\widehat{\\eta}$ & %4.3f & %4.3f \\\\ \n',data.b_level_year,data.b_diff_year);
fprintf(fid,' & (%4.3f) & (%4.3f) \\\\ \n',data.se_level_year,data.se_diff_year);
fprintf(fid,'Correlation & %4.3f & %4.3f \\\\ \n',data.rho_level_year,data.rho_diff_year);
fprintf(fid,'\\hline\n');
fprintf(fid,'\\end{tabular}');
fclose(fid);
