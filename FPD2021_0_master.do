********************************************************************************
* DESCRIPTION: Master file for Engbom, Moser, and Sauermann (2021).
*
* INPUTS:      User must update these items in the current file:
*                 - which sections to run
*                 - variable names
*                 - parameters
*                 - directories
*              User must create the following local directories:
*                 - DIR_MAIN
*                 - DIR_CODE
*                 - DIR_SUB
*              User must locate the following (executable or application) file:
*                 - FILE_MATLAB
*              User must update the following items in the MATLAB files
*              "FUN_CONNECTED.m", "FUN_AKM.m", "FUN_AKM_KSS.m":
*                 - input path variable DIR_INPUT
*                 - input path variable DIR_OUTPUT
*                 - input path variable DIR_LOG
*              Packages required to (efficiently) run programs:
*                 - confirmdir (Windows only)
*                 - reghdfe
*                 - gtools
*
* OUTPUTS:     See individual Stata and MATLAB code files for list of outputs.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 28, 2021.
********************************************************************************


********************************************************************************
* INITIAL HOUSEKEEPING
********************************************************************************
* clear memory and set system parameters -- NOTE: DO NOT CHANGE!
set more off
clear all
macro drop _all
cap log close _all
set type double
set excelxlsxlargefile on
set graphics off
set varabbrev off
set rmsg on
set matsize 11000
set linesize 100
set seed 12345
set maxvar 120000


********************************************************************************
* MACROS
********************************************************************************
* sections to run -- NOTE: NEED TO CHANGE!
global FPD2021_1_test_data = 0
global FPD2021_2_selection = 0
global FPD2021_3_connected = 0
global FPD2021_4_akm = 0
global FPD2021_5_panels = 0
global FPD2021_6_tables = 0
global FPD2021_7_comparison = 0
global FPD2021_8_time_series = 0
global FPD2021_9_lifecycle = 0
global FPD2021_X_additional_figures = 0

* linked employer-employee data variable names -- NOTE: NEED TO CHANGE! In SWE: persid (string), year (1985-2015), educ (1-7), age (18-64), exper (1-48), empid (string), estabsnr (string), earnings.
global year = "year" // year
global id_pers = "persid" // person ID
global id_firm = "empid" // firm ID
global id_est = "estabsnr" // establishment ID
global gender = "gender" // gender (1 = male, 2 = female)
global age = "age" // worker age (years)
global edu = "educ" // education level
global emp_dur = "emp_dur" // employment duration during current year (hours, months, or quarters)
global earn = "earnings" // monthly earnings (log real SEK)
global wage = "wage" // monthly wage (log real SEK)
global hours = "hours" // contractual weekly work time (hours)
global occ = "occ" // occupation code
global exp = "exper" // labor market experience since graduation
global tenure = "tenure" // tenure with current employer (years)

* firm financials data variable names -- NOTE: NEED TO CHANGE! In SWE: empid, estabsnr, fsize_estab, fsize, incorp, public, year, sector, ownerid, sales, va, ebit, profit, payroll, employees, assets, assets_st, cash, assets_lt, equipment, structures, debt, equity, investment, va_estab, payroll_estab, employees_estab, IDChangeFuture, IDChangePrevi~s, firm, estab, yob_firm, va_pw, va_pw_l, assets_pw, assets_pw_l, payroll_pw, payroll_pw_l, tfp_pw_l, tfp, tfp_l, va_pw_l_resid, va_pw_resid, va_resid, va_l_resid, location_estab
global f_size = "f_size" // firm size from RAMS (number of employees)
global f_size_fek = "employees" // firm size from FEK (number of employees)
global f_ind = "sector" // firm industry code
global f_yob "yob_firm" // firm year of birth
global f_age "f_age" // firm age (years)
global f_rev = "sales" // firm-level revenues (SEK)
global f_va = "va" // firm-level value added (SEK)
global f_ebit = "ebit" // firm-level EBIT (SEK)
global f_profit = "profit" // firm-level profits (SEK)
global f_assets = "assets" // firm-level assets (SEK)
global f_assets_lt = "assets_lt" // firm-level long-term assets (SEK)
global f_equip = "equipment" // firm-level equipment (SEK)
global f_struct = "structures" // firm-level structures (SEK)
global f_cash = "cash" // firm-level cash (SEK)
global f_debt = "debt" // firm-level debt (SEK)
global f_equity = "equity" // firm-level equity (SEK)
global f_inv = "investment" // firm-level investment (SEK)

* user-created variable names -- NOTE: DO NOT CHANGE!
global id_unique = "id_unique" // unique job ID
global pe = "pe" // AKM person FE
global f_fe = "fe" // AKM firm FE
global fy_fe = "fye" // AKM firm-year FE
global xb_year = "xb_year" // AKM (education-specific) year FEs
global xb_age = "xb_age" // AKM (education-specific) worker age FEs
global xb_emp_age = "xb_e_age" // AKM employer age FEs
global resid = "resid" // AKM residual

* user-set parameters -- NOTE: NEED TO CHANGE!
global test = 0 // 0 = real run on real data; 1 = test run on made up data
global sample_read = 0 // 0 = read full data, 1 = read sample data
global year_min_list = "1985   1985 2001   1985 1993 2001 2009   1985 1989 1993 1997 2001 2005 2009 2013   1985 1987 1989 1991 1993 1995 1997 1999 2001 2003 2005 2007 2009 2011 2013 2015" // list of start years of estimation period for CONNECTED and AKM (e.g., "1985   1985 2001   1985 1993 2001 2009   1985 1989 1993 1997 2001 2005 2009 2013   1985 1987 1989 1991 1993 1995 1997 1999 2001 2003 2005 2007 2009 2011 2013 2015") -- leave empty in order to run just once!
global year_max_list = "2016   2000 2016   1992 2000 2008 2016   1988 1992 1996 2000 2004 2008 2012 2016   1986 1988 1990 1992 1994 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016" // list of  end  years of estimation period for CONNECTED and AKM (e.g., "2016   2000 2016   1992 2000 2008 2016   1988 1992 1996 2000 2004 2008 2012 2016   1986 1988 1990 1992 1994 1996 1998 2000 2002 2004 2006 2008 2010 2012 2014 2016") -- leave empty in order to run just once!
global year_min = 1985 // start year of estimation period (e.g., 1985)
global year_max = 2016 // end year of estimation period (e.g., 2015)
global year_min_data = 1985 // start year of data, regardless of estimation period (e.g., 1985)
global year_max_data = 2016 // end year of data, regardless of estimation period (e.g., 2015)
global parallel_max = 1 // 0 = compute maximum number of parallel workers to use in MATLAB as = round(3 + 34/(year_max - year_min + 1)); > 0 = maximum number of parallel workers to use in MATLAB
global gtools = "g" // "g" = use gtools package for commands such as -collapse-, -egen-, etc.; "" = use Stata-native commands
global n_pool = 2 // number of pools to use in reghdfe -- higher number of pools is faster but requires more RAM
global gender_min = 1 // minimum gender (1 = men, 2 = women)
global gender_max = 2 // maximum gender (1 = men, 2 = women)
global age_min = 20 // minimum worker age
global age_max = 59 // maximum worker age
global demean_income_year = 1 // 0 = use raw income data; 1 = demean income by (demographic-specific) year
global demean_income_gender = 0 // 0 = do not demean income by gender interaction; 1 = also demean income by gender interaction
global demean_income_edu = 0 // 0 = do not demean income by education interaction; 1 = also demean income by education interaction
global akm_year_dummies = 0 // 0 = do not include year dummies; 1 = include year dummies
global akm_age_poly_order = 1 // 0 = do not include age terms in AKM regression; 1 = age dummies with income-age profile restricted to be flat from age ${akm_age_flat_min} to ${akm_age_flat_max}; 2 / 3 / etc. = include 2nd order term / 2nd and 3rd order terms / etc.
global akm_age_flat_min = 50 // minimum worker age for which income-age profile is restricted to be flat (only relevant if ${akm_age_poly_order} == 1)
global akm_age_flat_max = 51 // maximum worker age for which income-age profile is restricted to be flat (only relevant if ${akm_age_poly_order} == 1)
global akm_age_norm = 50 // age around which to normalize higher-order age polynomial terms in AKM estimation (relevant only if ${akm_age_poly_order} >= 2; coincides with where the age profile is assumed to be flat for interpretation of worker FEs and year FEs)
global akm_gender_inter = 1 // 0 = no gender interactions; 1 = interact gender with time trends and age profiles
global akm_edu_inter = 1 // 0 = no education interactions; 1 = interact education with time trends and age profiles
global akm_hours = 0 // 0 = do not control for contractual work hours; 1 = control for contractual work hours
global akm_occ = 0 // 0 = do not control for occupational categories; 1 = control for occupational categories
global akm_occ_coarsen = 0 // 0 = use original occupation categories as controls; 1 = use coarsened occupation categories (/10) as controls
global akm_tenure = 0 // 0 = do not control for years of tenure; 1 = control for years of tenure
global akm_tenure_top_coding = 0 // 0 = do not top-code tenure; 1 = top-code tenure at ${akm_tenure_max}
global drop_singletons_pers = 1 // AKM first stage regression: drop singleton workers?
global drop_singletons_firm = 0 // AKM first stage regression: drop singleton firms?
global drop_singletons_firm_year = 0 // AKM first stage regression: drop singleton firm-years?
global min_size_emp = 1 // AKM first stage regression: minimum firm size used for AKM estimation
global min_switchers = 1 // AKM first stage regression: minimum number of workers switching employers
global inc_concept_list = "${earn}" // list of worker income concepts ("${earn}", "${wage}", or "${earn} ${wage}")
global inc_default = "${earn}" // default worker income concept ("${earn}" or "${wage}")
global lower_perc = 5 // trim observations below this earnings percentile ("" = do not trim; "0.1" = trim lower 1/1000; "1" = trim lower 1/100; "5" = trim lower 5/100; "10" = trim lower 10/100)
global drop_below_mw = 1 // 0 = do not drop observations based on earnings relative to minimum wage, 1 = drop observations below the minimum wage (BRA)
global ${earn}_min = 0 // minimum worker earnings threshold (in levels, >= 0)
global ${earn}_max = . // maximum worker earnings threshold (in levels, <= .)
global ${wage}_min = 0 // minimum worker wage threshold (in levels, >= 0)
global ${wage}_max = . // maximum worker wage threshold (in levels, <= .)
global inc_growth_min = . // minimum income growth rate (in %) allowed for observation not to be dropped (. = do not drop any)
global inc_growth_max = . // maximum income growth rate (in %) allowed for observation not to be dropped (. = do not drop any)
global AKM_threshold_list = "1" // list of firm size thresholds for AKM estimation (1, 5, 10, 15, 25, 50, or 100)
global save_space = 1 // 0 = do not save space by first saving and later loading temporarily redundant variables in MATLAB; 1 = save space by first saving and later loading temporarily redundant variables in MATLAB

* test data parameters -- NOTE: DO NOT CHANGE!
global n_obs_year = 10000 // number of workers = number of observations per year
global n_firms = 100 // number of firms
global pe_mean = 0 // mean of person FEs in earnings equation
global pe_sd = 0.5^0.5 // std. dev. of person FEs in earnings equation
global f_fe_mean = 0 // mean of firm FEs in earnings equation
global f_fe_sd = 0.15^0.5 // std. dev. of firm FEs in earnings equation
global fy_fe_mean = 0 // mean of firm-year FEs in earnings equation
global fy_fe_sd = 0.05^0.5 // std. dev. of firm-year FEs in earnings equation
global eps_mean = 0 // mean of residual in earnings equation
global eps_sd = 0.05^0.5 // std. dev. of residual in earnings equation
global test_switcher_share = 0.1 // share of all simulated workers who ever switch employers during entire period (0 = none, 1 = all)

* automatically-set parameters -- NOTE: DO NOT CHANGE!
global dem_inter = (${akm_gender_inter} | ${akm_edu_inter})
if $demean_income_year | $demean_income_gender | $demean_income_edu global demeaned = "_demeaned"
else global demeaned = ""
if $test global test_str = "test_"
else global test_str = ""

if $test local n_bins_default = 20
else local n_bins_default = 50
global n_bins_comparison = `n_bins_default' // number of bins for pay concept comparison (default for non-testing run = 50)
global n_bins_changes = `n_bins_default' // number of bins for firm pay changes analysis (default for non-testing run = 50)
global n_bins_mob = `n_bins_default' // number of bins for firm pay mobility analysis (default for non-testing run = 50)
global n_bins_entry = `n_bins_default' // number of bins for analysis of firm pay heterogeneity at firm entry (default for non-testing run = 50)
global n_bins_cohort = 20 // number of bins for analysis of firm pay heterogeneity after firm entry (default for non-testing run = 20)
global n_bins_exit = `n_bins_default' // number of bins for analysis of firm pay heterogeneity around firm exit (default for non-testing run = 50)
global n_bins = `n_bins_default' // number of bins for AKM second stage (default for non-testing run = 50)



********************************************************************************
* DIRECTORIES
********************************************************************************
* use either Stata-native command -cap confirm file- (MacOSX, Unix) or user-written package -confirmdir- (Windows) to set the right command to check for directories
if inlist("`c(os)'", "MacOSX", "Unix") global cmd_confirm_dir = "cap confirm file"
else if "`c(os)'" == "Windows" global cmd_confirm_dir = "confirmdir"

* determine user
${cmd_confirm_dir} "P:/2019/107/"
if !_rc global user = "SWE" // "Nik" for Nik's local computer; "Chris_office" for Chris' office iMac Pro or MacBook; "Chris_home" for Chris' home iMac Pro; "SWE" for Uppsala server, "Grid" for Columbia GSB server
${cmd_confirm_dir} "/Users/cm3594/"
if !_rc global user = "Chris_office"
${cmd_confirm_dir} "/Users/economoser/"
if !_rc global user = "Chris_home"
${cmd_confirm_dir} "/Users/cmoser/"
if !_rc global user = "Chris_laptop"
${cmd_confirm_dir} "/Users/niklasengbom/"
if !_rc {
	global user = "Nik_home"
// 	global user = "Nik_office"
}
${cmd_confirm_dir} "/shared/share_cmoser/"
if !_rc global user = "Grid"

* set local directories and file names -- NOTE: NEED TO CHANGE!
if "${user}" == "SWE" { // when run on Swedish data server
	adopath + "P:/2018/144/programs"
	global DIR_MAIN = "P:/2019/107"
	global DIR_CODE = "${DIR_MAIN}/1 Code"
	global DIR_SUB = "${DIR_MAIN}/2 Data"
	global DIR_TEMP = "${DIR_SUB}"
	global DIR_DATA = "${DIR_TEMP}"
	global DIR_OUTPUT = "${DIR_TEMP}"
	global DIR_LOG = "${DIR_MAIN}/3 Logs"
	global DIR_EXPORT = "${DIR_MAIN}/4 Export"	
	global FILE_MATLAB = "C:/Program Files/MATLAB/R2020a/bin/matlab.exe"
	global RawDataDir = "P:/2018/144/1 Data/3 Clean"
}
else if inlist("${user}", "Chris_office", "Chris_home", "Chris_laptop", "Nik_home", "Nik_office", "Grid") { // when run on local machine by Chris or Nik
	if "${user}" == "Chris_office" global DIR_MAIN = "/Users/cm3594/Dropbox (CBS)/AKM Conference"
	else if "${user}" == "Chris_home" global DIR_MAIN = "/Users/economoser/Dropbox (CBS)/AKM Conference"
	else if "${user}" == "Chris_laptop" global DIR_MAIN = "/Users/cmoser/Dropbox (CBS)/AKM Conference"
	else if "${user}" == "Nik_home" global DIR_MAIN = "/Users/niklasengbom/Dropbox/AKM"
	else if "${user}" == "Nik_office" global DIR_MAIN = "/Users/niklasengbom/Dropbox/AKM"
	else if "${user}" == "Grid" global DIR_MAIN = "/shared/share_cmoser/20_AKM_Conference"
	global DIR_CODE = "${DIR_MAIN}/3_code/_work_in_progress"
	global DIR_SUB = "${DIR_MAIN}/5_results"
	global DIR_TEMP = "${DIR_SUB}/temp"
	global DIR_DATA = "${DIR_SUB}/data"
	global DIR_OUTPUT = "${DIR_SUB}/output"
	global DIR_LOG = "${DIR_SUB}/log"
	global DIR_EXPORT = "${DIR_SUB}/export"	
	if "${user}" == "Chris_office" {
		foreach f in ///
			"/Applications/MATLAB_R2019b.app/bin/matlab" ///
			"/Applications/MATLAB_R2020a.app/bin/matlab" ///
			"/Applications/MATLAB_R2020b.app/bin/matlab" ///
			"/Applications/MATLAB_R2021a.app/bin/matlab" ///
			"/Applications/MATLAB_R2021b.app/bin/matlab" ///
			{
			cap confirm file "`f'"
			if !_rc global FILE_MATLAB = "`f'"
		}
	}
	else if "${user}" == "Chris_home" global FILE_MATLAB = "/Applications/MATLAB_R2019b.app/bin/matlab"
	else if "${user}" == "Chris_laptop" global FILE_MATLAB = "/Applications/MATLAB_R2021b.app/bin/matlab"
	else if "${user}" == "Nik_home" global FILE_MATLAB = "/Applications/MATLAB_R2020b.app/bin/matlab"
	else if "${user}" == "Nik_office" global FILE_MATLAB = "/Applications/MATLAB_R2020a.app/bin/matlab"
	else if "${user}" == "Grid" global FILE_MATLAB = "/apps/MATLAB/R2021a/bin/matlab"
	if inlist("${user}", "Chris_office", "Chris_home") global RawDataDir = "/Users/cm3594/Data/RAIS/3_processed"
	else if "${user}" == "Chris_laptop" global RawDataDir = "/Users/cmoser/Data/RAIS/3_processed"
	else if "${user}" == "Grid" global RawDataDir = "/shared/share_cmoser/1_data/RAIS/3_processed"
}
cd "${DIR_SUB}"

* set server directories for server -- these need to be changed only if run on SWE data server (i.e., if "${user}" == "SWE")!
foreach dir in ///
	"${DIR_MAIN}/4 Export" ///
	"${DIR_TEMP}/_processed" ///
	{
	${cmd_confirm_dir} "`dir'"
	if !_rc global ExportDir = "`dir'"
}

* check that all directories and files exist -- NOTE: DO NOT CHANGE!
foreach dir_file in ///
	DIR_MAIN ///
	DIR_CODE ///
	DIR_SUB ///
	DIR_TEMP ///
	DIR_DATA ///
	DIR_OUTPUT ///
	DIR_LOG ///
	DIR_EXPORT ///
	FILE_MATLAB ///
	RawDataDir ///
	ExportDir ///
	{
	if "`dir_file'" == "FILE_MATLAB" cap confirm file "${FILE_MATLAB}"
	else ${cmd_confirm_dir} "${`dir_file'}"
	if _rc & inlist("`dir_file'", "DIR_MAIN") disp as error "USER ERROR: Directory does not exist -- breaking execution (`dir_file' = ${`dir_file'})."
	else if _rc & "`dir_file'" == "FILE_MATLAB" disp as error "USER ERROR: MATLAB file does not exist -- breaking execution (`dir_file' = ${`dir_file'})."
	else if _rc & inlist("`dir_file'", "DIR_CODE", "DIR_SUB", "DIR_TEMP", "DIR_DATA", "DIR_OUTPUT", "DIR_LOG", "DIR_EXPORT") {
		disp as error "USER WARNING: Directory does not exist -- creating it automatically (`dir_file' = ${`dir_file'})."
		!mkdir "${`dir_file'}"
	}
	if _rc & inlist("`dir_file'", "DIR_MAIN", "FILE_MATLAB") error 1
}
foreach subdir in ///
	csv ///
	eps ///
	pdf ///
	tex ///
	{
	${cmd_confirm_dir} "${DIR_EXPORT}/`subdir'"
	if _rc {
		disp as error "USER WARNING: Directory does not exist -- creating it automatically (${DIR_EXPORT}/`subdir')."
		!mkdir "${DIR_EXPORT}/`subdir'"
	}
}

* rebuild list of function libraries
mata: mata mlib index


********************************************************************************
* PROGRAMS
********************************************************************************
* user-written programs -- NOTE: DO NOT CHANGE!
program drop _all
program prog_comp_desc_sum_save // input: `1' = file name; output = compressed data plus screen output from commands -desc- and -sum-
	compress
	desc, fullnames
	sum, sep(0)
	save "`1'", replace
end
program graph_export_eps_pdf // input: `1' = file name; output = figures in EPS and PDF format
	graph export "${DIR_EXPORT}/eps/`1'.eps", replace
	graph export "${DIR_EXPORT}/pdf/`1'.pdf", replace
end
program graph_export_eps_pdf_png // input: `1' = file name; output = figures in EPS, PDF, and PNG format
	graph export "${DIR_EXPORT}/eps/`1'.eps", replace
	graph export "${DIR_EXPORT}/pdf/`1'.pdf", replace
	graph export "${DIR_EXPORT}/png/`1'.png", replace
end


********************************************************************************
* EXECUTE CODE
********************************************************************************
* run code sections -- NOTE: DO NOT CHANGE!
global sections: all globals
global FDP2021_sections = ""
foreach section of global sections {
	if substr("`section'", 1, 8) == "FPD2021_" global FDP2021_sections = "`section' ${FDP2021_sections}"
}
local section_connected = ""
local section_akm = ""
foreach section of global FDP2021_sections {
	if strpos(strlower("`section'"), "connected") local section_connected = "`section'"
	if strpos(strlower("`section'"), "akm") local section_akm = "`section'"
}
foreach section of global FDP2021_sections {
	if $`section' {
		local section_short = subinstr("`section'", "FPD2021_", "", .)
		log using "${DIR_LOG}/log_Stata_`section_short'_${year_min}_${year_max}", text name(`section') replace
		disp _newline(25)
		disp "********************************************************************************"
		disp "* STARTED `section'.do (${S_DATE}, ${S_TIME})"
		disp "********************************************************************************"
		timer on 1
		if "`section'" == "`section_connected'" & "${year_min_list}" != "" & "${year_max_list}" != "" { // for CONNECTED, loop through list of start/end years and call AKM from within
			local n_min: word count ${year_min_list}
			local n_max: word count ${year_max_list}
			cap assert `n_min' == `n_max'
			if _rc {
				disp as error "USER ERROR: Lists of minimum (\${year_min_list}) and maximum (\${year_max_list}) years are not of same length."
				error 1
			}
			forval n = 1/`n_min' {
				local year_min: word `n' of ${year_min_list}
				local year_max: word `n' of ${year_max_list}
				cap assert ${year_min_data} <= `year_min'
				display _newline(25)
				disp "*** CONNECTED/AKM for years `year_min'-`year_max'"
				if _rc {
					disp as error "USER ERROR: Must specify min. data year <= min. estimation year (year_min_data <= year_min)."
					error 1
				}
				if $`section_connected' do "${DIR_CODE}/`section_connected'.do" `year_min' `year_max'
				if $`section_akm' do "${DIR_CODE}/`section_akm'.do" `year_min' `year_max'
			}
		}
		else if "`section'" != "`section_akm'" | "${year_min_list}" == "" | "${year_max_list}" == "" { // if looping through list of start/end years, then skip AKM; for all other files, just run once with fixed start and end years
			cap assert ${year_min_data} <= ${year_min}
			if _rc {
				disp as error "USER ERROR: Must specify min. data year <= min. estimation year (year_min_data <= year_min)."
				error 1
			}
			do "${DIR_CODE}/`section'.do"
		}
		timer off 1
		timer list 1
		disp "********************************************************************************"
		disp "* FINISHED `section'.do (${S_DATE}, ${S_TIME}, RUNTIME = `=r(t1)' SECONDS)"
		disp "********************************************************************************"
		log close `section'
	}
	else if "$`section'" == "0" {
		disp _newline(10)
		disp "********************************************************************************"
		disp "* SKIPPING `section'.do (${S_DATE}, ${S_TIME})"
		disp "********************************************************************************"
	}
	else { // e.g., if "$`section'" == ""
		disp _newline(10)
		disp as error "********************************************************************************"
		disp as error "* COULD NOT FIND `section'.do (${S_DATE}, ${S_TIME})"
		disp as error "********************************************************************************"
	}
}

	
********************************************************************************
* FINAL HOUSEKEEPING
********************************************************************************
clear
cap log close _all
