********************************************************************************
* DESCRIPTION: Prepares data, calls MATLAB file FUN_CONNECTED.m to estimate
*              (regular and leave-one-out) connected sets, then restricts data
*              to most restrictive (leave-one-out firm-year FE) connected set
*              and outsheets restricted data in .csv format for AKM estimation.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 27, 2021.
********************************************************************************


*** read passed arguments
if "`1'" != "" & "`2'" != "" {
	global year_min_local = `1'
	global year_max_local = `2'
}
else {
	global year_min_local = ${year_min}
	global year_max_local = ${year_max}
}


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		* load data
		local vars_list = "${id_pers}* ${year} ${id_firm}* ${gender} ${edu} `inc_concept'"
// 		if (${akm_gender_inter} & (${akm_year_dummies} | ${akm_age_poly_order})) local vars_list = "`vars_list' ${gender}"
// 		if (${akm_edu_inter} & (${akm_year_dummies} | ${akm_age_poly_order})) local vars_list = "`vars_list' ${edu}"
		if ${akm_age_poly_order} local vars_list = "`vars_list' ${age}"
		if ${akm_hours} local vars_list = "`vars_list' ${hours}"
		if ${akm_occ} local vars_list = "`vars_list' ${occ}"
		if ${akm_tenure} & (!$akm_tenure_top_coding | ${year_min_local} > ${year_min_data} | inlist("${user}", "Chris_office", "Chris_home", "Grid")) local vars_list = "`vars_list' ${tenure}"
		use `vars_list' if inrange(${year}, ${year_min_local}, ${year_max_local}) using "${DIR_TEMP}/selection_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta", clear
		
		
		*** prepare data to be used by MATLAB routines
		* format year so it can be read by MATLAB
		replace ${year} = ${year} - ${year_min_local} + 1 // generate numerical year starting from 1
		rename ${year} year_akm
		label var year_akm "Year (AKM format)"
		
		* generate demographics variable to be interacted with year FEs and age FEs
		if $akm_gender_inter local var_gender = "${gender}"
		if $akm_edu_inter local var_edu = "${edu}"
		if $dem_inter {
			${gtools}egen int dem = group(`var_gender' `var_edu')
			label var dem "Demographic (`var_gender' `var_edu') group"
			tab dem, m
		}
		
		* prepare to outsheet necessary variables
		if $dem_inter global dem_outsheet = "dem"
		else global dem_outsheet = ""
		if $akm_age_poly_order global age_outsheet = "${age}" // if AKM estimation includes higher-order age terms
		else global age_outsheet = ""
		if $akm_hours global hours_outsheet = "${hours}"
		else global hours_outsheet = ""
		if $akm_occ global occ_outsheet = "${occ}"
		else global occ_outsheet = ""
		if $akm_tenure & (!$akm_tenure_top_coding | ${year_min_local} > ${year_min_data} | inlist("${user}", "Chris_office", "Chris_home", "Grid")) global tenure_outsheet = "${tenure}"
		else global tenure_outsheet = ""
		
		* normalize log income in each year to be mean zero by demographic (education x gender) subgroups
		local demean_dem = ""
		if $demean_income_year local demean_dem = "`demean_dem' year_akm"
		if $demean_income_gender local demean_dem = "`demean_dem' ${gender}"
		if $demean_income_edu local demean_dem = "`demean_dem' ${edu}"
		if "`demean_dem'" != "" {
			if "${gtools}" == "" bys `demean_dem': egen float `inc_concept'_mean = mean(`inc_concept')
			else gegen float `inc_concept'_mean = mean(`inc_concept'), by(`demean_dem')
			gen float `inc_concept'_demeaned = `inc_concept' - `inc_concept'_mean
			label var `inc_concept'_demeaned "Demeaned income for AKM decomp. (log)"
			drop `inc_concept'_mean
		}
		
		* format variables for file export
		foreach var of varlist * {
			assert !missing(`var') // make sure that all variables are nonmissing
		}
		local vars_format = "`inc_concept' ${id_pers}* ${id_firm}* year_akm"
		if $akm_age_poly_order local vars_format = "`vars_format' ${age}"
		if $dem_inter local vars_format = "`vars_format' dem"
		if $akm_hours local vars_format = "`vars_format' ${hours}"
		if $akm_occ local vars_format = "`vars_format' ${occ}"
		if $akm_tenure & (!$akm_tenure_top_coding | ${year_min_local} > ${year_min_data})  local vars_format = "`vars_format' ${tenure}"
		foreach var of varlist `vars_format' {
			cap confirm numeric var `var'
			if !_rc {
				sum `var', meanonly
				if "`var'" == "`inc_concept'" format `var' %`=ceil(max(log10(abs(r(min))),log10(abs(r(max))))) + 6'.6f
				else format `var' %`=ceil(max(log10(abs(r(min))),log10(abs(r(max)))))'.0f
			}
		}
		
		* save temp file
		prog_comp_desc_sum_save "${DIR_TEMP}/temp_akm_kss_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta"
		
		
		*** find (regular and leave-one-out) connected sets
		* outsheet list of wages, worker IDs, firm IDs, year IDs, etc. to find connected sets
		local vars_list = "${id_pers} ${id_firm} year_akm"
		keep `vars_list'
		compress
		desc
		sum `vars_list', sep(0)
		sort ${id_pers} year_akm
		export delim `vars_list' using "${DIR_TEMP}/tomatlab_connected.csv", delim(tab) novarnames nolabel replace
		
		* prepare file with parameters for computation of connected set
		clear
		set obs 1
		gen int year_min = ${year_min_local}
		gen int year_max = ${year_max_local}
		if "`inc_concept'" == "${earn}" gen byte inc_concept = 1
		else if "`inc_concept'" == "${wage}" gen byte inc_concept = 2
		gen int thresh = `thresh'
		format year_min year_max %2.0f
		format inc_concept %1.0f
		format thresh %8.0f
		compress
		desc
		sum, sep(0)
		export delim ///
			year_min year_max inc_concept thresh ///
			using "${DIR_TEMP}/parameters_connected.csv", delim(tab) novarnames nolabel replace
		clear
		
		* call MATLAB routine to compute (regular and leave-one-out) connected sets
		cap rm "${DIR_TEMP}/done_connected.txt"
		cap rm "${DIR_TEMP}/success_connected.txt"
		if inlist("`c(os)'", "MacOSX", "Unix") {
			!${FILE_MATLAB} -nosplash -noFigureWindows -nodesktop -nodisplay <"${DIR_CODE}/FUN_CONNECTED.m"
		}
		else if "`c(os)'" == "Windows" {
			!matlab -nosplash -noFigureWindows -batch -wait "run '${DIR_CODE}/FUN_CONNECTED.m'" // OLD: !matlab -nosplash -noFigureWindows -r "run '${DIR_CODE}/FUN_CONNECTED.m'; exit"
		}
		local done_connected = 0
		local first_sleep_cycle = 1
		while !`done_connected' {
			cap confirm file "${DIR_TEMP}/done_connected.txt"
			local done_connected = !_rc
			if !`done_connected' {
				if `first_sleep_cycle' disp "...entering sleep mode until FUN_CONNECTED.m is done."
				sleep 5000
			}
			local first_sleep_cycle = 0
		}
		rm "${DIR_TEMP}/done_connected.txt"
		rm "${DIR_TEMP}/tomatlab_connected.csv"
		cap confirm file "${DIR_TEMP}/success_connected.txt"
		if !_rc rm "${DIR_TEMP}/success_connected.txt"
		else {
			disp as error "USER ERROR: Execution of FUN_CONNECTED.m failed."
			error 1
		}
		
		* convert MATLAB output of (regular and leave-one-out) connected sets into Stata-readable format
		foreach m in "f" "fy" { // loop through models with firm FE ("f") and firm-year FEs ("fy")
			foreach c in "regular" "leave_one_out" { // loop through regular connected set ("regular") and leave-one-out connected set ("leave_one_out")
				import delim using "${DIR_TEMP}/connected_`m'_`c'.txt", clear
				rm "${DIR_TEMP}/connected_`m'_`c'.txt"
				save "${DIR_TEMP}/connected_`m'_`c'_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta", replace
			}
		}
		
		
		*** compute sizes and shares of (regular and leave-one-out) connected sets
		* load temporary data
		use "${DIR_TEMP}/temp_akm_kss_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta", clear
		
		* counts before merging with connected sets
		bys ${id_pers}: gen byte ind_worker = 1 if _n == 1
		label var ind_worker "Ind: Unique worker?"
		bys ${id_firm} year_akm: gen byte ind_firm_year = 1 if _n == 1
		label var ind_firm_year "Ind: Unique firm-year?"
		bys ${id_firm}: gen byte ind_firm = 1 if _n == 1
		label var ind_firm "Ind: Unique firm?"
		disp _newline(1)
		qui count
		local N_wy = r(N)
		local N_wy_disp: di %15.0fc r(N)
		qui count if ind_worker < .
		local N_w = r(N)
		local N_w_disp: di %15.0fc r(N)
		qui count if ind_firm_year < .
		local N_fy = r(N)
		local N_fy_disp: di %15.0fc r(N)
		qui count if ind_firm < .
		local N_f = r(N)
		local N_f_disp: di %15.0fc r(N)
		drop ind_worker ind_firm_year ind_firm
		disp "*** Counts for population before merging with connected sets or dropping singletons:"
		disp "--> Number (share) of worker-years = `N_wy_disp' (100.0%)"
		disp "--> Number (share) of workers      = `N_w_disp' (100.0%)"
		disp "--> Number (share) of firm-years   = `N_fy_disp' (100.0%)"
		disp "--> Number (share) of firms        = `N_f_disp' (100.0%)"
		
		* merge in (regular and leave-one-out) connected sets and mark overlap
		rename ${id_pers} id_pers
		foreach m in "f" "fy" { // loop through models with firm FE ("f") and firm-year FEs ("fy")
			foreach c in "regular" "leave_one_out" { // loop through regular connected set ("regular") and leave-one-out connected set ("leave_one_out")
				qui merge m:1 id_pers using "${DIR_TEMP}/connected_`m'_`c'_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta", keep(master match) gen(connected_`m'_`c')
// 				rm "${DIR_TEMP}/connected_`m'_`c'_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta"
			}
			assert connected_`m'_regular >= connected_`m'_leave_one_out // make sure that leave-one-out connected set is a subset of regular connected set
		}
		assert connected_f_regular >= connected_fy_regular // make sure that regular connected set for firm-year FE model is a subset of regular connected set for firm FE model
		assert connected_f_leave_one_out >= connected_fy_leave_one_out // make sure that leave-one-out connected set for firm-year FE model is a subset of leave-one-out connected set for firm FE model
		rename id_pers ${id_pers}
		
		* counts after merging with various connected sets for various models
		qui save "${DIR_TEMP}/temp_connected.dta", replace
		foreach m in "f" "fy" { // loop through models with firm FE ("f") and firm-year FEs ("fy")
			foreach c in "regular" "leave_one_out" { // loop through regular connected set ("regular") and leave-one-out connected set ("leave_one_out")
				use if connected_`m'_`c' == 3 using "${DIR_TEMP}/temp_connected.dta", clear
				drop connected_`m'_`c'
				qui bys ${id_pers}: gen byte ind_worker = 1 if _n == 1
				label var ind_worker "Ind: Unique worker?"
				qui bys ${id_firm} year_akm: gen byte ind_firm_year = 1 if _n == 1
				label var ind_firm_year "Ind: Unique firm-year?"
				qui bys ${id_firm}: gen byte ind_firm = 1 if _n == 1
				label var ind_firm "Ind: Unique firm?"
				qui count
				local N_wy_conn_`m'_`c' = r(N)
				local N_wy_conn_`m'_`c'_disp: di %15.0fc `=r(N)'
				qui count if ind_worker < .
				local N_w_conn_`m'_`c' = r(N)
				local N_w_conn_`m'_`c'_disp: di %15.0fc `=r(N)'
				qui count if ind_firm_year < .
				local N_fy_conn_`m'_`c' = r(N)
				local N_fy_conn_`m'_`c'_disp: di %15.0fc `=r(N)'
				qui count if ind_firm < .
				local N_f_conn_`m'_`c' = r(N)
				local N_f_conn_`m'_`c'_disp: di %15.0fc `=r(N)'
				drop ind_worker ind_firm_year ind_firm
				local share_wy_conn_`m'_`c': di %5.1f 100*`N_wy_conn_`m'_`c''/`N_wy'
				local share_w_conn_`m'_`c': di %5.1f 100*`N_w_conn_`m'_`c''/`N_w'
				local share_fy_conn_`m'_`c': di %5.1f 100*`N_fy_conn_`m'_`c''/`N_fy'
				local share_f_conn_`m'_`c': di %5.1f 100*`N_f_conn_`m'_`c''/`N_f'
				if "`m'" == "f" local m_str = "firm FEs"
				else if "`m'" == "fy" local m_str = "firm-year FEs"
				if "`c'" == "regular" local c_str = "regular"
				else if "`c'" == "leave_one_out" local c_str = "leave-one-out"
				disp _newline(1)
				disp "*** Counts after merging with `c_str' connected set for `m_str' model:"
				disp "--> Number (share) of worker-years = `N_wy_conn_`m'_`c'_disp' (`share_wy_conn_`m'_`c''%)"
				disp "--> Number (share) of workers      = `N_w_conn_`m'_`c'_disp' (`share_w_conn_`m'_`c''%)"
				disp "--> Number (share) of firm-years   = `N_fy_conn_`m'_`c'_disp' (`share_fy_conn_`m'_`c''%)"
				disp "--> Number (share) of firms        = `N_f_conn_`m'_`c'_disp' (`share_f_conn_`m'_`c''%)"
			}
		}
		drop connected_f_regular connected_f_leave_one_out connected_fy_regular
		rm "${DIR_TEMP}/temp_connected.dta"
		
		
		*** prepare estimation of AKM wage equation and KSS-corrected variance decomposition
		* restrict data to leave-one-out connected set
// 		keep if connected_fy_leave_one_out == 3
// 		drop connected_fy_leave_one_out
		
		* drop singleton worker observations
		local n_before = -1
		local n_after = 0
		while `n_before' != `n_after' {
			qui count
			local n_before = r(N)
			if $drop_singletons_pers qui bys $id_pers: keep if _N > 1
			if $drop_singletons_firm qui bys $id_firm: keep if _N > 1
			if $drop_singletons_firm_year qui bys $id_firm $year: keep if _N > 1
			qui count
			local n_after = r(N)
			if $drop_singletons_pers + $drop_singletons_firm + $drop_singletons_firm_year <= 1 local n_before = `n_after' // if only one criterin is applied, then one loop suffices
		}
		
		* create indicators for unique workers and unique firms
		qui bys ${id_pers}: gen byte ind_worker = 1 if _n == 1
		label var ind_worker "Ind: Unique worker?"
		qui bys ${id_firm} year_akm: gen byte ind_firm_year = 1 if _n == 1
		label var ind_firm_year "Ind: Unique firm-year?"
		qui bys ${id_firm}: gen byte ind_firm = 1 if _n == 1
		label var ind_firm "Ind: Unique firm?"
		qui count
		local N_wy_conn_nosingletons = r(N)
		local N_wy_conn_nosingletons_disp: di %15.0fc `=r(N)'
		qui count if ind_worker < .
		local N_w_conn_nosingletons = r(N)
		local N_w_conn_nosingletons_disp: di %15.0fc `=r(N)'
		qui count if ind_firm_year < .
		local N_fy_conn_nosingletons = r(N)
		local N_fy_conn_nosingletons_disp: di %15.0fc `=r(N)'
		qui count if ind_firm < .
		local N_f_conn_nosingletons = r(N)
		local N_f_conn_nosingletons_disp: di %15.0fc `=r(N)'
		drop ind_worker ind_firm_year ind_firm
		local share_wy_conn_nosingletons: di %5.1f 100*`N_wy_conn_nosingletons'/`N_wy'
		local share_w_conn_nosingletons: di %5.1f 100*`N_w_conn_nosingletons'/`N_w'
		local share_fy_conn_nosingletons: di %5.1f 100*`N_fy_conn_nosingletons'/`N_fy'
		local share_f_conn_nosingletons: di %5.1f 100*`N_f_conn_nosingletons'/`N_f'
		if "`m'" == "f" local m_str = "firm FE"
		else if "`m'" == "fy" local m_str = "firm-year FE"
		if "`c'" == "regular" local c_str = "regular"
		else if "`c'" == "leave_one_out" local c_str = "leave-one-out"
		disp _newline(1)
		disp "*** Counts after merging with leave-one-out connected set for firm-year FEs model and removing singletons:"
		disp "--> Number (share) of worker-years = `N_wy_conn_nosingletons_disp' (`share_wy_conn_nosingletons'%)"
		disp "--> Number (share) of workers      = `N_w_conn_nosingletons_disp' (`share_w_conn_nosingletons'%)"
		disp "--> Number (share) of firm-years   = `N_fy_conn_nosingletons_disp' (`share_fy_conn_nosingletons'%)"
		disp "--> Number (share) of firms        = `N_f_conn_nosingletons_disp' (`share_f_conn_nosingletons'%)"
		
		* outsheet list of wages, worker IDs, firm IDs, year IDs, etc. to estimate wage equation
		local vars_list = "`inc_concept'${demeaned} ${id_pers} ${id_firm} year_akm ${dem_outsheet} ${age_outsheet} ${hours_outsheet} ${occ_outsheet} ${tenure_outsheet}"
		keep `vars_list'
		compress
		desc
		sum `vars_list', sep(0)
		sort ${id_pers} year_akm
		export delim `vars_list' using "${DIR_TEMP}/tomatlab_akm_kss.csv", delim(tab) novarnames nolabel replace
		clear
	}
}
