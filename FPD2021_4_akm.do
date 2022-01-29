********************************************************************************
* DESCRIPTION: Calls MATALB file FUN_AKM.m to estimate and store wage components
*              from AKM (firm FE and firm-year FE) models, then calls MATLAB
*              file FUN_AKM_KSS.m to estimate KSS bias corrections for AKM (firm
*              FE and firm-year FE) models.
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


*** loop through different individual income concepts, minimum firm size thresholds, and models ("f" = firm FEs, "fy" = firm-year FEs)
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		foreach m in "fy" "f" {
			disp _newline(5)
			if "`m'" == "f" local m_str = "firm FE"
			else if "`m'" == "fy" local m_str = "firm-year FE"
			disp "--> inc_concept = `inc_concept', thresh = `thresh', model = `m_str'"
			
			* prepare file with parameters for estimation of AKM wage equation
			clear
			set obs 1
			gen int year_min = ${year_min_local}
			gen int year_max = ${year_max_local}
			gen byte akm_age_poly_order = ${akm_age_poly_order}
			gen int akm_age_flat_min = ${akm_age_flat_min}
			gen int akm_age_flat_max = ${akm_age_flat_max}
			gen int akm_age_norm = ${akm_age_norm}
			gen int age_min = ${age_min}
			gen int age_max = ${age_max}
			gen byte akm_year_dummies = ${akm_year_dummies}
			gen byte dem_inter = ${dem_inter}
			gen byte akm_hours = ${akm_hours}
			gen byte akm_occ = ${akm_occ}
			gen byte akm_tenure = (${akm_tenure} & (!$akm_tenure_top_coding | ${year_min_local} > ${year_min_data}))
			if "`m'" == "f" gen byte akm_model = 1
			else if "`m'" == "fy" gen byte akm_model = 2
			gen byte save_space = ${save_space}
			gen int parallel_max = ${parallel_max}
			format akm_age_poly_order dem_inter akm_hours akm_occ akm_tenure akm_model save_space %1.0f
			format akm_age_flat_min akm_age_flat_max akm_age_norm age_min age_max parallel_max %3.0f
			format year_min year_max %4.0f
			compress
			export delim ///
				year_min year_max akm_age_poly_order akm_age_flat_min akm_age_flat_max akm_age_norm age_min age_max akm_year_dummies dem_inter akm_hours akm_occ akm_tenure akm_model save_space parallel_max ///
				using "${DIR_TEMP}/parameters_akm_kss.csv", delim(tab) novarnames nolabel replace
			clear
			
			* call MATLAB via shell to estimate AKM wage equation (without KSS correction)
			cap rm "${DIR_TEMP}/done_akm.txt"
			cap rm "${DIR_TEMP}/success_akm.txt"
			if inlist("`c(os)'", "MacOSX", "Unix") {
				!${FILE_MATLAB} -nosplash -noFigureWindows -nodesktop -nodisplay <"${DIR_CODE}/FUN_AKM.m"
			}
			else if "`c(os)'" == "Windows" {
				!matlab -nosplash -noFigureWindows -batch -wait "run '${DIR_CODE}/FUN_AKM.m'" // OLD: !matlab -nosplash -noFigureWindows -r "run '${DIR_CODE}/FUN_AKM.m'; exit"
			}
			local done_akm = 0
			local first_sleep_cycle = 1
			while !`done_akm' {
				cap confirm file "${DIR_TEMP}/done_akm.txt"
				local done_akm = !_rc
				if !`done_akm' {
					if `first_sleep_cycle' disp "...entering sleep mode until FUN_AKM.m is done."
					sleep 5000
				}
				local first_sleep_cycle = 0
			}
			rm "${DIR_TEMP}/done_akm.txt"
			cap confirm file "${DIR_TEMP}/success_akm.txt"
			if !_rc rm "${DIR_TEMP}/success_akm.txt"
			else {
				disp as error "USER ERROR: Execution of FUN_AKM.m failed."
				error 1
			}
			
			* convert MATLAB output of components of AKM wage equation into Stata-readable format
			import delim using "${DIR_TEMP}/tostata_akm.txt", clear
			rm "${DIR_TEMP}/tostata_akm.txt"
			rename id_pers ${id_pers}
			label var ${id_pers} "Worker ID (numeric)"
			rename y year_akm
			label var year_akm "Year"
			label var pe "Predicted AKM worker FE"
			rename fe ${`m'_fe}
			label var ${`m'_fe} "Predicted AKM `m_str'"
			if $akm_year_dummies {
				rename xb_y xb_year
				if !$akm_edu_inter & !$akm_gender_inter label var xb_year "Predicted AKM year FE"
				else if ($akm_edu_inter | $akm_gender_inter) label var xb_year "Predicted AKM demographics-year FE"
			}
			if $akm_age_poly_order {
				rename xb_a xb_age
				if $akm_age_poly_order == 1 {
					if !$akm_edu_inter & !$akm_gender_inter label var xb_age "Predicted AKM age FE"
					else label var xb_age "Predicted AKM demographics-age FE"
				}
				else if $akm_age_poly_order >= 2 {
					if !$akm_edu_inter & !$akm_gender_inter label var xb_age "Predicted AKM higher-order age terms"
					else label var xb_age "Predicted AKM higher-order demographics-age terms"
				}
			}
			if $akm_hours {
				rename xb_h xb_hours
				label var xb_hours "Predicted AKM hours FE"
			}
			if $akm_occ {
				rename xb_o xb_occ
				label var xb_occ "Predicted AKM occupation FE"
			}
			if ${akm_tenure} & (!$akm_tenure_top_coding | ${year_min_local} > ${year_min_data}) {
				rename xb_t xb_tenure
				label var xb_tenure "Predicted AKM tenure FE"
			}
			
			* recast AKM estimates to float data types
			foreach var in pe ${`m'_fe} xb_year xb_age xb_hours xb_occ xb_tenure {
				cap confirm var `var', exact
				if !_rc recast float `var', force
			}
			
			* add other variables from temp file
			merge 1:1 ${id_pers} year_akm using "${DIR_TEMP}/temp_akm_kss_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta", keep(match master) nogen
			
			* clean up demographics variable and un-transform year variable
			if (${akm_gender_inter} | ${akm_edu_inter}) drop dem
			replace year_akm = year_akm + ${year_min_local} - 1
			rename year_akm ${year}
			
			* generate residual
			gen float resid = `inc_concept'${demeaned}
			foreach var in pe ${`m'_fe} xb_year xb_age xb_hours xb_occ xb_tenure {
				cap confirm var `var', exact
				if !_rc replace resid = resid - `var'
			}
			label var resid "Predicted AKM residual"
			
			* compute variance-covariance matrix and variance decomposition
			local cov_list = ""
			foreach var in `inc_concept'${demeaned} pe ${`m'_fe} xb_year xb_age xb_hours xb_occ xb_tenure resid {
				cap confirm var `var', exact
				if !_rc local cov_list = "`cov_list' `var'"
			}
			disp _newline(1)
			disp "corr pe ${`m'_fe}:"
			corr pe ${`m'_fe}
			local rho_pe_${`m'_fe}: di %4.3f r(rho)
			disp _newline(1)
			disp "corr `cov_list':"
			corr `cov_list'
			disp _newline(1)
			disp "corr pe ${`m'_fe}, cov:"
			corr pe ${`m'_fe}, cov			
			local cov_pe_${`m'_fe}: di %4.3f r(cov_12)
			disp _newline(1)
			disp "corr `cov_list', cov:"
			corr `cov_list', cov
			matrix C = r(C)
			local cov_counter = 1
			foreach var of local cov_list {
				local var_`var' = C[`cov_counter',`cov_counter']
				local var_`var': di %4.3f `var_`var''
				local ++cov_counter
				if "`var'" == "`inc_concept'${demeaned}" local cov = `var_`inc_concept'${demeaned}'
				else local cov = `cov' - `var_`var''
			}
			local cov: di %4.3f `cov'
			qui sum `inc_concept'
			local var_`inc_concept': di %4.3f r(Var)
			foreach var of local cov_list {
				local var_share_`var' = 100*`var_`var''/`var_`inc_concept'${demeaned}'
				local var_share_`var': di %4.1f `var_share_`var''
			}
			local var_share_cov = 100*`cov'/`var_`inc_concept'${demeaned}'
			local var_share_cov: di %4.1f `var_share_cov'
			local var_share_`inc_concept' = 100*`var_`inc_concept''/`var_`inc_concept'${demeaned}'
			local var_share_`inc_concept': di %4.1f `var_share_`inc_concept''
			local first = 1
			if "${demeaned}" == "" local `inc_concept'_not_demeaned = ""
			else local `inc_concept'_not_demeaned = "`inc_concept'"
			foreach var in ``inc_concept'_not_demeaned' `cov_list' cov resid2 {
				if `first' {
					disp _newline(1)
					disp "--> variance decomposition:"
				}
				if !inlist("`var'", "cov", "resid", "resid2") disp "Var(`var') = `var_`var'' (`var_share_`var''%)"
				else if "`var'" == "cov" disp "2*sum(Cov(.)) = `cov' (`var_share_cov'%)"
				else if "`var'" == "resid2" disp "Var(resid) = `var_resid' (`var_share_resid'%)"
				local first = 0
			}
			disp "Corr(worker FE, `m_str') = `rho_pe_${`m'_fe}'"
			disp "Cov(worker FE, `m_str') = `cov_pe_${`m'_fe}'"
			
			* save complete file, incl. AKM estimates
			foreach var in xb_year xb_age {
				cap confirm var `var', exact
				if !_rc local var_`var' = "`var'"
				else local var_`var' = ""
			}
			order ${id_pers}* ${year} ${id_firm}* ${gender} ${edu} ${age} `inc_concept' `inc_concept'${demeaned} pe ${`m'_fe} `var_xb_year' `var_xb_age' resid
			sort ${id_pers} ${year}
			prog_comp_desc_sum_save "${DIR_TEMP}/akm_estimates_`m'_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta"
			
			* call MATLAB via shell to estimate KSS-corrected variance decomposition based on AKM wage equation
			clear
			cap rm "${DIR_TEMP}/done_akm_kss.txt"
			cap rm "${DIR_TEMP}/success_akm_kss.txt"
// 			local dir_list: dir "${DIR_TEMP}/" dirs "ext_*"
// 			foreach d of local dir_list {
// 				rm "${DIR_TEMP}/`d'"
// 			}
			if inlist("`c(os)'", "MacOSX", "Unix") {
				!${FILE_MATLAB} -nosplash -noFigureWindows -nodesktop -nodisplay <"${DIR_CODE}/FUN_AKM_KSS.m"
			}
			else if "`c(os)'" == "Windows" {
				!matlab -nosplash -noFigureWindows -batch -wait "run '${DIR_CODE}/FUN_AKM_KSS.m'" // OLD: !matlab -nosplash -noFigureWindows -r "run '${DIR_CODE}/FUN_AKM_KSS.m'; exit"
			}
			local done_akm_kss = 0
			local first_sleep_cycle = 1
			while !`done_akm_kss' {
				cap confirm file "${DIR_TEMP}/done_akm_kss.txt"
				local done_akm_kss = !_rc
				if !`done_akm_kss' {
					if `first_sleep_cycle' disp "...entering sleep mode until FUN_AKM_KSS.m is done."
					sleep 5000
				}
				local first_sleep_cycle = 0
			}
			rm "${DIR_TEMP}/done_akm_kss.txt"
			cap confirm file "${DIR_TEMP}/success_akm_kss.txt"
			if !_rc rm "${DIR_TEMP}/success_akm_kss.txt"
			else {
				disp as error "USER ERROR: Execution of FUN_AKM_KSS.m failed."
				error 1
			}
		}
		rm "${DIR_TEMP}/temp_akm_kss_`inc_concept'_`thresh'_${year_min_local}_${year_max_local}.dta"
	}
}
