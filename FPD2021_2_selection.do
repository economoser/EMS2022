********************************************************************************
* DESCRIPTION: Loads (test, SWE, or BRA) data and makes preliminary sample
*              selection.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 27, 2021.
********************************************************************************


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		
		* load data
		local vars_list = "${id_pers} ${year} ${id_firm} ${gender} ${edu} ${age} `inc_concept'" // Note: other variables in SWE data include ${wage}, ${occ}, ${hours}, ${degree_year}.
		if $akm_hours local vars_list = "`vars_list' ${hours}"
		if $akm_occ local vars_list = "`vars_list' ${occ}"
		if $akm_tenure local vars_list = "`vars_list' ${tenure}"
		local vars_list_nogender = subinstr("`vars_list'", "${gender}", "", .) // Note: SWE data does not contain gender but has two separate files -- one for each gender.
		if $test use `vars_list' ${gender} using "${DIR_TEMP}/test_data.dta", clear
		else {
			if "${user}" == "SWE" {
				use `vars_list_nogender' using "${RawDataDir}/CombinedMainSpellWomen.dta", clear
				gen byte ${gender} = 2
				save "${DIR_TEMP}/temp_selection.dta", replace
				use `vars_list_nogender' using "${RawDataDir}/CombinedMainSpellMen.dta", clear
				gen byte ${gender} = 1
				append using "${DIR_TEMP}/temp_selection.dta"
				rm "${DIR_TEMP}/temp_selection.dta"
				label var ${gender} "Gender"
				label define gen_l 1 "Male" 2 "Female", replace
				label val ${gender} gen_l
				order `vars_list'
			}
			else if inlist("${user}", "Chris_office", "Chris_home", "Grid") {
				clear
				local vars_append = "persid year empid_firm gender edu age earn_mean_mw"
				if $akm_hours local vars_append = "`vars_append' hours"
				if $akm_occ local vars_append = "`vars_append' occ94_5"
				if $akm_tenure local vars_append = "`vars_append' tenure"
				forval y = $year_min_data/$year_max_data {
					if $sample_read append using "${RawDataDir}/`y'/sample_clean`y'.dta", keep(`vars_append')
					else append using "${RawDataDir}/`y'/clean`y'.dta", keep(`vars_append')
				}
				rename persid ${id_pers}
				rename year ${year}
				rename empid_firm ${id_firm}
				rename gender ${gender}
				rename edu ${edu}
				rename age ${age}
				rename earn_mean_mw `inc_concept'
				if $akm_hours rename hours ${hours}
				if $akm_occ rename occ94_5 ${occ}
				if $akm_tenure rename tenure ${tenure}
			}
		}
		
		* early selection criteria
		keep if ///
			inrange(${year}, ${year_min_data}, ${year_max_data}) ///
			& inrange(${gender}, ${gender_min}, ${gender_max}) ///
			& inrange(${age}, ${age_min}, ${age_max}) ///
			
		* make sure that key variables have nonmissing values
		foreach var in ${year} ${id_pers} ${id_firm} ${gender} ${age} ${edu} {
			cap confirm var `var', exact
			if !_rc keep if !missing(`var')
		}
		
		* make sure that some data was loaded into memory
		count
		if r(N) == 0 {
			disp as error "USER ERROR: No data was loaded into memory."
			error 1
		}
		foreach var in "${id_pers}" "${year}" "${id_firm}" "${gender}" "${edu}" "${age}" "`inc_concept'" {
			cap confirm var `var', exact
			if _rc {
				disp as error "USER ERROR: Variable `var' does not exist in data."
				error 1
			}
		}
		
		* recode time variables (year, quarter, month)
		local n_period_vars = 0
		foreach var in "${year}" "${quarter}" "${month}" {
			cap confirm var `var', exact
			if !_rc {
				if `n_period_vars' == 0 local period_vars = "`var'"
				else local period_vars = "`period_vars' `var'"
				local ++n_period_vars
			}
		}
		assert "`period_vars'" != ""
		
		* recode worker ID
		if $test | "${user}" == "SWE" {
			rename ${id_pers} ${id_pers}_original
			label var ${id_pers}_original "Worker ID (string, original)"
			${gtools}egen long ${id_pers} = group(${id_pers}_original)
			label var ${id_pers} "Worker ID (numeric)"
			keep if !missing(${id_pers}_original)
			local vars_list = subinstr("`vars_list'", "${id_pers}", "${id_pers} ${id_pers}_original", .)
		}
		
		* recode employer ID
		if $test | "${user}" == "SWE" {
			rename ${id_firm} ${id_firm}_original
			label var ${id_firm}_original "Firm ID (string, original)"
			${gtools}egen double ${id_firm} = group(${id_firm}_original)
			label var ${id_firm} "Firm ID (numeric)"
			keep if !missing(${id_firm}_original)
			cap confirm str var ${id_firm}_original
			if !_rc keep if ${id_firm}_original != "00000000"
			local vars_list = subinstr("`vars_list'", "${id_firm}", "${id_firm} ${id_firm}_original", .)
		}
		
		* recode education
		if "${user}" == "SWE" {
			recode ${edu} (1=2) (7=6) // pool primary school (1) with middle school (2), and pool Masters/PhDs (7) with Bachelors (6)
			label define edu_l 2 "Less than High School" 3 "High School" 4 "Some College" 5 "College" 6 "Postgraduate", replace
			label val ${edu} edu_l
			if "${gtools}" == "" bys ${id_pers}: egen byte edu_perm = max(${edu})
			else gegen byte edu_perm = max(${edu}), by(${id_pers})
			replace ${edu} = edu_perm
			drop edu_perm
		}
// 		else if inlist("${user}", "Chris_office", "Chris_home", "Grid") recode edu (1 = 2) // NOTE FROM 10/27/2021: May need this in order to avoid -ichol()- error during AKM estimation in MATLAB.
		
		* recode worker race
		
		* recode worker age
		
		* recode region
		
		* recode industry
		
		* recode occupation
		if $akm_occ & inlist("${user}", "Chris_office", "Chris_home", "Grid") {
			replace ${occ} = 99999 if ${occ} == .
			if $akm_occ_coarsen replace ${occ} = floor(${occ}/10) // NOTE FROM 10/05/2021: May need this in order for KSS corrections to avoid -ichol()- error during AKM estimation in MATLAB.
		}
		
		* recode hours
		if $akm_hours & inlist("${user}", "Chris_office", "Chris_home", "Grid") {
			replace ${hours} = 40 if hours == .
		}
		
		* recode income measures (earnings, wage, etc.)
		if "`inc_concept'" == "${earn}" & inlist("${user}", "Chris_office", "Chris_home", "Grid") replace ${earn} = 120 if ${earn} > 120 & ${earn} < .
		
		* compute tenure
		if $akm_tenure & "${user}" == "SWE" {
			bys ${id_pers} ${id_firm} ${year}: gen byte count_tenure = 1 if _n == 1
			bys ${id_pers} ${id_firm} (${year} count_tenure): gen byte ${tenure} = sum(count_tenure)
			label var ${tenure} "Tenure (years)"
			drop count_tenure
			if $akm_tenure_top_coding & ${year_min_local} > ${year_min_data} replace ${tenure} = min(${tenure}, ${year_min_local} - ${year_min_data} + 1) if ${tenure} < .
		}
		else if $akm_tenure & inlist("${user}", "Chris_office", "Chris_home", "Grid") {
			replace ${tenure} = ceil(${tenure}/12)
		}
		
		* deflate income measures (earnings, wage, etc.)
		
		* convert income measures (earnings, wage, etc.) from levels to natural logarithms
		rename `inc_concept' `inc_concept'_level
		gen float `inc_concept' = ln(`inc_concept'_level)
		if "${user}" == "SWE" {
			if "`inc_concept'" == "${earn}" label var ${earn} "Monthly earnings (log SEK)"
			else if "`inc_concept'" == "${wage}" label var ${wage} "Hourly wage (log SEK)"
		}
		else {
			if "`inc_concept'" == "${earn}" label var ${earn} "Earnings (log)"
			else if "`inc_concept'" == "${wage}" label var ${wage} "Wage (log)"
		}
		drop `inc_concept'_level
		
		* impose minimum and maximum worker income thresholds
		if "${lower_perc}" != "" & "${user}" == "SWE" {
			if inlist(${lower_perc}, 1, 5, 10, 25, 50, 75, 90, 95, 99) {
				forval y = $year_min_data/$year_max_data {
					qui sum `inc_concept' if ${year} == `y', d
					drop if `inc_concept' < r(p${lower_perc}) & ${year} == `y'
				}
			}
			else if "${gtools}" != "" {
				local n_quantiles = ceil(100/${lower_perc})
				gquantiles q = `inc_concept', xtile n(`n_quantiles') by(year)
				drop if q == 1
				drop q
			}
			else {
				disp as error "USER ERROR: Cannot recognize lower earnings percentile (lower_perc = P${lower_perc})."
				error 1
			}
		}
		if $drop_below_mw {
			if inlist("${user}", "Chris_office", "Chris_home", "Grid") & "`inc_concept'" == "${earn}" {
				drop if ${earn} < 0
			}
		}
		if ${`inc_concept'_min} > 0 keep if exp(`inc_concept') >= ${`inc_concept'_min}
		if ${`inc_concept'_max} < . keep if exp(`inc_concept') <= ${`inc_concept'_max}
		
		
		*** impose additional selection criteria
		* make sure that additional variables have nonmissing values
		foreach var in ${hours} ${occ} ${tenure} `inc_concept' {
			cap confirm var `var', exact
			if !_rc keep if !missing(`var')
		}
		
		* drop worker-firm matches that are potentially not observed for full time period
		
		* keep only observations with highest earnings for each person-year
		if $test | inlist("${user}", "Chris_office", "Chris_home", "Grid") { // Note: SWE data already has one observation per person-year preselected
			gen double rand = runiform()
			bys `period_vars' ${id_pers} (`inc_concept' rand): keep if _n == _N
			drop rand
		}
		
		* drop observations with excessive earnings growth
		if ${inc_growth_min} < . bys ${id_pers} (`period_vars'): keep if `inc_concept'[_n] - `inc_concept'[_n - 1] > ${inc_growth_min}
		if ${inc_growth_max} < . bys ${id_pers} (`period_vars'): keep if `inc_concept'[_n] - `inc_concept'[_n - 1] < ${inc_growth_max} | `inc_concept'[_n - 1] == .
		
		* create unique observation ID -- Note: at this point all data from all periods and all states must be loaded!
// 		gen double rand = runiform()
// 		sort `period_vars' ${id_firm} ${id_pers} `inc_concept' rand
// 		drop rand
// 		gen double ${id_unique} = _n
// 		label var ${id_unique} "Unique observation ID"
		
		* impose minimum firm size threshold
		if !${test} & `thresh' > 1 bys `period_vars' ${id_firm}: drop if _N < `thresh' // Note: Do not impose minimum firm size threshold on test data!

		* keep only firms with enough workers
		if $min_size_emp > 1 {
			bys ${id_firm} ${year}: gen long fsize = _N
			label var fsize "Firm size (number of employees)"
			keep if fsize >= $min_size_emp
			drop fsize
		}

		* keep only firms with enough switchers -- Note: Technically, this should be imposed iteratively after finding the connected set and before re-computing the connected set, etc., but this is approximately equivalent.
		if $min_switchers > 1 {
			bys ${id_pers} (${year}): gen byte ind_switcher = (${id_firm}[_n] != ${id_firm}[_n + 1] & !missing(${id_firm}[_n + 1])) | (${id_firm}[_n] != ${id_firm}[_n - 1] & !missing(${id_firm}[_n - 1]))
			label var ind_switcher "Ind: Switched firms between current and previous job spells?"
			bys ${id_firm}: egen long n_switcher = total(ind_switcher)
			label var n_switcher "Number of switchers at a given firm"
			keep if n_switcher >= $min_switchers
			drop n_switcher ind_switcher
		}
		
		
		*** save cleaned sample
		* save
		order `vars_list'
		prog_comp_desc_sum_save "${DIR_TEMP}/selection_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta"
	}
}
