# README for [Engbom, Moser, and Sauermann (2022)](https://ssrn.com/abstract=3531250)


## Guide to Files

The RAMS, LISA, LOUISE, and FEK datasets are confidential and therefore not uploaded as part of the files. The code to replicate all results has been included in the following files:

1. **[FPD2021_0_master.do](FPD2021_0_master.do)**: This is the master file. All other files should be executed from this file with the appropriate switches.
2. **[FPD2021_1_test_data.do](FPD2021_1_test_data.do)**: Simulate test data that can be used to test-run all other code.
3. **[FPD2021_2_selection.do](FPD2021_2_selection.do)**: Impose sample selection criteria.
4. **[FPD2021_3_connected.do](FPD2021_3_connected.do)**: Form connected set by calling [FUN_CONNECTED.m](FUN_CONNECTED.m).
5. **[FPD2021_4_akm.do](FPD2021_4_akm.do)**: Estimate AKM wage equation and store estimates by calling [FUN_AKM.m](FUN_AKM.m), then conduct leave-one-out bias correction due to [Kline, Saggio, and Solvsten (2020)](https://www.econometricsociety.org/publications/econometrica/2020/09/01/leave-out-estimation-variance-components) by calling [FUN_AKM_KSS.m](FUN_AKM_KSS.m).
6. **[FPD2021_5_panels.do](FPD2021_5_panels.do)**: Create panel datasets containing AKM estimates, firm financials data, etc.
7. **[FPD2021_6_tables.do](FPD2021_6_tables.do)**: Produces output tables to be included in the paper.
8. **[FPD2021_7_comparison.do](FPD2021_7_comparison.do)**: Compares firm-level mean earnings, firm FEs, and firm-year FEs.
9. **[FPD2021_8_time_series.do](FPD2021_8_time_series.do)**: Analyzes time series of the firm-year pay distribution.
10. **[FPD2021_9_lifecycle.do](FPD2021_9_lifecycle.do)**: Investigates life cycle profiles of firm pay.
11. **[FPD2021_X_additional_figures.do](FPD2021_X_additional_figures.do)**: Produces additional figures.
12. **[FUN_AKM_KSS.m](FUN_AKM_KSS.m)**: Conduct leave-one-out bias correction due to [Kline, Saggio, and Solvsten (2020)](https://www.econometricsociety.org/publications/econometrica/2020/09/01/leave-out-estimation-variance-components) by calling routines stored in [_LeaveOutTwoWay_3_02_Chris](_LeaveOutTwoWay_3_02_Chris).
13. **[FUN_AKM.m](FUN_AKM.m)**: Estimate AKM wage equation.
14. **[FUN_CONNECTED.m](FUN_CONNECTED.m)**: Form connected set.
15. **[TablesForPaper.m](TablesForPaper.m)**: Produce tables for paper.
16. **[comma_separator.m](comma_separator.m)**: Auxiliary file to produce comma-separated output files, based on code by user Image Analyst posted in December 2016 on the [MATLAB Answers Forum](https://www.mathworks.com/matlabcentral/answers/315519-how-to-display-data-with-commas-or-spaces).
17. **[csv2mat_block.m](csv2mat_block.m)**: Auxiliary file to produce a block matrix, based on code by S. Hsiang (smh2137@columbia.edu) written in May 2010.
18. **[csv2mat_numeric.m](csv2mat_numeric.m)**: Auxiliary file to produce a structure, based on code by S. Hsiang (smh2137@columbia.edu) written in May 2010.
19. **[_LeaveOutTwoWay_3_02_Chris](_LeaveOutTwoWay_3_02_Chris)**: Lightly edited version of the Github repository by [Raffaele Saggio](https://sites.google.com/site/raffaelesaggio/), available as [repository LeaveOutTwoWay by rsaggio87 on Github](https://github.com/rsaggio87/LeaveOutTwoWay).


## References

[Engbom, Niklas, Christian Moser, and Jan Sauermann. 2022. "Firm Pay Dynamics." Journal of Econometrics, forthcoming.](https://ssrn.com/abstract=3531250)
