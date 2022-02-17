# Haplogic validation pipeline (matlab code)

This set of scripts takes the HapLogic validation file output and generates ROC statistics including AUC, PPV, NPV as well as the expected vs. actual plots for all populations as well as all the broad race groups.


Filename | Description
---------|------------
run_all |  Start with this file. Calls the script run_validation.m to run validations for all broad race groups
run_validation.m | This file runs validations on a particular broad race group or "All" groups. Has been updated to only report 9of10 and 10of10 results
roc_AM.m | Calculates ROC statistics
Detailed_expected_vs_Predicted.m | Generates expected vs. actual plots
Weighted_correlation.m | estiamted weighted correltation between expected and actual results
Weighted_distance.m | estiamted weighted distance between expected and actual plots
partest.m | Used for the ROC calculations

## To run validation at the Matlab prompt type:
```run_all('test_validation.txt', 0.5, 0.01)```

## Output
For each broad race, the output is a .csv file with all the ROC statistics, weighted correlation and distance and std. error as well as an image plot of expected vs. actual 