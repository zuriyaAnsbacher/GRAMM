% example: run_all('test_validation.txt', 0.5, 0.01)

function run_all(filename, roc_thresh, roc_alpha)

run_validation(filename, 'ALL' ,roc_thresh, roc_alpha)
run_validation(filename, 'CAU' ,roc_thresh, roc_alpha)
run_validation(filename, 'HIS' ,roc_thresh, roc_alpha)
run_validation(filename, 'AFA' ,roc_thresh, roc_alpha)
run_validation(filename, 'API' ,roc_thresh, roc_alpha)
run_validation(filename, 'NAM' ,roc_thresh, roc_alpha)
