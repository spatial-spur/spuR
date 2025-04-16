% Save this as lvech_test.m in the same directory
% Test script to verify lvech implementation between MATLAB and R

% Create the same test matrix as in R
test_matrix = reshape(1:16, 4, 4);
disp('Test matrix:');
disp(test_matrix);

% Apply lvech function
matlab_result = lvech(test_matrix);
disp('MATLAB lvech result:');
disp(matlab_result);