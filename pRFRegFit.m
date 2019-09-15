function pRFRegFit(outputStruct, x, y, z)
%UNTITLED2 display model fits for single voxel
%   Detailed explanation goes here
figure

% The residuals are simply the difference between model and data:
p2_ = plot(outputStruct.t, outputStruct.data{x, y, z}.dataZ, 'k-', ...
           outputStruct.t, outputStruct.data{x, y, z}.lassoZ, 'm', ...
           outputStruct.t, outputStruct.data{x, y, z}.ridgeZ, 'c', ...
           outputStruct.t, outputStruct.data{x, y, z}.svmZ, 'b');

set(p2_, 'LineWidth',2);
legend(p2_, {'data', 'lasso', 'ridge', 'svm'})
end

