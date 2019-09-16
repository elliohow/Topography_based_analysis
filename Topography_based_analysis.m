%% Lee method for pRF analysis
% Dependencies:
% mrTools
% mgl
% MATLAB version R2018a

% GeneralRequirements:
% Motion compensated tSeries

% Requirements for ROI analysis:
% roiCoords which contains the coordinates for the ROI

%% Script setup

% Settings to speed up script
forceStimulusSetup = false; %Select true to always import stimulus file
forceTSeriesSetup = false; % Select true to always import time series

% HRF settings
hrfDisplay = false; % Select true to view HRF

% Stimuli settings
xFlip = true; % Select true to xFlip stimulus
yFlip = false; % Select true to yFlip stimulus
stimulusView = false; % Select true to view montage of stimuli

% Analysis settings
barStimOnly = true; % Select false to compare bar only results with all stimuli
voxelSelection = "roi"; % Select all, roi or single
singleVoxel = [34, 10, 6]; % Voxel to use if voxelSelection set to "single"
rng(42) % Random seed

% Output settings
pRF_2D_disp = true; % Select true to display pRF image for single voxel
pRF_fit_disp = true; % Select true to display pRF fit for single voxel
saveOutput = false; % Select true to save output as .mat file
clearWorkspace = false; % Select true to clear workspace after running script

%% Loop setup

% If barStimOnly is set to true, loop the script once using only bar
% stimuli and time series data, otherwise loop twice with all data.
if barStimOnly == true
    scriptCounter = 1;
else
    scriptCounter = 2;
end

%% Beginning of script loop

tic

while scriptCounter ~= 0
    
    if saveOutput == true
        fileName = input('What do you want to name the saved .mat file? ', 's');
    end
    
    %% Stimulus retrieval
    if forceStimulusSetup == true || barStimOnly == false || ~exist('stimuli', 'var')
        
        if scriptCounter == 2
            cd('Etc');
        else
            cd('BarStim')
        end
        
        tempStimFileList = dir('*.mat');
        tempStimFileNames = {tempStimFileList.name};
        indivStimuli = pRFGetStimImageFromStimfile(tempStimFileNames);
        
        cd ..
        %% HRF Creation
        
        TR = mean(diff(indivStimuli{1}.t));
        t_hrf = 0:TR:30;
        
        % hrf = fmribHRF(t_hrf); % Use this line instead to have an undershoot
        hrf = getGammaHRF(t_hrf); % hard-coded paramters from mrTools pRFFit.m
        hrf = hrf ./ sum(hrf(:)); % normalisation for convolution to make area under curve = 1
        
        if hrfDisplay == true
            plot(t_hrf, hrf, 'b', 'linewidth',2)
            xlabel('Time (s)');
            ylabel('HRF impulse response')
        end
        
        % other idea: once the code is running / stable.. think about using PARFOR
        
        %% Stimuli struct setup
        stimNumber = size(indivStimuli, 2);
        
        tempStimuli = squeeze(struct2cell(cat(1, indivStimuli{:})));
        
        stimuli.x =  tempStimuli{2, 1};
        stimuli.y =  tempStimuli{3, 1};
        stimuli.im = cat(3, tempStimuli{4, :});
        
        stimTimePoints = cat(2, tempStimuli{1, :});
        stimTimePoints = size(stimTimePoints, 2);
        stimuli.t = 1:TR:stimTimePoints*TR;
        
        % Vectors
        x = stimuli.x(:,1,1);
        y = squeeze(stimuli.y(1,:,1));
        
        % Mesh versions
        mx = stimuli.x;
        my = stimuli.y;
        
        %% Flip stimuli
        if xFlip == true
            stimuli.im = flipud(stimuli.im);
            for iC = 1:stimNumber
                indivStimuli{iC}.im = flip(indivStimuli{iC}.im, 1);
            end
            disp('(Lee_method) ==xFlipped stimuli.==')
        end
        
        if yFlip == true
            stimuli.im = fliplr(stimuli.im);
            for iC = 1:stimNumber
                indivStimuli{iC}.im = flip(indivStimuli{iC}.im, 2);
            end
            disp('(Lee_method) ==yFlipped stimuli.==')
        end
        
    end
    
    %% Stimulus view
    
    if stimulusView == true
        figure
        montage(permute(stimuli.im, [2,1,3]), ...
            'BorderSize', [2,2], ...
            'BackgroundColor',[1,1,1]*0.5)
    end
    
    %% K matrix setup
    
    % Response = K matrix * weights
    %
    % R = K * b
    %
    % R is measured, K is defined above, desired: best **b**
    K = cell(1, stimNumber);
    H = cell(1, stimNumber);
    
    % Calculate the K matrix for each stimulus block separately
    for iC = 1:stimNumber
        [K{iC}, H{iC}] = estimate_prf_linear_transform(indivStimuli{iC}, hrf);
    end
    
    %% Collect time series
    
    if forceTSeriesSetup == true || barStimOnly == false || ~exist('tSeries', 'var')
        
        if scriptCounter == 2
            cd('MotionComp/TSeries')
        else
            cd('BarTimeSeries/TSeries')
        end
        
        tempTSeriesList = dir('*.img');
        tempTSeriesNames = {tempTSeriesList.name};
        
        tSeries = cell(1, stimNumber);
        
        for iC = 1:stimNumber
            tSeries{iC} = mlrImageReadNifti(tempTSeriesNames{iC});
        end
        
        disp('(Lee_method) tSeries imported.')
        
        cd ../..
        
    end
    
    %% Single voxel selection
    
    tSeriesSize = [size(tSeries{1}, 1), size(tSeries{1}, 2), size(tSeries{1}, 3)];
    
    if voxelSelection == "single"
        tSeriesVoxel = cell(1, stimNumber);
        
        for iC = 1:stimNumber
            tSeriesVoxel{iC} = squeeze(tSeries{iC}(singleVoxel(1), ...
                singleVoxel(2), singleVoxel(3), :));
        end
        
        tSeriesSize = [1, 1, 1];
    end
    
    %% Output struct setup
    
    outputStruct = struct('Voxel_selection', voxelSelection,  ...
        'hrf', hrf, 'x', stimuli.x, 'y', stimuli.y, 't', stimuli.t, ...
        'image', stimuli.im, 'K', {K});
    
    if scriptCounter == 2
        outputStruct.stimuli = 'All_stimuli';
    else
        outputStruct.stimuli = 'Bar_stimuli';
    end
    
    if xFlip == true
        outputStruct.xFlip = true;
    end
    
    if yFlip == true
        outputStruct.yFlip = true;
    end
    
    %% Analysis
    
    % Pack params
    params = {tSeries, K, stimuli, hrf, stimNumber, ...
        voxelSelection, mx, outputStruct};
    
    if voxelSelection == "roi"
        params{8}.data = cell(1, size(ROICoords, 2));
        outputStruct = leeAnalysis(params, ROICoords);
        
    else
        params{8}.data = cell(tSeriesSize(1), tSeriesSize(2), tSeriesSize(3));
        
        % Run the analysis
        if voxelSelection == "all"
            outputStruct = leeAnalysis(params);
            
        elseif voxelSelection == "single"
            outputStruct = leeAnalysis(params, singleVoxel);
        end
    end
    
    %% Display fits
    
    if pRF_fit_disp == true && voxelSelection == "single"
        pRFRegFit(outputStruct, 1, 1, 1)
    end
    
    if pRF_2D_disp == true && voxelSelection == "single"
        pRFDisplay(outputStruct, 1, 1, 1)
    end
    
    %% Save output
    
    if saveOutput == true
        
        cd Output
        save(fileName, 'outputStruct')
        cd ..
        
    end
    
    %% Script loop decrease
    
    scriptCounter = scriptCounter - 1;
end

%% Clean up
if clearWorkspace == true
    % Retain variables used to speed up script or for further analysis
    clearvars -except outputStruct x y stimuli tSeries stimNumber ...
        indivStimuli hrf mx roiCoords
end

toc

%% getGammaHRF
function fun = getGammaHRF(t)

fun = thisGamma(t, 1, 1.0, 0.0, 0.6, 6.0)/100;

end

%% thisGamma
function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
    gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end

gammafun = (amplitude*gammafun+offset);
end

%% Analysis function
function outputStruct = leeAnalysis(params, varargin)

% Unpack params
[tSeries, K, stimuli, hrf, stimNumber,  ...
    voxelSelection, mx, outputStruct] = params{:};

loopCount = 0;

voxels = size(outputStruct.data, 1) * size(outputStruct.data, 2) ...
    * size(outputStruct.data, 3);

if voxelSelection == "roi"
    roiCoords = varargin{1};
elseif voxelSelection == "single"
    voxelCoords = varargin{1};
end

for xC = 1:size(outputStruct.data, 1)
    for yC = 1:size(outputStruct.data, 2)
        for zC = 1:size(outputStruct.data, 3)
            
            tSeriesVoxel = cell(1, stimNumber);
            
            if voxelSelection == "all" || voxelSelection == "roi"
                loopCount = loopCount + 1;
                disp(['(Lee_method) Calculating pRF for voxel ', num2str(loopCount), '/', num2str(voxels)])
            end
            
            for iC = 1:stimNumber
                
                if voxelSelection == "roi"
                    roiX = roiCoords(1, yC);
                    roiY = roiCoords(2, yC);
                    roiZ = roiCoords(3, yC); %fix this to work more efficiently
                    tSeriesVoxel{iC} = squeeze(tSeries{iC}(roiX, roiY, roiZ, :));
                    
                elseif voxelSelection == "single"
                    tSeriesVoxel{iC} = squeeze(tSeries{iC}(voxelCoords(1), ...
                        voxelCoords(2), ...
                        voxelCoords(3), :));
                else
                    tSeriesVoxel{iC} = squeeze(tSeries{iC}(xC, yC, zC, :));
                end
                
                % Remove any time series that contain NaN
                if sum(isnan(tSeriesVoxel{iC}(:))) > 0
                    outputStruct.data{xC, yC, zC} = nan;
                    break
                end
                
            end
            
            if isnan(outputStruct.data{xC, yC, zC})
                continue
            end
            
            %% Regularisation
            % regularisation choices... maybe look at different choices for the
            % regularisation parameter... and compare results./
            
            dataPoints = size(mx, 1) * size(mx, 2);
            
            Bpinv = nan(dataPoints, stimNumber);
            Blasso = nan(dataPoints, stimNumber);
            Bridge = nan(dataPoints, stimNumber);
            Brsvm = nan(dataPoints, stimNumber);
            
            for iC = 1:stimNumber
                % Apply temporal filtering to the time series
                tSeriesVoxel{iC} = highpass(tSeriesVoxel{iC}, 0.01, 0.6667);
                tSeriesVoxel{iC} = lowpass(tSeriesVoxel{iC}, 0.25, 0.6667);
                
                % Linear regression model fitting with different
                % regularization techniques
                tempBpinv = pinv(K{iC}) * tSeriesVoxel{iC};
                tempBlasso = fitrlinear(K{iC}, tSeriesVoxel{iC}, 'Regularization', 'lasso');
                tempBridge = fitrlinear(K{iC}, tSeriesVoxel{iC}, 'Regularization', 'ridge');
                tempBrsvm =  fitrsvm(K{iC}, tSeriesVoxel{iC});
                
                tempBlasso = tempBlasso.Beta;
                tempBridge = tempBridge.Beta;
                tempBrsvm = tempBrsvm.Beta;
                
                % Normalise the time series using the zscore method
                Bpinv(:, iC) = normalize(tempBpinv);
                Blasso(:, iC) = normalize(tempBlasso)';
                Bridge(:, iC) = normalize(tempBridge)';
                Brsvm(:, iC) = normalize(tempBrsvm)';
            end
            
            % Average between the scans
            Bpinv = mean(Bpinv, 2);
            Blasso = mean(Blasso, 2);
            Bridge = mean(Bridge, 2);
            Brsvm = mean(Brsvm, 2);
            
            %% function for reshaping vector -> image
            reshape_estimate = @(x) reshape(x, size(mx));
            
            % Turn the 1d representations back into the shape of stimulus space            
            estimatedPRFpinv = reshape_estimate(Bpinv);
            estimatedPRFlasso = reshape_estimate(Blasso);
            estimatedPRFridge = reshape_estimate(Bridge);
            estimatedPRFsvm = reshape_estimate(Brsvm);
                     
            %% Calculate model prediction
            %  substituting back in to the model
            
            pinvOut = calculate_prf_response(stimuli,...
                estimatedPRFpinv, 'reshape', hrf);
            
            lassoOut = calculate_prf_response(stimuli,...
                estimatedPRFlasso, 'reshape', hrf);
            
            ridgeOut = calculate_prf_response(stimuli,...
                estimatedPRFridge, 'reshape', hrf);
            
            svmOut = calculate_prf_response(stimuli,...
                estimatedPRFsvm, 'reshape', hrf);
            
            %% Z score normalisation and regularisation correction
            % https://stats.stackexchange.com/questions/48045/can-the-bias-introduced-by-lasso-change-the-sign-of-a-coefficient
            % http://statweb.stanford.edu/~tibs/lasso/lasso.pdf
            % check sign of correlation of OLS and LASSO
            
            tSeriesVoxel = cat(1, tSeriesVoxel{:});
            
            dataZ = zscore(tSeriesVoxel);
            
            pinvZ = zscore(pinvOut);
            lassoZ = zscore(lassoOut);
            ridgeZ = zscore(ridgeOut);
            svmZ = zscore(svmOut);
            
            lassoSign = sign(corr(dataZ, lassoZ));
            lassoZ = lassoSign .* lassoZ;
            
            ridgeSign = sign(corr(dataZ, ridgeZ));
            ridgeZ = ridgeSign .* ridgeZ;
            
            %% Output struct construction
            
            estimatedPRFpinv = fliplr(estimatedPRFpinv);
            estimatedPRFlasso =  fliplr(estimatedPRFlasso);
            estimatedPRFridge = fliplr(estimatedPRFridge);
            estimatedPRFsvm = fliplr(estimatedPRFsvm);
            
            outputStruct.data{xC, yC, zC}.Lasso_pRF = estimatedPRFlasso;
            outputStruct.data{xC, yC, zC}.Lasso_best_fit = lassoOut;
            outputStruct.data{xC, yC, zC}.Ridge_pRF = estimatedPRFridge;
            outputStruct.data{xC, yC, zC}.Ridge_best_fit = ridgeOut;
            outputStruct.data{xC, yC, zC}.svm_pRF = estimatedPRFsvm;
            outputStruct.data{xC, yC, zC}.svm_best_fit = svmOut;
            
            outputStruct.data{xC, yC, zC}.PInv_Rsqr = corr(dataZ, pinvZ) ^ 2;
            outputStruct.data{xC, yC, zC}.Lasso_Rsqr = corr(dataZ, lassoZ) ^ 2;
            outputStruct.data{xC, yC, zC}.Ridge_Rsqr = corr(dataZ, ridgeZ) ^ 2;
            outputStruct.data{xC, yC, zC}.svm_Rsqr = corr(dataZ, svmZ) ^ 2;
            
            %% Additional struct construction dependent on analysis type
            
            if voxelSelection == "single"
                outputStruct.Voxel = voxelCoords;
                
                outputStruct.data{1}.pinv_pRF = estimatedPRFpinv;
                outputStruct.data{1}.pinv_best_fit = pinvOut;
                
                outputStruct.data{1}.tSeries = tSeriesVoxel;
                
                outputStruct.data{1}.dataZ = dataZ;
                outputStruct.data{1}.pinvZ = pinvZ;
                outputStruct.data{1}.lassoZ = lassoZ;
                outputStruct.data{1}.ridgeZ = ridgeZ;
                outputStruct.data{1}.svmZ = svmZ;
                
            elseif voxelSelection == "roi"
                outputStruct.data{xC, yC, zC}.Voxel = [roiX, roiY, roiZ];
                
            end
        end
    end
end

end
