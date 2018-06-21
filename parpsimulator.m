classdef parpsimulator
    %PARPSIMULATOR  Class to simulate free diffusion of PARP molecules
    %
    % parpsimulator Properties:
    %   pxSize - Physical units of the simulated grid in microns/pixel
    %   deltaT - Physical units of each time step in seconds/time step
    %   diffusionCoeff - Diffusion coefficient in microns^2/seconds
    %   mparpFrac - Percent fraction of mobile PARP particles
    %   unbleachedFrac - Percent fraction of photobleached particles
    %   numParticles - Number of particles to simulate
    %   numSteps - Number of time steps to simulate
    %   outputMovie - If true, a movie will be saved
    %
    % parpsimulator Methods:
    %   simulate - Run the simulation using specified parameters
    %   parameterSweep - Search for best fit diffusion coefficient and mobile PARP fraction
    %
    % Example:
    %   %Create an instance of the class
    %   sim = parpsimulator
    %
    %   %Run the simulation with default settings
    %   
    
    properties
        
        %Simulation parameters
        pxSize = 0.08677;      %Physical size of grid in microns/pixel
        deltaT = 0.19;         %Physical size of time step seconds/time step
        
        %Parameters to fit
        meanStepSize = 16.80;  %Mean particle step size in pixels/time step
        mparpFrac = 25;        %Percent fraction of mobile PARP
        unbleachedFrac = 46;     %Percent fraction of photobleached particles
        
        %Simulation parameters
        numParticles = 12000;  %Number of simulated particles
        numSteps = 500;        %Number of simulation time steps
        particleMotion = 'lattice'; %Type of particle motion ('2D' or 'circular')
        initialPositions = 'random'; %Initial position of particles
        
        outputMovie = false;   %If true, an AVI of the simulation will be saved      
        
    end
    
    methods
        
        function output = example(obj, varargin)
            %EXAMPLE  Run an example simulation
            %
            %  I = EXAMPLE(OBJ) runs the simulation using an example
            %  dataset. 
            %
            %  Example:
            %    %Create a new parpsimulator object
            %    S = parpsimulator;
            %
            %    %Run the example script
            %    I = example(S);
            %
            %    %Plot the resulting intensity in the trapped region
            %    %(normalized to the first time simulated timepoint)
            %    plot(I.Time, I.ROI0./I.ROI0(1))
            %
            %  I = EXAMPLE(OBJ, TYPE) runs the simulation using other
            %  test scenarios to validate the model. TYPE can be:
            %    'nucleus' - (Default) Mask is from a nucleus
            %    'circle' - Mask is a circle
            %    'point' - Simulate diffusion from a single point
            %     
            
            if isempty(varargin)
                simType = 'nucleus';
            else
                simType = varargin{1};
            end
            
            switch lower(simType)
                
                case 'nucleus'
                    %Load test mask and trap ROI
                    mask = 'testdata\2.15.18_2min_GFPP1_Hela_005NuclMask.txt';
                    roi = 'testdata\2.15.18_2min_GFPP1_Hela_005ROI.txt';
                    
                    output = simulate(obj, mask, roi, 'verbose');
                    
                case 'point'
                    
                    mask = true(1024);
                    roi = false(1024);
                    
                    S = parpsimulator;
                    S.numParticles = 100;
                    S.meanStepSize = 2;
                    S.unbleachedFrac = 0;
                    S.mparpFrac = 100;
                    S.initialPositions = 'center';
                    
                    output = simulate(S, mask, roi);
                    
            end
            
        end
        
        function varargout = simulate(obj, mask, roi, varargin)
            %SIMULATE  Run the free diffusion simulation
            %
            %  S = SIMULATE(OBJ, MASK, ROI) will run the free diffusion
            %  simulation using the current simulation parameters, as
            %  defined in the object properties. MASK specifies the cell
            %  nucleus and ROI specifies the trapping/irradiation region.
            %  The simulation grid is set to be the same dimensions as
            %  MASK.
            % 
            %  Simulation data is returned as a struct S, which has the
            %  following fields: 'Time' the simulation time, 'ROI0' the
            %  intensity in the trapping region, 'ParticleData' which
            %  contains particle positions, and 'isMobile' which specified
            %  whether the particle is mobile or stationary. 'ParticleData'
            %  is an N-by-2-by-T matrix representing the [X, Y] coordinates
            %  of the particles at each simulated time step. N is the
            %  number particles and T is the total number of timesteps in
            %  the simulation.
            %
            %  SIMULATE(OBJ, ..., FILENAME) will save the simulated data as
            %  both a MAT-file containing all particle locations, and a CSV
            %  file with the simulation time and number of particles in
            %  ROI.
            %
            %  The diffusion of PARP molecules is simulated using a
            %  random walk. Particle movement is constrained within the
            %  cell nucleus, as well as within the 'trap region'. Particle
            %  movement can be constained along a lattice or allowed to
            %  diffuse freely in a circle around its current position.
            %  
            %  The cell nucleus should be defined using MASK. MASK can
            %  either be the char array with the path to a text file with
            %  the nuclear parameter coordinates, or a logical matrix
            %  containing the nuclear mask. If a text file is used, the X
            %  and Y coordinates of each point along the nuclear boundary
            %  should be specified, with new points in individual rows.
            %  Tracing the nuclear boundary can be done for example by
            %  using 'bwtraceboundary'.
            %
            %  The irradiation region is specified by ROI. ROI can either
            %  be a char array with the path to a text file containing the
            %  ROI coordinates [top, left, width, height] or a logical
            %  matrix containing the ROI mask. This region will be used to
            %  compute the intensity increase, as well as the trapping
            %  region. The trapping region is ROI expanded along its
            %  longest dimension to fill the nuclear mask.
            %
            %  See also: bwtraceboundary
            
            %--Output data structure--%
            %  D.ROI0
            %  D.Time
            %
            %  D.ParticleData (rows = particle, columns = [X, Y], z = time)
            %  D.isMobile
            
            %Parse optional inputs
            outputFN = '';
            isVerbose = true;
            
            iArg = 1;
            while iArg <= numel(varargin)
                
                switch lower(varargin{iArg})
                    
                    case 'verbose'
                        isVerbose = true;   %Display information about the simulation after running
                        
                    case 'silent'
                        isVerbose = false;  %Do not display information
                        
                    otherwise
                        %Parse filename
                        if ischar(varargin{iArg})
                            %Parse the directory and filename
                            [fpath, fname] = fileparts(varargin{iArg});
                            
                            %Check that directory exists, if not make it
                            if ~isempty(fpath) && ~exist(fpath, 'dir')
                                mkdir(fpath);
                            end
                            
                            outputFN = fullfile(fpath, fname);
                            
                        else
                             error('parpsimulator:simulate:InvalidParameter',...
                                'Expected parameter name to be a char array.')
                        end
                end
                iArg = iArg + 1;
            end
            
            nuclMask = parpsimulator.loadMask(mask);
            
            [roiMask, trapROI] = parpsimulator.loadROI(roi);
            
            %The simulation computes data units of pixel and time step.
            %Conversion between the discretized to physical units occurs
            %before returning the final data.
                        
            %Initialize storage matrix for output data
            dataOut.Time = ((1:obj.numSteps) * obj.deltaT);
            dataOut.ROI0 = zeros(1, obj.numSteps);
            
            %--- Start simulation ---%
            tStart = tic;  %Start the timer
            
            %---Time step 1: Initialization---%
            
            %Generate a set of particles within the cell
            switch lower(obj.initialPositions)
                
                case 'random'
                    %Select particle positions randomly within the nuclear mask
                    particleXY = parpsimulator.randomPoints(nuclMask, obj.numParticles);
                    
                case 'center'
                    %Set particle positions to the center of the image
                    particleXY = repmat(size(nuclMask)/2,obj.numParticles, 1);
                    
            end
            
            %Determine if the particle is mobile or stationary
            mInd = randsample(obj.numParticles, round(obj.numParticles .* (obj.mparpFrac/100)));
            isMobile = false(obj.numParticles, 1);
            isMobile(mInd) = true;
            %isMobile = rand(obj.numParticles, 1) <= (obj.mparpFrac/100);
            
            %Determine if particle is in trapped region or not
            isTrapped = parpsimulator.inROI(trapROI, particleXY);
            
            %Record number of particles in trapped region initially
            numTrappedInit = nnz(isTrapped);
            dataOut.ROI0(1) = nnz(isTrapped);
            
            %For particles originally in the trapped region, randomly
            %remove a number of particles which are bleached
            bleachedInd = randsample(find(isTrapped), floor(nnz(isTrapped) * (1 - obj.unbleachedFrac/100)));
            
%             isBleached = false(obj.numParticles,1);
%             isBleached(bleachedInd) = true;
            
            %For proteins in the trapped region, simulate bleaching by removing a
            %fraction of them
            particleXY(bleachedInd,:) = [];
            isMobile(bleachedInd) = [];
            isTrapped(bleachedInd) = [];
            
            numRemainingParticles = size(particleXY,1);
            
            %Store particle positions
            dataOut.ParticlePos = zeros(size(particleXY,1), 2, obj.numSteps);
            dataOut.ParticlePos(:,:,1) = particleXY;
            dataOut.isMobile = isMobile;            
            
            %---Timestep 2:N---%
            for iStep = 2:obj.numSteps
                
                %Move each particle according to the step statistics. The
                %particles are assumed to be able to move about a circle
                %around their current position. The step size is determined
                %by a normal distribution with a mean determined by the
                %diffusion coefficient, and a standard deviation of 0.001.
                switch lower(obj.particleMotion)
                    
                    case {'circle', 'continuous'}
                        mvmt = (randn(numRemainingParticles, 2) + obj.meanStepSize) .* isMobile;
                        angle = rand(numRemainingParticles, 2) * 2 * pi;
                        newXY = particleXY + [mvmt(:, 1) .* cos(angle(:, 1)), mvmt(:, 2) .* sin(angle(:,2))];
                        
                    case {'2d', 'lattice'}
                        
                        mvmt = [randsample([-1, 0, 1], numRemainingParticles, true)', randsample([-1, 0, 1], numRemainingParticles, true)'] .* (randn(numRemainingParticles, 2) * 0.001 + obj.meanStepSize) .* isMobile;
                        newXY = particleXY + mvmt;
                        
                end
                
                %Check if new step is allowed (i.e. within the cell mask and not migrating out of the trapped region)
                validStep = parpsimulator.inROI(nuclMask, newXY) & ~isTrapped;
                
                %Update the particle position if new position is valid
                particleXY(validStep,:) = newXY(validStep,:);
                
                %Update trapped status
                isTrapped = parpsimulator.inROI(trapROI, particleXY);
                
                %Compute number of particles in the trapped region
                dataOut.ROI0(iStep) = nnz(isTrapped);
                
                %Save particle position to output data structure
                dataOut.ParticlePos(:,:,iStep) = particleXY;
                
            end
            
            %---Print simulation statistics---
            
            if isVerbose
                toc(tStart); %Stop timer and display time taken to simulate
                
                %Compute the step size from the diffusion coefficient
                %diffConstInPxPerTime = (obj.diffusionCoeff / (obj.pxSize)^2) * obj.deltaT;
                scaledD = parpsimulator.calcScaledDiffCoeff(obj.meanStepSize, obj.pxSize, obj.deltaT);
                
                fprintf('Effective diffusion constant = %.3g micron^2/s\n', scaledD);
                
                disp(['Actual mobile fraction = ', num2str(numel(mInd)/(obj.numParticles))])
                
                %Display number of bleached particles
                disp(['Actual bleached fraction = ', num2str(numel(bleachedInd)/(numTrappedInit))])
            end
            
            %---Output movie---%
            %Save data
            if ~isempty(outputFN)
               
                %Save movie
                if obj.outputMovie
                    parpsimulator.exportAVI([outputFN, '.avi'], nuclMask, roiMask, dataOut)
                end
                
                %Save MAT-file
                dataOut.ParticlePosInMicrons = dataOut.ParticlePos .* obj.pxSize;
                save([outputFN,'.mat'], 'dataOut');
                
                %Save simulation parameters
                fid = fopen([outputFN, '_params.txt'], 'w');
                
                objProp = properties(obj);
                for iProp = 1:numel(objProp)
                    fprintf(fid, '%s = %f\n', objProp{iProp}, obj.(objProp{iProp}));
                end
                fclose(fid);
                
                %Save output time and intensity as CSV
                fid = fopen([outputFN, '_Intensity.txt'], 'w');
                for iS = 1:obj.numSteps
                    fprintf(fid, '%f, %f\n', dataOut.Time(iS), dataOut.ROI0(iS));
                end                
                fclose(fid);
                
                fprintf('Data saved as %s\n', outputFN);
            end
            
            if nargout > 0
                varargout{1} = dataOut;
            end
            
        end
        
        function [bestFit, Rsq, parpRAMP, stepSizeRAMP] = parameterSweep(obj, mask, roi, inputData, varargin)
            %PARAMETERSWEEP  Find parameters to fit simulation to data
            %
            %  S = PARAMETERSWEEP(OBJ, MASK, ROI, TESTDATA) will compute
            %  the so-called fitness landscape of the simulation compared
            %  to TESTDATA, by varying the mobile PARP fraction and the
            %  mean step size parameters. The best fit parameters and data
            %  are returned as a struct S (i.e. S.mparpFrac and
            %  S.meanStepSize).
            %
            %  By default, the PARP fraction is varied from 1% to 30% in
            %  steps of 1%, and the mean step size is varied from 1 to 20
            %  pixels. MASK and ROI should be the nuclear mask and ROI of
            %  the irradiation region as used in the 'simulate' function.
            %
            %  To improve the accuracy of the final parameter, it could be
            %  helpful to carry out two sweeps: first, a rough sweep should
            %  be carried out to identify approximate values for the mobile
            %  PARP fraction and the mean step size. A second, finer sweep
            %  should then carried out the identify the best-fit values
            %  within the tolerance specified.
            %
            %  S = PARAMETERSWEEP(..., 'Parameter', Value) allows the
            %  default range and sweep step size to be changed. The
            %  following lists the parameters available with default values
            %  in paratheses:
            %
            %     offsetFrames - number of pre-irradiation frames (6)
            %     mparpfracrange - range of mobile PARP fraction ([10, 35])
            %     mparpprecision - step size of mobile PARP fraction (1)
            %     stepsizerange - range of mean step size values ([8 15])
            %     stepsizeprecision - precision of step sizes (0.5)
            %     showplots - true if plots should be shown (true)
            %
            %  [S, R, X, Y] = PARAMETERSWEEP(...) will also return the
            %  matrix of r-squared values in R, the mobile PARP fraction
            %  values in X and the mean step sizes in Y.
            %
            %  See also: parpsimulator.simulate
            
            tic;
            %Parse the variable inputs
            numPreIrradFrames = 6;
            mparpFracRange = [10, 35];
            mparpSweepSize = 1;
            meanStepSizeRange = [8, 15];
            meanStepSizeSize = 0.5;
            showPlots = true;
            
            iP = 1;
            while iP <= numel(varargin)
                
                switch lower(varargin{iP})
                    
                    case 'offsetframes'
                        
                        numPreIrradFrames = varargin{iP + 1};
                        
                    case 'mparpfracrange'
                        
                        mparpFracRange = varargin{iP + 1};
                        
                    case 'mparpprecision'
                        
                        mparpSweepSize = varargin{iP + 1};
                        
                    case 'stepsizerange'
                        meanStepSizeRange = varargin{iP + 1};
                        
                    case 'stepsizeprecision'
                        meanStepSizeSize = varargin{iP + 1};
                        
                    case 'showplots'
                        showPlots = varargin{iP + 1};
                        
                end
                
                iP = iP + 2;
            end
            
            
            
            %Load the data to be fitted
            if ischar(inputData)
                
                data = csvread(inputData, 4, 0);
                fitTime = data(:,1);
                fitData = data(:,12);

                %Normalize the data
                fitData = fitData./fitData(1);
                
            else
                %Check that the input data has two columns
                if size(inputData,2) ~= 2
                    error('parpsimulator:parameterSweep:IncorrectSizeInput',...
                        'Expected input data to have two columns: for time and intensity.')
                end
                
                fitTime = inputData(:,1);
                fitData = inputData(:,2);
            end
             
            %Determine length of time to simulate to match input time
            obj.numSteps = ceil(fitTime(end) ./ obj.deltaT);
                        
            %Initial rough sweep
            parpRAMP = mparpFracRange(1):mparpSweepSize:mparpFracRange(2);
            stepSizeRAMP = meanStepSizeRange(1):meanStepSizeSize:meanStepSizeRange(2);
            
            %Initialize output variable
            Rsq = zeros(numel(parpRAMP), numel(stepSizeRAMP));
            simResults = cell(numel(parpRAMP), numel(stepSizeRAMP));
            
            %Perform the sweep
            for iMPARP = 1:numel(parpRAMP)
                for iMeanStepSize = 1:numel(stepSizeRAMP)
                    
                    obj.mparpFrac = parpRAMP(iMPARP);
                    obj.meanStepSize = stepSizeRAMP(iMeanStepSize);
                    
                    I = simulate(obj, mask, roi, 'silent');
                    
                    simTime = [0; I.Time' + fitTime(numPreIrradFrames)];
                    simInt = [1; (I.ROI0./I.ROI0(1))'];
                    
                    simResults{iMPARP, iMeanStepSize} = [simTime simInt];
                    
                    Rsq(iMPARP, iMeanStepSize) = parpsimulator.computeError([simTime simInt], [fitTime fitData]);
                    
                end
            end
            
            %Find the best fit (i.e. the r-squared value closest to 1)
            [~, bestInd] = min(abs(Rsq(:) - 1));
            
            [fitX, fitY] = ind2sub(size(Rsq), bestInd);
            
            bestFit.mparpFrac = parpRAMP(fitX);
            bestFit.meanStepSize = stepSizeRAMP(fitY);
            bestFit.equivDiffCoeff = parpsimulator.calcScaledDiffCoeff(bestFit.meanStepSize, obj.pxSize, obj.deltaT);
            bestFit.Rsq = Rsq(bestInd);
            bestFit.simTime = simResults{fitX, fitY}(:,1);
            bestFit.simIntensity = simResults{fitX, fitY}(:,2);
            
            toc;
            
            if showPlots
                
                %Make a best fit plot
                figure;
                plot(fitTime, fitData, 'bo', simResults{fitX, fitY}(:,1), simResults{fitX, fitY}(:,2), 'r')
                xlabel('Time (s)');
                ylabel('Normalized Intensity');
                legend('Experimental data', 'Simulated fit')
                
                figure;
                %Plot the R-squared values
                pcolor(stepSizeRAMP,parpRAMP,Rsq)
                xlabel('Mean step size(\mum)')
                ylabel('Mobile PARP fraction (%)')
                shading flat
                grid off
                
            end
        end
        
    end
    
    methods (Static, Hidden)
        
        function exportAVI(filename, nuclMask, roiMask, dataStruct)
            %DRAWPARTICLES  Draw particles to generate a movie
            %
            %  IMG = DRAWPARTICLES(MASK, PARTICLEDATA)
                        
            vidOut = VideoWriter(filename);
            vidOut.FrameRate = 25;
            vidOut.Quality = 100;
            open(vidOut);
            
            %Generate the cell image - the nuclear outline will be grey and
            %the background will be white. The outline is thickened to be
            %more distinct.
            nuclImage = uint8(imdilate(bwperim(nuclMask),strel('diamond', 1)));
            nuclImage(nuclImage == 0) = 255;
            nuclImage = nuclImage .* 120;            
            
            %Add the trap ROI
            nuclImage(imdilate(bwperim(roiMask), strel('diamond', 1))) = 0;
            nuclImage = repmat(nuclImage, 1, 1, 3);
            
            for frame = 1:size(dataStruct.ParticlePos,3)
                
                %Plot the non-mobile fraction
                imgOut = insertShape(nuclImage, 'filledcircle', [dataStruct.ParticlePos(~dataStruct.isMobile,:,frame), ones(nnz(~dataStruct.isMobile),1) .* 0.5], 'color', [180 180 180]);
                
                %Plot the mobile fraction
                imgOut = insertShape(imgOut, 'filledcircle', [dataStruct.ParticlePos(dataStruct.isMobile,:,frame), ones(nnz(dataStruct.isMobile),1) .* 0.5], 'color', 'blue');
                              
                %Insert time
                imgOut = insertText(imgOut, [size(nuclImage, 2)/2, size(nuclImage, 1)] - [0, 30], sprintf('T = %0.2f s', dataStruct.Time(frame)), 'AnchorPoint', 'Center', 'BoxOpacity', 0, 'TextColor', 'black');
                
                %Write frame to AVI file
                writeVideo(vidOut, imgOut);

            end
            close(vidOut);
            
        end
        
        function isInROI = inROI(mask, points)
            %INROI  Returns true if the specified point is in the region of interest
            %
            %  L = INROI(MASK, POINT) returns a value of true if the specified point
            %  is located in the MASK. POINT should be a 1-by-2 vector with the XY
            %  coordinates of the point. Note that POINT is specified in pixels.
            %
            %  L = INROI(MASK, POINTS) returns a logical vector with the same number of
            %  elements as rows in POINTS. POINTS should be an N-by-2 matrix specifying
            %  a series of points, with each row corresponding to a different location.
            %
            %  Example:
            %    %Generate linear indices for the image
            %    XX = 1:100;
            %    [XX, YY] = meshgrid(XX, XX);
            %
            %    %Generate a circular mask, with a radius of 8
            %    mask = false(100);
            %    mask((XX - 50).^2 + (YY - 50).^2 <= 8^2) = true;
            %
            %    %Select a point inside and outside the circle
            %    points = [50, 50; 1, 1];
            %
            %    %Plot the mask and selected point
            %    imshow(mask)
            %    hold on
            %    plot(points(:,2), points(:,1), 'x')
            %
            %    %Check whether the points are in the ROI
            %    L = INROI(mask, points);
            
            points = round(points);
            
            ind = sub2ind(size(mask), points(:,2), points(:,1));
            
            isInROI = mask(ind);
            
        end
        
        function XYout = randomPoints(mask, numPoints)
            %RANDOMPOINTS  Generate a set of randomly chosen points within a region
            %
            %  XY = RANDOMPOINTS(MASK, NUMPTS) generates a set of randomly chosen
            %  points within the region specified by the MASK. MASK should be a logical
            %  matrix with a value of true (1) for the region of interest.
            %
            %  The function works by randomly selecting points within the horizontal
            %  and vertical extent of the MASK. The code then checks and accepts points
            %  within the MASK.
            %
            %  Example:
            %    %Generate the image grid
            %    XX = linspace(-10, 10, 100);
            %    [XX,YY] = meshgrid(XX, XX);
            %
            %    %Generate a circular mask, with a radius of 8
            %    mask = false(100);
            %    mask(XX.^2 + YY.^2 <= 8^2) = true;
            %
            %    %Select 200 random points within the circle
            %    pts = RANDOMPOINTS(mask, 200);
            %
            %    %Plot the mask and show the selected points
            %    pcolor(mask)
            %    shading flat; grid off; axis image
            %    hold on
            %    plot(pts(:, 1), pts(:, 2), 'r.')
            %
            %  Note: The function uses the 'rand' function. To make the process
            %  repeatable, use the built-in MATLAB functions 'rng' to set the seed and
            %  generator type.
            %
            %  See also: rng
            
            %Validate the inputs
            if ~islogical(mask)
                %Try to binarize the mask
                temp = false(size(mask));
                temp(mask > 0) = 1;
                mask = temp;
            end
            
            %Compress the mask to find limits of random number generation
            maskHorz = find(any(mask, 1));
            maskVert = find(any(mask, 2));
            
            xLim = [min(maskHorz), max(maskHorz)];
            yLim = [min(maskVert), max(maskVert)];
            
            %Generate points
            XYout = zeros(numPoints, 2);  %Storage matrix
            generatedPts = 0;
            
            while generatedPts < numPoints
                
                %Compute number of poitns to generate
                nPtsToGenerate = numPoints - generatedPts;
                
                randPts = [rand(nPtsToGenerate,1) * (xLim(2) - xLim(1)) + xLim(1), ...
                    rand(nPtsToGenerate, 1) * (yLim(2) - yLim(1)) + yLim(1)];
                
                %Accept points that are within the mask
                chkPt = round(randPts);
                chkInd = sub2ind(size(mask), chkPt(:,2), chkPt(:,1));
                ptIsValid = mask(chkInd);
                XYout((generatedPts + 1):(generatedPts + nnz(ptIsValid)),:) = randPts(ptIsValid,:);
                
                %Update number of generated points
                generatedPts = generatedPts + nnz(ptIsValid);
                
            end
            
        end
        
        function scaledD = calcScaledDiffCoeff(meanStepSize, pxSize, timeStep)
            %CALCSCALEDDIFFCOEFF  Calculate scaled diffusion coefficient
            %
            %  SCALED_D = CALCSCALEDDIFFCOEFF(meanStepSize, pxSize,
            %  timestep) computes the scaled diffusion coefficient from the
            %  mean step size (in pixels), the pixel size (in
            %  microns/pixel) and the size of the timestep (in seconds).
            %  The equation is 
            %
            %    scaledD = (1/3) * meanStepSize^2 * pxSize^2 * (1/ timeStep)
            
            scaledD = (1/3) .* meanStepSize.^2 .* pxSize.^2 .* (1./ timeStep);
            
        end
        
        function rsquared = computeError(testData, expData)
            %COMPUTEERROR  Compute the R-square error
            %
            %  R = COMPUTEERROR(TESTDATA, EXPDATA) computes the R-squared
            %  coefficient between the TESTDATA and EXPDATA. TESTDATA and
            %  EXPDATA should be a N-by-2 matrix, with column 1 being the
            %  timepoints and column 2 being the intensity.
            
            %Interpolate to match timepoints
            interpROI0 = interp1(testData(:,1), testData(:,2), expData(:,1), 'linear', 'extrap');
            
            %Normalize the data to the first point
            interpROI0 = interpROI0./interpROI0(1);
            
            %Compute R-squared error
            SStot = sum((expData(:,2) - mean(expData(:,2))).^2);
            SSres = sum((expData(:,2) - interpROI0).^2);
            
            rsquared = 1 - SSres./SStot;
            
        end
        
        function maskOut = loadMask(maskFile)
            
            %Parse inputs
            if ischar(maskFile)
                %Load the nuclear mask
                maskPerim = dlmread(maskFile,',');
                maskOut = false(512);
                maskOut(sub2ind([512 512], maskPerim(:,1)', maskPerim(:,2)')) = true;
                maskOut = imfill(maskOut,'holes');
            else
                maskOut = maskFile;
            end
            
        end
        
        function [roiMask, trapROI] = loadROI(roiFile)
            
            if ischar(roiFile)
                %Load the ROI coordinates
                roiCoords = dlmread(roiFile,',');
                roiMask = false(512);
                roiMask(roiCoords(2):(roiCoords(2) + roiCoords(4)), roiCoords(1):(roiCoords(1) + roiCoords(3))) = true;
                
                %Expand the trap ROI to fill the image
                trapROI = false(size(roiMask));
                
                if roiCoords(4) > roiCoords(3)
                    trapROI(1:end, roiCoords(1):(roiCoords(1) + roiCoords(3))) = true;
                else
                    trapROI(roiCoords(2):(roiCoords(2) + roiCoords(4)), 1:end) = true;
                end
                
            else
                roiMask = roiFile;
                trapROI = roiFile;
            end
            
        end
    end
        
end
