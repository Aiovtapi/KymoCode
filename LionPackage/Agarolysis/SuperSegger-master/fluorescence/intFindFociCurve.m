function [ data ] = intFindFociCurve( data, CONST, channelID )
% intFindFociCurve: finds the location of foci and assigns them to the cells.
% The foci are found using the curve (L-) filtered image. The foci are the
% identified using watershed. The foci regions are assigned to the closest cell
% and the location of foci is found in subpixel resolution, using a gaussian fit.
%
% INPUT : 
%       data : cell/regions file (err file)
%       CONST : segmentation constants
%       channelID : fluorescence channel number
%
% OUTPUT : 
%       data : updated data with sub-pixel fluorescence model
%       a locus[channelID] field is added with the following structure :
%              r : sub-pixel global coordinates for foci location (x y)
%              fitSigma : sigma of gaussian fit
%              intensity : maximum intensity of gaussian fit
%              fitScore : score of gaussian fitting to the foci (0 - 1)
%              normIntensity : normalized intensity (divided by std(cell
%              fluor)
%              score : locus intensity / std(cell fluorescence) *  gaussian fit score;
%              shortaxis : sub-pixel local coordinates for foci location (x y)
%              longaxis : sub-pixel local coordinates for foci location (x y)
%
% Copyright (C) 2016 Wiggins Lab 
% Written by Connor Brennan, Stella Stylianidou & Paul Wiggins.
% University of Washington, 2016
% This file is part of SuperSegger.
% 
% SuperSegger is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SuperSegger is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SuperSegger.  If not, see <http://www.gnu.org/licenses/>.

DEBUG_FLAG = false;
fieldname = ['locus', num2str(channelID)];
options =  optimset('MaxIter', 1000, 'Display', 'off', 'TolX', 1/10);

% Get images out of the structures.
originalImage = double(data.(['fluor',num2str(channelID)]));

% Subtract gaussian blurred image 
% to get rid of big structure
hg = fspecial( 'gaussian' , 210, 30 );
highPassImage = originalImage - imfilter( originalImage, hg, 'replicate' );
cytoplasmicFlourescenceSTD = std(double(highPassImage(data.mask_bg)));
normalizedImage = originalImage/cytoplasmicFlourescenceSTD; % normalization so that intensities;

%Take only pixels above 1 std (noise reduction?)
normalizedImage = normalizedImage - 1;
normalizedImage(normalizedImage < 0) = 0; % logical mask of foci found

filter_width = 1.5;
[~,~,fluorFiltered] = curveFilter(normalizedImage,filter_width);
data.(['fluor',num2str(channelID),'_filtered']) = fluorFiltered;

mask_mod = bwmorph (data.mask_bg, 'dilate', 1);

fociWatershed = watershed(-fluorFiltered); % watershed to identify foci
maskedFociWatershed = logical(double(fociWatershed).*double(mask_mod));

fociRegionLabels = bwlabel(maskedFociWatershed);
props = regionprops( fociRegionLabels, {'BoundingBox'} );
numFociRegions = numel(props);
imsize = size(highPassImage);

% initialize focus fields
focusInit.r = [nan,nan];
focusInit.score = nan;
focusInit.intensity = nan;
focusInit.normIntensity = nan;
focusInit.shortaxis = nan;
focusInit.longaxis = nan;
focusInit.fitSigma = nan;
focusInit.fitScore = nan;

fociData = [];
cellIDs = [];

if DEBUG_FLAG
    figure(2);
    clf;
end

for ii = 1:numFociRegions
    tempData = focusInit;
    
    [xPad,yPad] = getBBpad( props(ii).BoundingBox, imsize, 3 );
    [meshX,meshY] = meshgrid(xPad, yPad);
    
    maskToFit = (fociRegionLabels(yPad, xPad) == ii); % foci region
    imageToFit  = fluorFiltered(yPad, xPad);  % filtered image 
    imageToFit = imageToFit .* double(maskToFit);
    
    [~, maxIndex] = max(imageToFit(maskToFit));
    
    tempImage = imageToFit(maskToFit);
    tempData.intensity = tempImage(maxIndex);
   
    tempImage = meshX(maskToFit);
    fociX = tempImage(maxIndex);
    tempData.r(1) = fociX;
   
    tempImage = meshY(maskToFit);
    fociY = tempImage(maxIndex);
    tempData.r(2) = fociY;

    % figure out which cell the focus belongs to
    maskSize = [numel(yPad),numel(xPad)];

   
    cellsLabel = data.regs.regs_label(yPad,xPad);
    cellsMask = logical(cellsLabel);
    tempMask = zeros(maskSize);
    
    tempMask(fociY - yPad(1)+1, fociX - xPad(1)+1 ) = 1;
    distanceToFoci = bwdist( tempMask );

    cellIDList = cellsLabel(cellsMask);
    [~, minDistanceIndex] = min(distanceToFoci(cellsMask));
    bestCellID = cellIDList(minDistanceIndex);
    
    if ~isempty( bestCellID )
        croppedImage = highPassImage(data.CellA{bestCellID}.yy, data.CellA{bestCellID}.xx);
        cellFlouresenseSTD = std(croppedImage(data.CellA{bestCellID}.mask));  
        
        if tempData.intensity >0
            %Initialize parameters
            backgroundIntensity = 0;
            gaussianIntensity = fluorFiltered(fociY, fociX) - backgroundIntensity;
            sigmaValue = 1;

            parameters(1) = fociX;
            parameters(2) = fociY;
            parameters(3) = gaussianIntensity;
            parameters(4) = sigmaValue;
            %parameters(5) = backgroundIntensity;

            [parameters] = fminsearch( @doFit, parameters, options);

            gaussianApproximation = makeGassianTestImage(meshX, meshY, parameters(1), parameters(2), parameters(3), backgroundIntensity, parameters(4));

            %Crop out fit gaussian from original image
            croppedImage = imageToFit;
            croppedImage(gaussianApproximation < 0.1 * max(max(gaussianApproximation))) = 0;

            imageTotal = sqrt(sum(sum(croppedImage)));
            guassianTotal = sqrt(sum(sum(gaussianApproximation)));

            fitScore = sum(sum(sqrt(croppedImage) .* sqrt(gaussianApproximation))) / (imageTotal * guassianTotal);

            if DEBUG_FLAG
                figure(1);
                clf;
                subplot(2, 2, 1);
                imshow(imageToFit, []);
                subplot(2, 2, 2);        
                imshow(croppedImage, []);    
                subplot(2, 2, 3);        
                imshow(gaussianApproximation, []);   
                subplot(2, 2, 4);          
                title(['Score: ', num2str(fitScore)]);
                keyboard;
            end

            tempData.r(1) = parameters(1);
            tempData.r(2) = parameters(2);
            tempData.fitSigma = parameters(4);
            tempData.intensity = parameters(3);
            tempData.fitScore = fitScore;

            %Calculate scores        
            tempData.normIntensity = tempData.intensity;
            tempData.score = tempData.intensity / (sqrt(pi)*tempData.fitSigma) * tempData.fitScore;
            tempData.shortaxis = ...
                (tempData.r-data.CellA{bestCellID}.coord.rcm)*data.CellA{bestCellID}.coord.e2;
            tempData.longaxis = ...
                (tempData.r-data.CellA{bestCellID}.coord.rcm)*data.CellA{bestCellID}.coord.e1;


            %Assign to array
            cellIDs(ii) = bestCellID;  
            focusData(ii) = tempData;
            if DEBUG_FLAG
               figure(2);
               hold on;
               plot(fociX, fociY, '.r' );
               text(fociX, fociY, num2str( tempData.intensity, '%1.2g' ));
            end
        end
    end
end


% assign to cells
for ii = 1:data.regs.num_regs

    fociIndex = find(cellIDs == ii);
    tempData = focusData(fociIndex);
    
    [~, order] = sort( [tempData.intensity], 'descend' );
    sortedFoci = tempData(order);
    
    focus = focusInit;
    if numel(sortedFoci) > 0
        maxIndex = find([sortedFoci.intensity] > 0.333 * sortedFoci(1).intensity);
        if numel(maxIndex) > CONST.trackLoci.numSpots(channelID)
            maxIndex = maxIndex(1:CONST.trackLoci.numSpots(channelID));
        end

        numFoci = numel(maxIndex);

        for jj = 1:numFoci
            focus(jj) = sortedFoci(jj);
        end 
    end
    
    scores = [focus(:).score];
    focus = focus(~isnan(scores));  
    data.CellA{ii}.(fieldname) = focus;
    
    % creates a filtered image of the cell
    xPad = data.CellA{ii}.xx;
    yPad = data.CellA{ii}.yy;
    data.CellA{ii}.(['fluor',num2str(channelID),'_filtered'])=fluorFiltered( yPad, xPad );
end

    function error = doFit(parameters)
        % doFit : does the gaussian fit to the foci and calculates the
        % error.
        %parameters are : 
        %parameter(1) - Sub-pixel resolution of foci position X
        %parameter(2) - Sub-pixel resolution of foci position Y
        %parameter(3) - Intensity of the gaussian
        %parameter(4) - sigma of gaussian
        %parameter(5) - background intensity
        gaussian = makeGassianTestImage(meshX, meshY, parameters(1), parameters(2), parameters(3), backgroundIntensity, parameters(4));
        tempImage = (double(imageToFit) - gaussian);
        error = sum(sum(tempImage.^2));
    end

    function testImage = makeGassianTestImage(meshX, meshY, fociX, fociY, gaussianIntensity, backgroundIntensity, sigmaValue)
        testImage = backgroundIntensity + gaussianIntensity * exp( -((meshX - fociX).^2 + (meshY - fociY).^2)/(2 * sigmaValue^2) );
    end

end