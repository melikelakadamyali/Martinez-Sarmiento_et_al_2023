% oniRawBincropper.m
% Function to crop the a tiff file and overlay the corresponding bin file
% outputs a .Dax and a .bin file
% outputs a file for each cropped region

% scale parameters for the plot object
pxtonm = 117; % scale for nanometers to pixel units
upscale = 10;  % upscale factor of super res data to image

% load a tif file and upscale it
[tifName, tifLoc, ~] = uigetfile('.tif');
[~,~,ext] = fileparts(fullfile(tifLoc,tifName));
ext_len = length(ext);
tifPath = fullfile(tifLoc,tifName);
tifData = Tiff(tifPath,'r');
imageData = read(tifData);
osz = size(imageData);
bigImage = imresize(imageData,upscale);
imsz = size(bigImage); 
bigImage = bigImage - min(min(bigImage));
bgIThresh = prctile(bigImage(bigImage>0),100);
bigImage(bigImage>bgIThresh) = bgIThresh;

% % load a bin file
% [binName, binLoc, ~] = uigetfile('.bin');
% LL = Insight3( fullfile( binLoc,binName ) );
% yx = LL.getXYcorr;
% dat = LL.getData();
% yxBig = yx*upscale;
% yxPix = ceil(yxBig);
% load a csv file
[csvName, csvLoc, ~] = uigetfile('.csv');
q=csvread(fullfile(csvLoc,csvName),1,0);
yx = q(:,[4,3])/pxtonm;
yxBig = yx*upscale;
yxPix = ceil(yxBig);

% create the image from scaled coordinates (skip neg coordinates)
Img = zeros(imsz+1);
for ii = 1:size(yxPix,1)
    if any(yxPix(ii,:)<=1) || any(yxPix(ii,:)>imsz)
        continue
    end
    % simulate a rough gaussian (only 1 std out)
   x = yxPix(ii,2);
   y = yxPix(ii,1);
   Img(y,x) = Img(y,x)+41;
   Img(y+1,x) = Img(y+1,x)+26;
   Img(y,x+1) = Img(y,x+1)+26;
   Img(y+1,x+1) = Img(y+1,x+1)+16;
   Img(y-1,x) = Img(y-1,x)+26;
   Img(y,x-1) = Img(y,x-1)+26;
   Img(y-1,x-1) = Img(y-1,x-1)+16;
   Img(y+1,x-1) = Img(y+1,x-1)+16;
   Img(y-1,x+1) = Img(y-1,x+1)+16;   
end
Img(end,:) = [];
Img(:,end) = [];
% saturate Img at 95th prctile
cdThresh = prctile(Img(Img>0),95);

scaledPts = Img;
scaledPts(Img>cdThresh) = cdThresh;

C = imfuse(scaledPts,bigImage);
%C = imfuse(scaledPts,bigImage(:,1:imsz(2)/2));
% permuting Img to match FIJI display
% class taken from mathworks submission by Jonas Reber
roiwindow = CROIEditor(C);

% wait for roi to be assigned 
waitfor(roiwindow,'roi'); 
if ~isvalid(roiwindow) 
disp('you closed the window without applying a ROI, exiting...'); 
return 
end 

% get ROI information, like binary mask, labelled ROI, and number of 
% ROIs defined 
[mask, labels, n] = roiwindow.getROIData; 
delete(roiwindow); 

% generate cell array of images and coordinates given the cropping
croppedImages = cell(n,1);
croppedCoordinates = cell(n,1);

%% crop images and coordinates
for ii = 1:n
    [row, col] = find(labels(:,:,1)==ii);
    % get standard pixel coordinates
    leftBound = max(floor(min(col)/upscale),1);
    rightBound = min(ceil(max(col)/upscale),osz(2));
    topBound = max(floor(min(row)/upscale),1);
    bottomBound = min(ceil(max(row)/upscale),osz(1));
    % create cropped images for smaller .dax files
    croppedImages{ii} = imageData(topBound:bottomBound,leftBound:rightBound);
    % create cropped coordinates to associate .bin to above .dax files
    valids = yx(:,1) >= topBound & yx(:,1) <= bottomBound ...
        & yx(:,2) >= leftBound & yx(:,2) <= rightBound;
    croppedCoordinates{ii} = q(valids,:);
    croppedCoordinates{ii}(:,3) = croppedCoordinates{ii}(:,3)/pxtonm-leftBound+1;
    croppedCoordinates{ii}(:,4) = croppedCoordinates{ii}(:,4)/pxtonm-topBound+1; 
end

%% save out mask
maskSaveName = [tifPath(1:end-4) 'maskData.mat'];
save(maskSaveName,'mask','labels','n');

%% write out cropped images to .dax, write out cropped coordinates to .bin
for ii = 1:n
    % make the dax
    dax = DAX(fullfile(tifLoc,[tifName(1:end-ext_len) '_' num2str(ii) '.dax']));
    dax.forceFileOverwrite(1);
    dax.width = 0;
    % write into the DAX
    dax.write( croppedImages{ii} )
    dax.createInf() % creates the .inf file
    % make the bin
    i3 = Insight3();
    i3.setData(croppedCoordinates{ii}(:,[3,4]));
    i3.forceFileOverwrite(true);
    saveloc = fullfile(tifLoc,[tifName(1:end-ext_len) '_' num2str(ii) '.bin']);
    i3.write(saveloc)
end
