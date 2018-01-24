close all; clear vars; clc
%Face tracking and iPPG, using Anton M. Unakafov "best" techniques
loop=1;

%%Variable and Object set up

% Create the face detector object.

faceDetector = vision.CascadeObjectDetector;

% Create the point tracker object.
pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

% Navigate to the file
[FileName, PathName] = uigetfile({'.mp4'},'C:\Users\Jon\OneDrive - Loughborough University\final Year Project\VideoExamples');
        
% Read the file         
file_name=char(fullfile(PathName, FileName));
mov = VideoReader(file_name);

% Get file properties
vidHeight = mov.Height;
vidWidth = mov.Width;
vidLength = mov.Duration;
frame_rate = mov.FrameRate;
numFrames=floor(vidLength*frame_rate);

% Capture one frame to get its size.
videoFrame =readFrame(mov);
frameSize = size(videoFrame);

% Create the video player object. 
videoPlayer = vision.VideoPlayer('Position', [50 50 [frameSize(2), frameSize(1)]+30]);
videoPlayer2 = vision.VideoPlayer('Position', [10 10 [frameSize(2), frameSize(1)]+30]);
sectionCount=0;
frameCount = 0;
 %5 seconds of frames
sectionSample=floor(10*frame_rate);
MeanRM=zeros(sectionSample,1);
MeanGM=zeros(sectionSample,1);
MeanBM=zeros(sectionSample,1);
numPts = 0;
while loop==1
runLoop = true;




frameDif=0;


%%face Tracking
mov.CurrentTime=(sectionCount*sectionSample)/frame_rate;

while runLoop && frameDif < sectionSample

    mov.CurrentTime=(1/frame_rate)*frameCount;

    % Get the next frame.
    videoFrame = readFrame(mov);
    isempty(videoFrame);

    %videoFrame = imresize(readFrame(mov),0.5);

    

    videoFrameGray = rgb2gray(videoFrame);
    %videoFrameGray = videoFrame;
    frameCount = frameCount + 1;
    frameCount;
    sectionCount;
    sectionSample;
    frameDif=frameCount-((sectionCount*sectionSample));

     if numPts < 10
        % Detection mode.
        bbox = faceDetector.step(videoFrameGray);
        size(bbox)
        bbox(1,:)
        [m,n] = size(bbox);
        
        cropFrame = imcrop(videoFrame,bbox(1,:));
        

        if ~isempty(bbox)
            % Find corner points inside the detected region.
            facepoints = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :), 'MinQuality', 0.01);
            

            % Re-initialize the point tracker.
            xyPoints = facepoints.Location;
            numPts = size(xyPoints,1);
            release(pointTracker);
            initialize(pointTracker, xyPoints, videoFrameGray);

            % Save a copy of the points.
            oldPoints = xyPoints;

            % Convert the rectangle represented as [x, y, w, h] into an
            % M-by-2 matrix of [x,y] coordinates of the four corners. This
            % is needed to be able to transform the bounding box to display
            % the orientation of the face.
            bboxPoints = bbox2points(bbox(1, :));
           

            % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4]
            % format required by insertShape.
            bboxPolygon = reshape(bboxPoints', 1, []);
            
            cropFrame = imcrop(videoFrame,bbox(1,:));

            % Display a bounding box around the detected face.
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            
            % Display detected corners.
            videoFrame = insertMarker(videoFrame, xyPoints, '+', 'Color', 'white');
        
        else
         loop=0;
         fprintf('empty');
        end

     else
        
        % Tracking mode.
        [xyPoints, isFound] = step(pointTracker, videoFrameGray);
        visiblePoints = xyPoints(isFound, :);
        oldInliers = oldPoints(isFound, :);

        numPts = size(visiblePoints, 1);
        isempty(bbox);

        
            if ~isempty(bbox)
            % Estimate the geometric transformation between the old points
            % and the new points.
            [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);

            % Apply the transformation to the bounding box.
            bboxPoints = transformPointsForward(xform, bboxPoints);
            bbox2(1)=bboxPoints(1,1);
            bbox2(2)=bboxPoints(1,2);
            bbox2(3)=abs(bboxPoints(1,1)-bboxPoints(2,1));
            bbox2(4)=abs(bboxPoints(1,2)-bboxPoints(3,2));
            cropFrame = imcrop(videoFrame,bbox2);

            % Convert the box corners into the [x1 y1 x2 y2 x3 y3 x4 y4]
            % format required by insertShape.
            bboxPolygon = reshape(bboxPoints', 1, []);
            
            
            % Display a bounding box around the face being tracked.
           
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);
            
            % Display tracked points.
            videoFrame = insertMarker(videoFrame, visiblePoints, '+', 'Color', 'white');

            % Reset the points.
            oldPoints = visiblePoints;
            setPoints(pointTracker, oldPoints);
        
        end

     end

    % call createMask function to get the mask and the filtered image
    [BW,TEST,maskedRGBImage,OUT] = createMask(cropFrame);
    
    MeanR=mean2(maskedRGBImage(:,:,1));
    MeanRM(frameDif)=MeanR;
    MeanG=mean2(maskedRGBImage(:,:,2));
    MeanGM(frameDif)=MeanG;
    MeanB=mean2(maskedRGBImage(:,:,3));
    MeanBM(frameDif)=MeanB;
    
    %POS method
    x1t(frameDif)=MeanG-MeanB;
    x2t(frameDif)=MeanG+MeanB-(2*MeanR);
%     
%         % call createMask function to get the mask and the filtered image
%     [BW,TEST,maskedRGBImage,OUT] = createMask(cropFrame);
%     MeanR=mean2(maskedRGBImage(:,:,1));
%     MeanRM(frameCount)=MeanR;
%     MeanG=mean2(maskedRGBImage(:,:,2));
%     MeanGM(frameCount)=MeanG;
%     MeanB=mean2(maskedRGBImage(:,:,3));
%     MeanBM(frameCount)=MeanB;
%     
%     %POS method
%     x1t(frameCount)=MeanG-MeanB;
%     x2t(frameCount)=MeanG+MeanB-(2*MeanR);
    
    
    
    % plot the original image, mask and filtered image all in one figure
    figure(1);
    cla(subplot(2,3,4));
    subplot(2,3,1);imshow(cropFrame);title('Original Image');
    subplot(2,3,2);imshow(BW);title('Mask');
    subplot(2,3,3);imshow(maskedRGBImage);title('Filtered Image');
    
    subplot(2,3,4);str=num2str(frameCount);text(0.5,0.5,str);title('frame count');
    
    subplot(2,3,5);imshow(OUT);title('Outliers');
    subplot(2,3,6);imshow(TEST);title('Combined');
    
    % Display the annotated video frame using the video player object.
       step(videoPlayer, videoFrame);
%     step(videoPlayer2, cropFrame);
%     release(videoPlayer2);

    % Check whether the video player window has been closed.
    %runLoop = isOpen(videoPlayer);
 
end

% Clean up.
clear cam;

figure(2*(sectionCount+1));
subplot(4,1,1);
plot(MeanRM);
title('Raw');
% Y=fft(MeanRM);
% L=length(MeanRM);
% P2=abs(Y/L);
% P1=P2(1:L/2+1);
% P1(2:end-1)=2*P1(2:end-1);
% f=(1/30)*(0:(L/2))/L;
% figure
% plot(f,P1);
delta1=movstd(x1t,9);
delta2=movstd(x2t,9);
    
iPPG=x1t+((delta1/delta2)*x2t);
subplot(4,1,2);
plot (iPPG);
title('iPPG');

%% Filter raw signals
fc_lp = 4.0; % high cut-off
fc_hp = 0.7; % low cut-off
fs = frame_rate;

Wn = [fc_hp/(fs/2) fc_lp/(fs/2)]; % normalise with respect to Nyquist frequency

[b,a] = fir1(255, Wn, 'bandpass'); 

iPPG_filt = filter(b,a,iPPG(:));
subplot(4,1,3);
plot (iPPG_filt);
title('Filtered iPPG');
% hr=cwt(iPPG_filt, 'amor');
% figure(3);
% plot(abs(hr));
% title('Wavlet Transform');

y=fft(iPPG_filt);
f=(0:length(y)-1)*fs/length(y);
figure(2*(sectionCount+1));
subplot(4,1,4);
plot(f,abs(y));
title('FFT HR');
[~,position]=max(abs(y));
peak_f=f(position);
HR=round(peak_f*60)
sectionCount=sectionCount+1
frameDif=1;
if frameCount>=numFrames-1
    loop=0;
end
end
release(videoPlayer);
release(videoPlayer2);
release(pointTracker);
release(faceDetector);
function [BW,TEST,maskedRGBImage,OUT] = createMask(RGB) 
% Convert RGB image to HSV image
I = rgb2hsv(RGB);
% Define thresholds for 'Hue'. Modify these values to filter out different range of colors.
channel1Min = 1/360;
channel1Max = 360/360;
% Define thresholds for 'Saturation'
channel2Min = 1/255;
channel2Max = 132/255;
% Define thresholds for 'Value'
channel3Min = 88/255;
channel3Max = 255/255;
% Create mask based on chosen histogram thresholds
BW = ( (I(:,:,1) >= channel1Min) & (I(:,:,1) <= channel1Max) ) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);

M=uint8(mean2(RGB));
Sd=uint8(std2(RGB));

%OUT = ( abs((M-RGB(:,:,:) <= M-Sd) ));
% OUT = ( (abs((M-RGB(:,:,1) >= M-Sd)) & abs((M-RGB(:,:,1) <= M+Sd)))  & ...
%     (abs((M-RGB(:,:,2) >= M-Sd)) & abs((M-RGB(:,:,2) <= M+Sd))) & ...
%     (abs((M-RGB(:,:,3) >= M-Sd)) & abs((M-RGB(:,:,3) <= M+Sd))));
Gamma=0.2;
OUT1=lt(abs(M-RGB(:,:,:)),(Gamma*Sd));
size(OUT1);
OUT=((OUT1(:,:,1)) | (OUT1(:,:,2)) | (OUT1(:,:,3)));
size(BW);
size(OUT);
TEST=bitand(BW , OUT);
% Initialize output masked image based on input image.
maskedRGBImage = RGB;
% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~TEST,[1 1 3])) = 0;
end

