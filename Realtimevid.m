close all; clear all; clc;

% Create the face detector object.
%faceDetector = vision.CascadeObjectDetector('MinSize', []);

% Create the point tracker object.
pointTracker = vision.PointTracker('MaxBidirectionalError', 2);

% Navigate to the file
[FileName, PathName] = uigetfile({'.mp4'},'C:\Users\Jon\OneDrive - Loughborough University\final Year Project\VideoExamples');
        
% Read the file         
file_name=char(fullfile(PathName, FileName));
mov = VideoReader(file_name);


% Create the webcam object.
%cam = webcam();

% Capture one frame to get its size.
vidHeight = mov.Height;
vidWidth = mov.Width;
vidLength = mov.Duration;
frame_rate = mov.FrameRate;
numFrames=int16(vidLength*frame_rate);

videoFrame = readFrame(mov);
frameSize = size(videoFrame);
frame_rate = mov.FrameRate;
frameMax=numFrames-10;
mov.CurrentTime=0;
%newVideoFrame=zeros(frameMax, frameSize(1), frameSize(2),frameSize(3));
% Create the face detector object.
faceDetector = vision.CascadeObjectDetector('MinSize', [vidHeight/3 vidWidth/3]);


% Create the video player object.
videoPlayer = vision.VideoPlayer('Position', [100 100 [frameSize(2), frameSize(1)]+30]);
writerObj = VideoWriter('myVideo2.avi');
open(writerObj);
runLoop = true;
numPts = 0;
frameCount = 0;
timerVal=tic;
tic;
while runLoop && frameCount < frameMax

    % Get the next frame.
    videoFrame = readFrame(mov);
    videoFrameGray = rgb2gray(videoFrame);
    frameCount = frameCount + 1;

    if numPts < 10
        % Detection mode.
        bbox = faceDetector.step(videoFrameGray);

        if ~isempty(bbox)
            % Find corner points inside the detected region.
            points = detectMinEigenFeatures(videoFrameGray, 'ROI', bbox(1, :));

            % Re-initialize the point tracker.
            xyPoints = points.Location;
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

            % Display a bounding box around the detected face.
            videoFrame = insertShape(videoFrame, 'Polygon', bboxPolygon, 'LineWidth', 3);

            % Display detected corners.
            videoFrame = insertMarker(videoFrame, xyPoints, '+', 'Color', 'white');
        end

    else
        % Tracking mode.
        [xyPoints, isFound] = step(pointTracker, videoFrameGray);
        visiblePoints = xyPoints(isFound, :);
        oldInliers = oldPoints(isFound, :);

        numPts = size(visiblePoints, 1);

        if numPts >= 10
            % Estimate the geometric transformation between the old points
            % and the new points.
            [xform, oldInliers, visiblePoints] = estimateGeometricTransform(...
                oldInliers, visiblePoints, 'similarity', 'MaxDistance', 4);

            % Apply the transformation to the bounding box.
            bboxPoints = transformPointsForward(xform, bboxPoints);

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
    
    min_x = min(bboxPoints(:,1));
    max_x = max(bboxPoints(:,1));
    min_y = min(bboxPoints(:,2));
    max_y = max(bboxPoints(:,2));
    cropVideoFrame=imcrop(videoFrame,[min_x min_y 500 500]);
    Red(frameCount,:,:)=(cropVideoFrame(:,:,1));
    Green(frameCount,:,:)=(cropVideoFrame(:,:,2));
    Blue(frameCount,:,:)=(cropVideoFrame(:,:,3));
%     meanRed(frameCount)=mean2(cropVideoFrame(:,:,1));
%     meanGreen(frameCount)=mean2(cropVideoFrame(:,:,2));
%     meanBlue(frameCount)=mean2(cropVideoFrame(:,:,3));
    
    writeVideo(writerObj,cropVideoFrame);
    %mov = immovie(cropVideoFrame);

    % Display the annotated video frame using the video player object.
    step(videoPlayer, videoFrame);

    % Check whether the video player window has been closed.
    runLoop = isOpen(videoPlayer);
    toc;
    test1(frameCount)=toc;
end

% Clean up.
clear cam;
close(writerObj);
release(videoPlayer);
release(pointTracker);
release(faceDetector);
% a=diff(test1);
% mean(a)
% frameRate=1/mean(a)
fs = frame_rate;
% Red(1:400,:)=[];
% Green(1:400,:)=[];
% Blue(1:400,:)=[];
[a1,a2,a3]=size(Red);
for a = 1:a1
    
[BW,TEST,maskedRGBImage,OUT] = createMask(Red(a,:,:),Green(a,:,:),Blue(a,:,:));

Masked(a,:,:,:)=maskedRGBImage;

end
sumRed=sum(sum(Red,2),3);
sumGreen=sum(sum(Green,2),3);
sumBlue=sum(sum(Blue,2),3);


x=(0:length(sumRed)-1);

stdRed=std(sumRed);
stdGreen=std(sumGreen);
stdBlue=std(sumBlue);
mnRed=mean(sumRed);
mnGreen=mean(sumGreen);
mnBlue=mean(sumBlue);
X(1,:)=(sumRed-mnRed)/stdRed;
X(2,:)=(sumBlue-mnBlue)/stdBlue;
X(3,:)=(sumGreen-mnGreen)/stdGreen;
figure(1);
plot(x,sumRed, 'r',x,sumGreen, 'g',x,sumBlue, 'b');
meanRed=mean(mean(Red,2),3);
meanGreen=mean(mean(Green,2),3);
meanBlue=mean(mean(Blue,2),3);
figure(2);
plot(x,meanRed, 'r',x,meanGreen, 'g',x,meanBlue, 'b');
figure(3);
plot(x,X(1,:), 'r',x,X(2,:), 'g',x,X(3,:), 'b');

yR=fft(X(1,:));
yG=fft(X(2,:));
yB=fft(X(3,:));
f=(0:length(yR)-1)*fs/length(yR);
figure(4);
subplot (3,1,1);
plot(f,abs(yR));
subplot (3,1,2);
plot(f,abs(yG));
subplot (3,1,3);
plot(f,abs(yB));
% subplot(3,1,1);
% plot(meanRed);
% subplot(3,1,2);
% plot(meanGreen);
% subplot(3,1,3);
% plot(meanBlue);

% yR=fft(meanRed);
% yG=fft(meanGreen);
% yB=fft(meanBlue);
% T=1/frameRate;
% f=(0:frameMax-1)/T*frameMax;
% figure(2);
% subplot (3,1,1);
% plot(f,abs(yR));
% subplot (3,1,2);
% plot(f,abs(yG));
% subplot (3,1,3);
% plot(f,abs(yB));

%% Filter raw signals
fc_lp = 4.0; % high cut-off
fc_hp = 0.7; % low cut-off
%fs = frameRate;

Wn = [fc_hp/(fs/2) fc_lp/(fs/2)]; % normalise with respect to Nyquist frequency

%[b,a] = fir1(255, Wn, 'bandpass'); 
%bpFilt = designfilt('bandpassiir','FilterOrder',100, ...
%         'CutoffFrequency1',fc_hp,'CutoffFrequency2',fc_lp, ...
%         'SampleRate',fs);
bpFilt = designfilt('bandpassiir', 'FilterOrder', 100, ...
                    'HalfPowerFrequency1', 0.7, 'HalfPowerFrequency2', ...
                    4, 'SampleRate', 25);
   
iPPG_filtR = filter(bpFilt,X(1,:));
iPPG_filtG = filter(bpFilt,X(2,:));
iPPG_filtB = filter(bpFilt,X(3,:));
figure(5);
x=(0:length(iPPG_filtR)-1);
plot(x,iPPG_filtR, 'r',x,iPPG_filtG, 'g',x,iPPG_filtB, 'b');
% subplot (3,1,1);
% plot(iPPG_filtR);
% subplot (3,1,2);
% plot(iPPG_filtG);
% subplot (3,1,3);
% plot(iPPG_filtB);

yR=fft(iPPG_filtR);
yG=fft(iPPG_filtG);
yB=fft(iPPG_filtB);
f=(0:length(yR)-1)*fs/length(yR);
figure(6);
subplot (3,1,1);
plot(f,abs(yR));
subplot (3,1,2);
plot(f,abs(yG));
subplot (3,1,3);
plot(f,abs(yB));

% yR=fft(meanRed);
% yG=fft(meanGreen);
% yB=fft(meanBlue);
% f=(0:length(yR)-1)*fs/length(yR);
% figure(5);
% subplot (3,1,1);
% plot(f,abs(yR));
% subplot (3,1,2);
% plot(f,abs(yG));
% subplot (3,1,3);
% plot(f,abs(yB));

[~,position]=max(abs(yR));
peak_f=f(position);
HR_R=(peak_f*60)
%HR_R=round(peak_f*60)
[~,position]=max(abs(yG));
peak_f=f(position);
HR_G=(peak_f*60)
[~,position]=max(abs(yB));
peak_f=f(position);
HR_B=(peak_f*60)

function [BW,TEST,maskedRGBImage,OUT] = createMask(R,G,B) 
% Convert RGB image to HSV image
RGB(:,:,1)=R;
RGB(:,:,2)=G;
RGB(:,:,3)=B;
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