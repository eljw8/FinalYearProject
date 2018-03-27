close all; clc;
%clear all;

%Dual camera test

Height=720;
Width=960;
Resolution=[Height,Width];

% Create the face detector object.
faceDetector = vision.CascadeObjectDetector('MinSize', [Width/4 Height/4]);


frameMax = 8*15;

%%record video with webcam


% [Left,Right,M,t]=record(frameMax,Resolution);
% 
% plot(t);
% 
% mkdir('Video');
% num=4;
% L=sprintf('Video\\Left%d.avi',num);
% R=sprintf('Video\\Right%d.avi',num);
% vL = VideoWriter(L,'Motion JPEG AVI');
% vR = VideoWriter(R,'Motion JPEG AVI');
% 
% Save(Left,Right,M,vL,vR)

% vL = VideoReader('Video\Left3.avi');
% vR = VideoReader('Video\Right3.avi');
% [Left,Right] = videoImport(vL,vR);

% Detect feature points
imagePoints1 = detectMinEigenFeatures(rgb2gray(squeeze(Left(1,:,:,:))), 'MinQuality', 0.0001);

% Visualize detected points
figure
imshow(squeeze(Left(1,:,:,:)), 'InitialMagnification', 50);
title('150 Strongest Corners from the First Image');
hold on
plot(selectStrongest(imagePoints1, 150));

% Create the point tracker
tracker = vision.PointTracker('MaxBidirectionalError', 1, 'NumPyramidLevels', 5);

% Initialize the point tracker
imagePoints1 = imagePoints1.Location;
initialize(tracker, imagePoints1, squeeze(Left(1,:,:,:)));

% % Track the points
% [imagePoints2, validIdx] = step(tracker, squeeze(Right(1,:,:,:)));
% matchedPoints1 = imagePoints1(validIdx, :);
% matchedPoints2 = imagePoints2(validIdx, :);

% % Visualize correspondences
% figure
% showMatchedFeatures(squeeze(Left(1,:,:,:)), squeeze(Right(1,:,:,:)), matchedPoints1, matchedPoints2);
% title('Tracked Features');
frame=1;
    I1 = undistortImage(squeeze(Left(frame,:,:,:)), stereoParams.CameraParameters1);
    I2 = undistortImage(squeeze(Right(frame,:,:,:)), stereoParams.CameraParameters2);
    
    grayI1=rgb2gray(I1);
    grayI2=rgb2gray(I2);
    bbox1=faceDetector(grayI1);
    bbox2=faceDetector(grayI2);
    points1=detectMinEigenFeatures(grayI1,'ROI',bbox1);
    points2=detectMinEigenFeatures(grayI2,'ROI',bbox2);
    
    pointImage1 = insertMarker(I1,points1.Location,'+','Color','white');
    pointImage2 = insertMarker(I2,points2.Location,'+','Color','white');
    figure(1);
    imshowpair(I1, I2, 'montage')
    title('Detected interest points');
    
    tracker = vision.PointTracker('MaxBidirectionalError',1);
    initialize(tracker,points1.Location,objectFrame);
    
for frame = 2:frameMax
    
    I1 = undistortImage(squeeze(Left(frame,:,:,:)), stereoParams.CameraParameters1);
    I2 = undistortImage(squeeze(Right(frame,:,:,:)), stereoParams.CameraParameters2);
    
    grayI1=rgb2gray(I1);
    grayI2=rgb2gray(I2);
    bbox1=faceDetector(grayI1);
    bbox2=faceDetector(grayI2);
    if ~isempty(bbox1) && ~isempty(bbox2)
        cropI1=imcrop(I1,bbox1);
        cropI2=imcrop(I2,bbox2);
        [imagePoints1, validIdx] = tracker(I1);
        points1 = detectMinEigenFeatures(grayI1, 'ROI', bbox1,'MinQuality', 0.001);
        points2 = detectMinEigenFeatures(grayI2, 'ROI', bbox2,'MinQuality', 0.001);
        I1=insertMarker(I1,selectStrongest(points1,50), 'color', 'red');
        I1=insertMarker(I1,imagePoints1, 'color', 'white');
        I2=insertMarker(I2,selectStrongest(points2,50));
        figure (2);
        imshowpair(I1, I2, 'montage');
        hold on;
    end
end


% fprintf('Detecting Face\n');
% for frame = 1:frameMax
%   figure(2);
% %   a=reshape(video1(frame,:,:,:),480,640,3);
%   a=squeeze(Left(frame,:,:,:));
%   grayVideo1=rgb2gray(a);
%   bbox1 = faceDetector.step(grayVideo1);
%   if ~isempty(bbox1)
%   bboxPoints1 = bbox2points(bbox1(1, :));
%   % Detect feature points in the face region.
%   points1 = detectMinEigenFeatures(grayVideo1, 'ROI', bbox1);
%   bboxPolygon1 = reshape(bboxPoints1', 1, []);
% %   cropLeft(frame,:,:,:)=imcrop(a,bbox1(1, :));
%   a = insertShape(a, 'Polygon', bboxPolygon1, 'LineWidth', 3);
%   subplot(1,2,1);
%   imshow(a);
%   hold on
%   plot(points1.selectStrongest(10));
%   end
%   b=squeeze(Right(frame,:,:,:));
%   grayVideo2=rgb2gray(b);
%   bbox2 = faceDetector.step(grayVideo2);
%   if ~isempty(bbox2)
%   bboxPoints2 = bbox2points(bbox2(1, :));
%   bboxPolygon2 = reshape(bboxPoints2', 1, []);
%   b = insertShape(b, 'Polygon', bboxPolygon2, 'LineWidth', 3);
%   subplot(1,2,2);
%   imshow(b);
%   end
%   %mov(frame)=im2frame(a);
% end
% fprintf('Detection complete\n')

% function [Left,Right,M,t] = record(frameMax,resolution)
% %%record dual camera footage using webcams
% %frameMax = total number of frames to be recorded
% %resolution = [height,width] array of resolution that is being recorded
% 
%     cam1= webcam(1);  %% right as user sees, left as camera sees
%     cam2= webcam(2);  %% left as user sees, right as camera sees
%     txt=sprintf('%dx%d',resolution(2),resolution(1));
%     cam1.Resolution = txt;
%     cam2.Resolution = txt;
%     video1=uint8(zeros(frameMax,resolution(1),resolution(2),3));
%     video2=uint8(zeros(frameMax,resolution(1),resolution(2),3));
%     fprintf('filming\n');
%         for frameCount = 1 : frameMax
%             tic
%             video1(frameCount,:,:,:) = snapshot(cam1);
%             %imshow(squeeze(video1(frameCount,:,:,:)));
%             video2(frameCount,:,:,:) = snapshot(cam2);
%             t(frameCount)=toc;
%         end
%     clear cam1;
%     clear cam2;
%     fprintf('Filming stopped\n');
%     t(1)=[];
%     Mean=mean(t);
%     M=1/Mean;
%     Left=video1;
%     Right=video2;
% 
% end

% function []=Save(Left,Right,frameRate,vL,vR)
% %%function to save two camera streams Left and Right
% %Left and Right must be 4d arrays, (frame number, height, width, 3(colour))
% %frameRate= frame rate of video to be written
% %num = suffix of file names being saved
% %vL,vR =VideoWriter objects
% 
% 
% vL.FrameRate=frameRate;
% 
% vR.FrameRate=frameRate;
% open(vL);
% open(vR);
% [a,~,~,~] = size(Left);
% fprintf('Saving\n');
% for frame = 1:a
%     aL=squeeze(Left(frame,:,:,:));
%     aR=squeeze(Right(frame,:,:,:));
%     writeVideo(vL,aL);
%     writeVideo(vR,aR);
% end
% close(vL);
% close(vR);
% 
% fprintf('Saving complete\n');
% end

% function[Left, Right] = videoImport(vL,vR)
% %%imports two videos, one as left one as right
% %vL = video object cameras left
% %vR = video object cameras Right
% fprintf('Importing\n');
%     totalFrameL=vL.NumberOfFrames;
%     totalFrameR=vR.NumberOfFrames;
%     heightL=vL.Height;
%     heightR=vR.Height;
%     widthL=vL.Width;
%     widthR=vR.Width;
%     
%     Left=uint8(zeros(totalFrameL,heightL,widthL,3));
%     Right=uint8(zeros(totalFrameR,heightR,widthR,3));
% 
%     
%     for frameNum = 1:totalFrameL
%         Left(frameNum,:,:,:)=read(vL,frameNum);
%     end
%     for frameNum = 1:totalFrameR
%         Right(frameNum,:,:,:)=read(vR,frameNum);
%     end
% fprintf('Importing complete\n');
% end
