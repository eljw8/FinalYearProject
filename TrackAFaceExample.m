clear all; close all; clc;
load Functions\matlab.mat
%% Track a Face in Scene
%
%% 
% Create System objects for reading and displaying video and for drawing a bounding box of the object. 
videoFileReaderL = vision.VideoFileReader('Video\Left4.avi');
videoFileReaderR = vision.VideoFileReader('Video\Right4.avi');
videoPlayerL = vision.VideoPlayer();
videoPlayerR = vision.VideoPlayer();
faceDetector = vision.CascadeObjectDetector('MinSize', [720/4 960/4]);

%% 
% Read the first video frame, which contains the object, define the region.
objectFrameL = step(videoFileReaderL);
bboxL=faceDetector(objectFrameL);
objectRegionL = bboxL;
objectFrameR = step(videoFileReaderR);
bboxR=faceDetector(objectFrameR);
objectRegionR = bboxR;
%% 
% As an alternative, you can use the following commands to select the object region using a mouse. The object must occupy the majority of the region. 
% figure; imshow(objectFrame);
% objectRegion=round(getPosition(imrect))
%% 
% Show initial frame with a red bounding box.
objectImageL = insertShape(objectFrameL,'Rectangle',objectRegionL,'Color','blue'); 
objectImageR = insertShape(objectFrameR,'Rectangle',objectRegionR,'Color','red'); 
figure;
imshowpair(objectImageL, objectImageR, 'montage');
title('box shows object region');
%% 
% Detect interest points in the object region.
pointsL = detectMinEigenFeatures(rgb2gray(objectFrameL),'ROI',objectRegionL);
pointsR = detectMinEigenFeatures(rgb2gray(objectFrameR),'ROI',objectRegionR);
%% 
% Display the detected points.
pointImageL = insertMarker(objectFrameL,pointsL.Location,'+','Color','white');
pointImageR = insertMarker(objectFrameR,pointsR.Location,'+','Color','white');
figure;
imshowpair(pointImageL, pointImageR, 'montage');
title('Detected interest points');

%% 
% Matching features
[featuresL,valid_pointsL] = extractFeatures(rgb2gray(objectFrameL),pointsL);
[featuresR,valid_pointsR] = extractFeatures(rgb2gray(objectFrameR),pointsR);
indexPairs = matchFeatures(featuresL,featuresR);
matchedPointsL = valid_pointsL(indexPairs(:,1),:);
matchedPointsR = valid_pointsR(indexPairs(:,2),:);
undistortedPointsL = undistortPoints(matchedPointsL.Location,stereoParams.CameraParameters1);
undistortedPointsR = undistortPoints(matchedPointsR.Location,stereoParams.CameraParameters2);
figure; showMatchedFeatures(objectImageL,objectImageR,matchedPointsL,matchedPointsR,'montage');
figure; showMatchedFeatures(objectImageL,objectImageR,undistortedPointsL,undistortedPointsR,'montage');
%% 
% Create a tracker object.
tracker = vision.PointTracker('MaxBidirectionalError',1);
%% 
% Initialize the tracker.
initialize(tracker,matchedPointsL.Location,objectFrameL);
initialize(tracker,matchedPointsR.Location,objectFrameR);
%% 
% Read, track, display points, and results in each video frame.
while ~isDone(videoFileReaderL)
      frameL = step(videoFileReaderL);
      frameR = step(videoFileReaderR);
      [pointsL, validityL] = step(tracker,frameL);
      [pointsR, validityR] = step(tracker,frameR);
      outL = insertMarker(frameL,pointsL(validityL, :),'+','color','blue');
      outR = insertMarker(frameR,pointsR(validityR, :),'+','color','red');
      step(videoPlayerL,outL);
      step(videoPlayerR,outR);
end
%% 
% Release the video reader and player.
release(videoPlayerL);
release(videoPlayerR);
release(videoFileReaderL);
release(videoFileReaderR);