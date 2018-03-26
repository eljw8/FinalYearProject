close all; clear all; clc;

%Dual camera test

Height=720;
Width=960;
% Create the face detector object.
faceDetector = vision.CascadeObjectDetector('MinSize', [Width/4 Height/4]);

cam1= webcam(1);  %% right as user sees, left as camera sees
cam2= webcam(2);  %% left as user sees, right as camera sees
% cam1.AvailableResolutions
% cam2.AvailableResolutions
txt=sprintf('%dx%d',Width,Height);
cam1.Resolution = txt;
cam2.Resolution = txt;



% preview(cam1);
% preview(cam2);

%frameCount = 1;
frameMax = 8*15;
video1=uint8(zeros(frameMax,Height,Width,3));
%Video1=gpuArray(video1);
video2=uint8(zeros(frameMax,Height,Width,3));
%Video2=gpuArray(video2);

%test=snapshot(cam1);
fprintf('filming\n');
for frameCount = 1 : frameMax
    tic
    video1(frameCount,:,:,:) = snapshot(cam1);
    %imshow(squeeze(video1(frameCount,:,:,:)));
    video2(frameCount,:,:,:) = snapshot(cam2);
    t(frameCount)=toc;
end
clear cam1;
clear cam2;
fprintf('Filming stopped\n');



% video1=gather(Video1);
% video2=gather(Video2);
plot(t);
M=mean(t);
frameRate=1/M;
sd=std(t);


vL = VideoWriter('Video\Left.avi','Motion JPEG AVI');
vR = VideoWriter('Video\Right.avi','Motion JPEG AVI');
vL.FrameRate=frameRate;
vR.FrameRate=frameRate;
open(vL);
open(vR);
fprintf('Saving\n');
for frame = 1:frameMax
    aL=squeeze(video1(frame,:,:,:));
    aR=squeeze(video2(frame,:,:,:));
    writeVideo(vL,aL);
    writeVideo(vR,aR);
end
close(vL);
close(vR);

fprintf('Saving complete\n');

fprintf('Detecting Face\n');
for frame = 1:frameMax
  figure(2);
%   a=reshape(video1(frame,:,:,:),480,640,3);
  a=squeeze(video1(frame,:,:,:));
  grayVideo1=rgb2gray(a);
  bbox1 = faceDetector.step(grayVideo1);
  bboxPoints1 = bbox2points(bbox1(1, :));
  bboxPolygon1 = reshape(bboxPoints1', 1, []);
  a = insertShape(a, 'Polygon', bboxPolygon1, 'LineWidth', 3);
  subplot(1,2,1);
  imshow(a);
  b=squeeze(video2(frame,:,:,:));
  grayVideo2=rgb2gray(b);
  bbox2 = faceDetector.step(grayVideo2);
  bboxPoints2 = bbox2points(bbox2(1, :));
  bboxPolygon2 = reshape(bboxPoints2', 1, []);
  b = insertShape(b, 'Polygon', bboxPolygon2, 'LineWidth', 3);
  subplot(1,2,2);
  imshow(b);
  %mov(frame)=im2frame(a);
end

