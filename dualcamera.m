close all; clear all; clc;

%Dual camera test

cam1= webcam(1);
cam2= webcam(2);
% cam1.AvailableResolutions
% cam2.AvailableResolutions
cam1.Resolution = '960x720';
cam2.Resolution = '960x720';



% preview(cam1);
% preview(cam2);

%frameCount = 1;
frameMax = 8*15;
video1=uint8(zeros(frameMax,720,960,3));
%Video1=gpuArray(video1);
video2=uint8(zeros(frameMax,720,960,3));
%Video2=gpuArray(video2);

test=snapshot(cam1);
'filming'
for frameCount = 1 : frameMax
    tic
    video1(frameCount,:,:,:) = snapshot(cam1);
    %imshow(squeeze(video1(frameCount,:,:,:)));
    video2(frameCount,:,:,:) = snapshot(cam2);
    t(frameCount)=toc;
end
'done'
% video1=gather(Video1);
% video2=gather(Video2);
plot(t);
mean(t)
1/mean(t)
std(t)

for frame = 1:frameMax
  figure(2);
%   a=reshape(video1(frame,:,:,:),480,640,3);
  a=squeeze(video1(frame,:,:,:));
  grayVideo1=rgb2gray(a);
  subplot(1,2,1);
  imshow(a);
  b=squeeze(video2(frame,:,:,:));
  grayVideo2=rgb2gray(b);
  subplot(1,2,2);
  imshow(b);
  %mov(frame)=im2frame(a);
end

