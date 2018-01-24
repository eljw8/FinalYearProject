close all; clear vars; clc


v = VideoWriter('newfile2.avi');
[FileName, PathName] = uigetfile();

% Read the file         
file_name=char(fullfile(PathName, FileName));
mov = VideoReader(file_name);

% Get file properties
vidHeight = mov.Height;
vidWidth = mov.Width;
vidLength = mov.Duration;
frame_rate = mov.FrameRate;
numFrames=int16(vidLength*frame_rate);

% Capture one frame to get its size.
% videoFrame =readFrame(mov);
% frameSize = size(videoFrame);

faceDetector = vision.CascadeObjectDetector('MinSize', [vidHeight/2 vidWidth/2]);

%videoPlayer = vision.VideoPlayer();
forehead=zeros(numFrames,3);
frame=1;
dropped=0;
open(v)
while hasFrame(mov)
    video = readFrame(mov);
    bbox = faceDetector(video);
    if ~isempty(bbox())
%     IFaces = insertObjectAnnotation(video, 'rectangle', bbox, 'Face');   
%     figure(1), imshow(IFaces), title('Detected faces');
    cropFrame = imcrop(video,bbox(1,:));
    writeVideo(v,cropFrame)
    xvalue=(0.5*bbox(1,3));
    yvalue=0.25*bbox(1,4);
    forehead(frame,1)=cropFrame(xvalue, yvalue, 1); 
    forehead(frame,2)=cropFrame(xvalue, yvalue, 2); 
    forehead(frame,3)=cropFrame(xvalue, yvalue, 3); 
    else
        dropped=dropped+1;    
    end
    frame=frame+1;
end
close(v)
foreheadnew = forehead(any(forehead,2),:);
figure(1)
subplot(3,1,1);plot(foreheadnew(:,1));
subplot(3,1,2);plot(foreheadnew(:,2));
subplot(3,1,3);plot(foreheadnew(:,3));

%Perform FFT 
red_FFT = abs(fft(foreheadnew(:,1)));
green_FFT = abs(fft(foreheadnew(:,2)));
blue_FFT = abs(fft(foreheadnew(:,3)));

f_axis = (0:length(foreheadnew(:,2))-1)/(length(foreheadnew(:,2))*frame_rate);

figure('Name','FFT of filtered signals'); 

subplot(4,1,1);
plot(f_axis, red_FFT,'r'); 
title('FFT of filtered signal, red');

subplot(4,1,2);
plot(f_axis, green_FFT,'g'); 
title('FFT of filtered signal, green');

subplot(4,1,3);
plot(f_axis, blue_FFT,'b'); 
subplot(4,1,4);
plot(blue_FFT);
title('FFT of filtered signal, blue');
xlabel('Frequency, Hz')
mean1=mean(foreheadnew);
mean1=mean1/256;
hsv = rgb2hsv(mean1);