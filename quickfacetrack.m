close all; clear vars; clc;


%v = VideoWriter('newfile2.avi');
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

faceDetector = vision.CascadeObjectDetector('MinSize', [vidHeight/3 vidWidth/3]);

%videoPlayer = vision.VideoPlayer();
forehead=zeros(numFrames,3);
frame=1;
dropped=0;
%open(v)
while hasFrame(mov)
    video = readFrame(mov);
    grayvideo=rgb2gray(readFrame(mov));
    bbox = faceDetector(grayvideo);
    if ~isempty(bbox())
%     IFaces = insertObjectAnnotation(video, 'rectangle', bbox, 'Face');   
%     figure(1), imshow(IFaces), title('Detected faces');
    cropFrame = imcrop(video,bbox(1,:));
    %%writeVideo(v,cropFrame)
    xvalue=abs(0.5*bbox(1,3));
    yvalue=abs(0.25*bbox(1,4));
    forehead(frame,1)=cropFrame(xvalue, yvalue, 1); 
    forehead(frame,2)=cropFrame(xvalue, yvalue, 2); 
    forehead(frame,3)=cropFrame(xvalue, yvalue, 3); 
    RGB=[cropFrame(xvalue, yvalue, 1),cropFrame(xvalue, yvalue, 2),cropFrame(xvalue, yvalue, 3)];
    RGB2=double(RGB)/256;
    HSV=rgbtohsv(RGB2(1),RGB2(2),RGB2(3));
    
    
    [BW,TEST,maskedRGBImage,OUT] = createMask(HSV(2),HSV(3),cropFrame);
    figure(10);
    cla(subplot(2,3,4));
    subplot(2,3,1);imshow(cropFrame);title('Original Image');
    subplot(2,3,2);imshow(BW);title('Mask');
    subplot(2,3,3);imshow(maskedRGBImage);title('Filtered Image');
    
    %%subplot(2,3,4);str=num2str(frameCount);text(0.5,0.5,str);title('frame count');
    
    subplot(2,3,5);imshow(OUT);title('Outliers');
    subplot(2,3,6);imshow(TEST);title('Combined');
    else
        cropFrame = video;
    %%writeVideo(v,cropFrame)
    xvalue=(0.5*bbox(1,3));
    yvalue=0.25*bbox(1,4);
    forehead(frame,1)=cropFrame(xvalue, yvalue, 1); 
    forehead(frame,2)=cropFrame(xvalue, yvalue, 2); 
    forehead(frame,3)=cropFrame(xvalue, yvalue, 3); 
    RGB=[cropFrame(xvalue, yvalue, 1),cropFrame(xvalue, yvalue, 2),cropFrame(xvalue, yvalue, 3)];
    RGB2=double(RGB)/256;
    HSV=rgbtohsv(RGB2(1),RGB2(2),RGB2(3));
    
    
    [BW,TEST,maskedRGBImage,OUT] = createMask(HSV(2),HSV(3),cropFrame);
    figure(10);
    cla(subplot(2,3,4));
    subplot(2,3,1);imshow(cropFrame);title('Original Image');
    subplot(2,3,2);imshow(BW);title('Mask');
    subplot(2,3,3);imshow(maskedRGBImage);title('Filtered Image');
    
    %%subplot(2,3,4);str=num2str(frameCount);text(0.5,0.5,str);title('frame count');
    
    subplot(2,3,5);imshow(OUT);title('Outliers');
    subplot(2,3,6);imshow(TEST);title('Combined');
        dropped=dropped+1;    
    end
    frame=frame+1;
end
%close(v)
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
mean1=mean(foreheadnew)
mean1=mean1/256;
mean2=mean(nonzeros(foreheadnew))
mean3=gather(mean1);
hsv = rgb2hsv(mean3);

%% Filter raw signals
fc_lp = 4.0; % high cut-off
fc_hp = 0.7; % low cut-off
fs = frame_rate;

Wn = [fc_hp/(fs/2) fc_lp/(fs/2)]; % normalise with respect to Nyquist frequency

[b,a] = fir1(255, Wn, 'bandpass'); 

iPPG_filtR = filter(b,a,foreheadnew(:,1));
iPPG_filtG = filter(b,a,foreheadnew(:,2));
iPPG_filtB = filter(b,a,foreheadnew(:,3));
figure(3);
subplot (3,1,1);
plot(iPPG_filtR);
subplot (3,1,2);
plot(iPPG_filtG);
subplot (3,1,3);
plot(iPPG_filtB);

yR=fft(iPPG_filtR);
yG=fft(iPPG_filtG);
yB=fft(iPPG_filtB);
f=(0:length(yR)-1)*fs/length(yR);
figure(4);
subplot (3,1,1);
plot(f,abs(yR));
subplot (3,1,2);
plot(f,abs(yG));
subplot (3,1,3);
plot(f,abs(yB));
% [~,position]=max(abs(yR));
% peak_f=f(position);
% HR_R=round(peak_f*60)
% [~,position]=max(abs(yG));
% peak_f=f(position);
% HR_G=round(peak_f*60)
% [~,position]=max(abs(yB));
% peak_f=f(position);
% HR_B=round(peak_f*60)


function [BW,TEST,maskedRGBImage,OUT] = createMask(sat,val,RGB) 
%sat=gather(sat);
%val=gather(val);
%RGB=gather(RGB);
% Convert RGB image to HSV image
I = rgbtohsv(RGB(:,:,1),RGB(:,:,2),RGB(:,:,3));
%I = rgbtohsv(RGB);
% Define thresholds for 'Hue'. Modify these values to filter out different range of colors.
channel1Min = 0;
channel1Max = 1;
% Define thresholds for 'Saturation'
satmin=sat-0.025;
if satmin<0; satmin=0;end;
satmax=sat+0.025;
if satmax>1; satmax=1;end;

% Define thresholds for 'Value'
valmin=val-0.035;
if valmin<0;valmin=0;end;
valmax=val+0.035;
if valmax>1;valmax=1;end;

% Create mask based on chosen histogram thresholds
% BW = ( (I(:,:,1) >= channel1Min) & (I(:,:,1) <= channel1Max) ) & ...
%     (I(:,:,2) >= satmin ) & (I(:,:,2) <= satmax) & ...
%     (I(:,:,3) >= valmin ) & (I(:,:,3) <= valmax);
BW = ( (I(1) >= channel1Min) & (I(1) <= channel1Max) ) & ...
    (I(2) >= satmin ) & (I(2) <= satmax) & ...
    (I(3) >= valmin ) & (I(3) <= valmax);

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

function [I]=rgbtohsv(r,g,b)
A=[r g b];

var_min=min(A);
var_max=max(A);
del_max=var_max - var_min;

v = var_max;

if del_max == 0
    h=0;
    s=0;
else
    s=del_max/var_max;
    del_R = ( ( ( var_max - r ) / 6 ) + ( del_max / 2 ) ) / del_max;
    del_G = ( ( ( var_max - g ) / 6 ) + ( del_max / 2 ) ) / del_max;
    del_B = ( ( ( var_max - b ) / 6 ) + ( del_max / 2 ) ) / del_max;
    if r ==var_max
        h=del_B - del_G;
    elseif g == var_max
        h=(1/3)+del_R-del_B;
    elseif b == var_max
        h=(2/3)+del_G-del_R;         
    end
  if h<0
      h=h+1;
  end
  if h>1
      h=h-1;
  end
            
end
I=[h,s,v];
end

