%taking photos, to use for calibration function
%this is opened by calling stereocameraCalibrator function

close all; clear all; clc;

cam1= webcam(1);
cam2= webcam(2);

frameMax=20;

for frameCount = 1 : frameMax
    tic
    a =snapshot(cam1);
    b =snapshot(cam2);
    s0='Pictures\';
    sA='A\';
    s1='a';
    s2=num2str(frameCount);
    s3='.png';
    s=strcat(s0,sA,s1,s2,s3);
    imwrite(a,s); 
    
    sB='B\';
    s4='b';
    s5=strcat(s0,sB,s4,s2,s3);
    imwrite(b,s5);
    t(frameCount)=toc;
    disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
    pause;
end

clear cam1;
clear cam2;

stereoCameraCalibrator(strcat(s0,sA),strcat(s0,sB),20)


c=imread('Pictures\A\a10.png');
d=imread('Pictures\B\b10.png');
[J1,J2] = rectifyStereoImages(c,d,stereoParams,'OutputView','valid');
figure(2);
subplot(2,2,1);
imshow(c);
subplot(2,2,2); 
imshow(d);
subplot(2,2,3);
imshow(J1);
subplot(2,2,4);
imshow(J2);
figure(3);
imshow(stereoAnaglyph(J1,J2));