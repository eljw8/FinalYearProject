%taking photos

close all; clear all; clc;

cam1= webcam(1);
cam2= webcam(2);

frameMax=20;

for frameCount = 1 : frameMax
    tic
    a =snapshot(cam1);
    b =snapshot(cam2);
    s1='a';
    s2=num2str(frameCount);
    s3='.jpeg';
    s=strcat(s1,s2,s3);
    imwrite(a,s);
    
    
    s4='b';
    s5=strcat(s4,s2,s3);
    imwrite(b,s5);
    t(frameCount)=toc;
end