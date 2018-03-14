clear all; close all; clc;

%recreating method from Qi Zhang et al. 
[FileName, PathName] = uigetfile();
file_name=char(fullfile(PathName, FileName));
v = VideoReader(file_name);
vidHeight = v.Height;
vidWidth = v.Width;

faceDetector = vision.CascadeObjectDetector('MinSize', [vidHeight/3 vidWidth/3]);

dropped=0;
test=0;
one=0;
two=0;
three=0;
grater=0;

while hasFrame(v)
    video=readFrame(v);
    %YCBCR=rgb2ycbcr(video);
    GRAY=rgb2gray(video);
%     figure(1)
%     subplot(2,4,1);
%     imshow(video);
%     subplot(2,4,2);
%     imshow(video(:,:,1));
%     title('Red')
%     subplot(2,4,3);
%     imshow(video(:,:,2));
%     title('Green')
%     subplot(2,4,4);
%     imshow(video(:,:,3));
%     title('Blue')
%     subplot(2,4,6);
%     imshow(YCBCR(:,:,1));
%     title('Y')
%     subplot(2,4,7);
%     imshow(YCBCR(:,:,2));
%     title('Cb')
%     subplot(2,4,8);
%     imshow(YCBCR(:,:,3));
%     title('Cr')
    %bbox = faceDetector(YCBCR(:,:,1));
    bbox = faceDetector(GRAY);
    IFaces = insertObjectAnnotation(video, 'rectangle', bbox, 'Face');  
    %subplot(2,4,5);
    %A=[v.CurrentTime, bbox];
    if ~isempty(bbox())
        [m,n] = size(bbox);
        if m==1
            one=one+1;
        elseif m==2
            two=two+1;
        elseif m==3
            three=three+1;
        else
            grater=grater+1;
        end
        test=test+1;
    else
       dropped=dropped+1; 
    end
    figure(1);
    imshow(IFaces);
    title('detect');
    
%     MASK=(YCBCR(:,:,2)>=77 & YCBCR(:,:,2)<=127 & YCBCR(:,:,3)>=133 & YCBCR(:,:,3)<=173);
%     subplot(2,4,5);
%     %test=bitand(video,MASK);
%     
%     video(repmat(~MASK,[1 1 3])) = 0;
%     
%     imshow(video);
%     title('mask')
%     hold on;
end
