function [Left,Right,M,t] = record(frameMax,resolution)
%record dual camera footage using webcams
%frameMax = total number of frames to be recorded
%resolution = [height,width] array of resolution that is being recorded

    cam1= webcam(1);  %% right as user sees, left as camera sees
    cam2= webcam(2);  %% left as user sees, right as camera sees
    txt=sprintf('%dx%d',resolution(2),resolution(1));
    cam1.Resolution = txt;
    cam2.Resolution = txt;
    video1=uint8(zeros(frameMax,resolution(1),resolution(2),3));
    video2=uint8(zeros(frameMax,resolution(1),resolution(2),3));
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
    t(1)=[];
    Mean=mean(t);
    M=1/Mean;
    Left=video1;
    Right=video2;

end