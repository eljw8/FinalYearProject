function []=Save(Left,Right,frameRate,vL,vR)
%%function to save two camera streams Left and Right
%Left and Right must be 4d arrays, (frame number, height, width, 3(colour))
%frameRate= frame rate of video to be written
%num = suffix of file names being saved
%vL,vR =VideoWriter objects


vL.FrameRate=frameRate;

vR.FrameRate=frameRate;
open(vL);
open(vR);
[a,~,~,~] = size(Left);
fprintf('Saving\n');
for frame = 1:a
    aL=squeeze(Left(frame,:,:,:));
    aR=squeeze(Right(frame,:,:,:));
    writeVideo(vL,aL);
    writeVideo(vR,aR);
end
close(vL);
close(vR);

fprintf('Saving complete\n');
end
