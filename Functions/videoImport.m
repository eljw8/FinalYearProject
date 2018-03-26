function[Left, Right] = videoImport(vL,vR)
%%imports two videos, one as left one as right
%vL = video object cameras left
%vR = video object cameras Right
fprintf('Importing\n');
    totalFrameL=vL.NumberOfFrames;
    totalFrameR=vR.NumberOfFrames;
    heightL=vL.Height;
    heightR=vR.Height;
    widthL=vL.Width;
    widthR=vR.Width;
    
    Left=uint8(zeros(totalFrameL,heightL,widthL,3));
    Right=uint8(zeros(totalFrameR,heightR,widthR,3));

    
    for frameNum = 1:totalFrameL
        Left(frameNum,:,:,:)=read(vL,frameNum);
    end
    for frameNum = 1:totalFrameR
        Right(frameNum,:,:,:)=read(vR,frameNum);
    end
fprintf('Importing complete\n');
end