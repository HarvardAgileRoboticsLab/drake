clear; clc; close all; 

images = dir('*.tiff'); 
 writerObj = VideoWriter('Trot_1Hz_175V.avi');
 writerObj.FrameRate = 30;
 % open the video writer
 open(writerObj);
 % write the frames to the video
 for i = 1:length(images)
     % convert the image to a frame
     frame = im2frame(imread(images(i).name));
     writeVideo(writerObj, frame);
 end
 % close the writer object
 close(writerObj);