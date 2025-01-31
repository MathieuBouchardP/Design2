close all
clear all
clc

F(length(zoom)) = struct('cdata',[],'colormap',[]);
writerObj = VideoWriter('TemperatureDistribution2D.avi');
open(writerObj);

