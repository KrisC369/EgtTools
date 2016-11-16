function [SOL ] = RunEvol2(A)
%RUNEVOL2 This function produces a phase plots for the 2 action replicator
%dynamics specified by payoff values in A.
%
Size = 20;
runs = 200;
%lin = linspace(0.25,0.26,runs);
lin = linspace(0,1,runs);
%lin = linspace(0.2,0.35,runs);
clear t;
clear y;
clear toPlotX1;
clear toPlotX2;
clear nchoosekCache;
%global nchoosekCache 
nchoosekCache = nan(15+1,15+1);
    function [] = update(nchoosekCacheUpdate)
        nchoosekCache=nchoosekCacheUpdate;
    end
callbackfcn = @(nchoosekCacheUpdate)update(nchoosekCacheUpdate);
for i = 1:runs
    [~,y] = ode45(@(t,y) flexANArg(t,y,A,nchoosekCache,callbackfcn),[0 2], [lin(i) 1-lin(i)]);
    [s, ~]=size(y(:,1));
    toPlotX1(i,1:s)=y(:,1);
    %toPlotX2(i,1:s)=y(:,2);
    i
end

toPlotX1(toPlotX1==0) = nan;
%toPlotX2(toPlotX2==0) = nan;

for i = 2:runs-1
    %plot(toPlotX1(i,:),toPlotX2(i,:));
    [~,b]= size(toPlotX1);
    plot(linspace(0,b,b),toPlotX1(i,:));
    %plot(linspace(0,b,b),toPlotX2(i,:),'r');
    hold on
end
hold off
SOL = solveN(A);
end