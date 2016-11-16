function [SOL ] = RunEvol3(A)
%matlabpool(4);
%Size = 20;
runs = 200;
lin = linspace(0,1,runs);
clear t;
clear y;
clear toPlotX1;
clear toPlotX2;
clear nchoosekCache;
nchoosekCache = nan(6+1,6+1);
    function [] = update(nchoosekCacheUpdate)
        nchoosekCache=nchoosekCacheUpdate;
    end
callbackfcn = @(nchoosekCacheUpdate)update(nchoosekCacheUpdate);

i=1;
[~,y] = ode45(@(t,y) flexANArg(t,y,A,nchoosekCache,callbackfcn),[0 2], [lin(i) 1-lin(i)]);
[s, ~]=size(y(:,1));
toPlotX1=nan(200,s);
toPlotX1(i,1:s)=y(:,1);
parfor i = 2:runs
    [~,y] = ode45(@(t,y) flexANArg(t,y,A,nchoosekCache,callbackfcn),[0 2], [lin(i) 1-lin(i)]);
    [~, ~]=size(y(:,1));
    toPlotX1(i,:)=y(:,1);
    i
end

toPlotX1(toPlotX1==0) = nan;

for i = 2:runs-1
    [~,b]= size(toPlotX1);
    plot(linspace(0,b,b),toPlotX1(i,:));
    hold on
end
hold off
SOL = solveN(A);

end