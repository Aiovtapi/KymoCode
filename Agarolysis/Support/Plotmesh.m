close all

plotbacs = 1:80;


figure
imagesc(double(flimg))
hold on


for i = plotbacs;
    plot(Bettermesh{i,1}(:,1),Bettermesh{i,1}(:,2),'r')
    plot(Bettermesh{i,1}(:,3),Bettermesh{i,1}(:,4),'r')
end

flimgsize = size(flimg);
bound = 2;

nflimg = zeros(flimgsize(1)+2*bound,flimgsize(2)+2*bound);
nflimg(bound+1:end-bound,bound+1:end-bound)=flimg;

figure
imagesc(double(nflimg))
hold on


for i = plotbacs;
    plot(BBettermesh{i,1}(:,1),BBettermesh{i,1}(:,2),'r')
    plot(BBettermesh{i,1}(:,3),BBettermesh{i,1}(:,4),'r')
end




imagedb = double(flimg);


for j =1:80;
    lowx = round(Cellbox(j,1,1));
    lowy = round(Cellbox(j,1,3));
    
    mask = Bacmask{j,1};
    msize = size(mask);
    
    highx = lowx + msize(2) - 1;
    highy = lowy + msize(1) - 1;
    
    imagedb(lowy:highy,lowx:highx) = imagedb(lowy:highy,lowx:highx).*mask;
end
    
figure
imagesc(imagedb)







