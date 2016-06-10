function matout=matoverlay_sub(mat1,mat2,x,y)
%Tristan Ursell
%Matrix Overlay function
%March 2011
%
%matout=matoverlay(mat1,mat2,x,y);
%
%This function takes an input matrix 'mat1' and creates an image of the
%matrix 'mat2' at the position (x,y), padded by zeros.  The output matrix
%will have the same size at mat1.
%
%Essentially this is performing a sparse, fully valid convolution of the
%point (x,y) and mat2 with the output size of mat1.  The point (x,y) uses
%the imaging convention for the coordinate axes. For best results, the size
%of mat2 should be odd (but it can be even).
%
%Example:
%
%N=200;
%x=1+499*rand(1,N);
%y=1+499*rand(1,N);
%
%mat1=zeros(500,500);
%
%mat2=mat2gray(fspecial('gaussian',[11,11],3));
%
%I0=zeros(size(mat1));
%ints=rand(1,N);
%for i=1:N
%   I0=I0+ints(i)*matoverlay_sub(mat1,mat2,x(i),y(i));
%end
%
%figure;
%colormap(hot)
%subplot(1,2,1)
%imagesc(mat2)
%axis equal
%axis tight
%title('mat2')
%
%subplot(1,2,2)
%hold on
%imagesc(I0)
%plot(x-1/2,y-1/2,'bo')
%axis equal
%axis tight
%title('Sparse Convolution of mat1 and mat2')
%


%get sizes of input matrices
[sy1,sx1]=size(mat1);
[sy2,sx2]=size(mat2);

%check the values of x and y
if or(floor(x)<1,floor(y)<1)
    error('The coordinate point must round to a value greater than (1,1).')
end
if or(ceil(x)>sx1,ceil(y)>sy1)
    error('The coordinate point must round to a valid value.')
end

%find indices of mat2 center
xmid=round(sx2/2);
ymid=round(sy2/2);

%turn on to revert to matoverlay
%x=round(x);
%y=round(y);

if floor(x)==x
    x=x+1e-10;
end
if floor(y)==y
    y=y+1e-10;
end

%calculate weights (starting from top left, clockwise)
w1=(x-floor(x))*(y-floor(y));
w2=(ceil(x)-x)*(y-floor(y));
w3=(ceil(x)-x)*(ceil(y)-y);
w4=(x-floor(x))*(ceil(y)-y);

%construct weighted kernel
kernel=zeros(size(mat2));
kernel(2:end,2:end)=kernel(2:end,2:end)+w1*mat2(1:end-1,1:end-1);
kernel(2:end,:)=kernel(2:end,:)+w2*mat2(1:end-1,:);
kernel=kernel+w3*mat2;
kernel(:,2:end)=kernel(:,2:end)+w4*mat2(:,1:end-1);

%create temporary new mat1 with larger mat2 sized borders.
mattemp=zeros(sy1+2*ymid+2,sx1+2*xmid+2);

%place down green's function
mattemp(floor(y)+1:floor(y)+sy2,floor(x)+1:floor(x)+sx2)=kernel;

%set output matrix
matout=mattemp(ymid+1:ymid+sy1,xmid+1:xmid+sx1);


