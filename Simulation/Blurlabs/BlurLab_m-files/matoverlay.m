function matout=matoverlay(mat1,mat2,x,y)
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
%mat1=(rand(500,500)<0.001).*rand(500,500);
%mat2=mat2gray(fspecial('gaussian',[11,11],3));
%
%[y,x]=find(mat1>0);
%
%I0=zeros(size(mat1));
%for i=1:length(x)
%   I0=I0+mat1(y(i),x(i))*matoverlay(mat1,mat2,x(i),y(i));
%end
%
%I1=conv2(mat1,mat2,'same');
%
%figure;
%colormap(hot)
%subplot(1,3,1)
%imagesc(mat1)
%axis equal
%axis tight
%title('mat1')
%
%subplot(1,3,2)
%imagesc(mat2)
%axis equal
%axis tight
%title('mat2')
%
%subplot(1,3,3)
%imagesc(I0)
%axis equal
%axis tight
%title('Sparse Convolution of mat1 and mat2')
%


%get sizes of input matrices
[sy1,sx1]=size(mat1);
[sy2,sx2]=size(mat2);

%check the values of x and y
if or(x<0.5,y<0.5)
    error('The coordinate point must round to a value greater than (1,1).')
end
if or(x>=(sx1+1/2),y>=(sy1+1/2))
    error('The coordinate point must round to a valid value.')
end

%convert coordinate point to indices
xpt=round(x);
ypt=round(y);

%find indices of mat2 center
xmid=round(sx2/2);
ymid=round(sy2/2);

%create temporary new mat1 with larger mat2 sized borders.
mattemp=zeros(sy1+2*ymid+2,sx1+2*xmid+2);
%mattemp(ymid+1:ymid+sy1,xmid+1:xmid+sx1)=mat1;
leftpt=xpt+1;
botpt=ypt+1;
mattemp(botpt:botpt+sy2-1,leftpt:leftpt+sx2-1)=mat2;

%set output matrix
matout=mattemp(ymid+1:ymid+sy1,xmid+1:xmid+sx1);


