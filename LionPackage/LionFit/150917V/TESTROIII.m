i=3;
j=1;

figure(1)
hold on
A=imagesc(ydatacrpdR1{i,j});
plot(x{j}(i,2)+3,x{j}(i,4),'+g');
hold off
figure(2)
B=imagesc(ROI{i,j});