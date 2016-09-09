figure(1)
subplot(1,3,1);
h1=histogram(abs(Lvalnorm{1}),25,'FaceColor','c','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
xlabel('Long Axis Position (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Dif spot location (norm) along long axis');
axis([0 1 0 0.1])
subplot(1,3,2);
h2=histogram(abs(Lvalnorm{2}),25,'FaceColor','y','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
xlabel('Long Axis Position (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('Tus spot location (norm) along long axis')
axis([0 1 0 0.1])
subplot(1,3,3);
h3=histogram(abs(Lvalnorm{3}),25,'FaceColor','r','facealpha',.5,'Normalization','probability');
set(gca,'XTick',[0 1 2 3 4 5 6])
set(gca,'FontSize',12)
xlabel('Long Axis Position (-)','FontWeight','bold');
ylabel('Probability (-)','FontWeight','bold');
title('DnaN spot location (norm) along long axis')
axis([0 1 0 0.1])
