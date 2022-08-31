function Plot_error(r, tl)

 figure;
 set(gcf,'Position',[400,200,800,600]);

 subplot('Position',[0.1,0.6,0.75,0.3]); 
 plot(r(1:300),tl(1:300), 'm-', 'Linewidth',2);
 xlabel('Range (m)'); ylabel('Error (dB)');
 set(gca,'FontSize',16,'LineWidth',2,'box','on');

 
 subplot('Position',[0.1,0.15,0.75,0.3]); 
 plot(r(300:end),tl(300:end), 'm-', 'Linewidth',2);
 xlabel('Range (m)'); ylabel('Error (dB)');
 axis([300 3000 0 0.1])
 set(gca,'FontSize',16,'LineWidth',2,'box','on');
 
end