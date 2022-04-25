function Plot(kr, psi)

    figure; 
    plot(real(kr), abs(psi), 'r-',  'LineWidth', 3);
    xlabel('kr (1/m)'); ylabel('Magnitude');
    legend('WISpec','box','off','FontSize',16);view(0,90);
    set(gca,'Position',[0.1,0.15,0.8,0.75]);
    set(gcf,'Position',[400,200,800,600]);
    set(gca,'FontSize',16,'LineWidth',2,'box','on');

end