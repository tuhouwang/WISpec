function Plot(kr, psi)

    figure; hold on;
    plot(real(kr), abs(psi), 'r-',  'LineWidth', 3);
    xlabel('kr (1/m)'); ylabel('Magnitude');
    legend('z=36 m','box','off');
    set(gca,'Position',[0.1,0.15,0.8,0.75]);
    set(gcf,'Position',[400,200,1200,600]);
    set(gca,'FontSize',16,'LineWidth',2,'box','on');

end