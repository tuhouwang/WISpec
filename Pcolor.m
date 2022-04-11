function Pcolor(r, z, tl, casename, tlmin, tlmax)

    figure;
    pcolor(r, z, tl);
    colormap(flipud(jet));caxis([tlmin tlmax]); 
    shading flat; view(0, -90);title(casename);
    xlabel('Range (m)'); ylabel('Depth (m)');
    colorbar('YDir','Reverse','FontSize', 16);
    set(gca,'Position',[0.1,0.15,0.75,0.75],'FontSize',16);
    set(gcf,'Position',[100,100,800,800]); 

end