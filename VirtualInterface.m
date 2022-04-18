function [dep, c, rho, alpha, Layers, Ns, R, s] = VirtualInterface(dep, c, rho, alpha, zs, Layers, Ns)

    for s = 1 : Layers
         if(zs > dep{s}(1) && zs < dep{s}(end) )  
             break;
         end
    end
    
    for m = 1 : length(dep{s})
         if(zs >= dep{s}(m) && zs < dep{s}(m+1) )  
             break;
         end
    end
    
    dep1   = cell(Layers+1,1);
    rho1   = cell(Layers+1,1);
    c1     = cell(Layers+1,1); 
    alpha1 = cell(Layers+1,1);
    
    if(s ~= 1)
        dep1  (1:s-1) = dep  (1:s-1);
        rho1  (1:s-1) = rho  (1:s-1);
        c1    (1:s-1) = c    (1:s-1);
        alpha1(1:s-1) = alpha(1:s-1);     
    end
    
    if(zs == dep{s}(m))
       dep1  (s) = {dep{s}(1:m)};
       rho1  (s) = {rho{s}(1:m)};
       c1    (s) = {c{s}(1:m)};
       alpha1(s) = {alpha{s}(1:m)};
       dep1  (s+1) = {dep{s}(m:end)};
       rho1  (s+1) = {rho{s}(m:end)};
       c1    (s+1) = {c{s}(m:end)};
       alpha1(s+1) = {alpha{s}(m:end)};
    else
       dep1  (s) = {[dep{s}(1:m);zs]};
       rho1  (s) = {interp1(dep{s},rho{s},  dep1{s},'linear','extrap')};
       c1    (s) = {interp1(dep{s},c{s},    dep1{s},'linear','extrap')};
       alpha1(s) = {interp1(dep{s},alpha{s},dep1{s},'linear','extrap')};
       dep1  (s+1) = {[zs;dep{s}(m+1:end)]};
       rho1  (s+1) = {interp1(dep{s},rho{s},  dep1{s+1},'linear','extrap')};
       c1    (s+1) = {interp1(dep{s},c{s},    dep1{s+1},'linear','extrap')};
       alpha1(s+1) = {interp1(dep{s},alpha{s},dep1{s+1},'linear','extrap')};
    end  
    
    dep1  (s+2:end) = dep  (s+1:end);
    rho1  (s+2:end) = rho  (s+1:end);
    c1    (s+2:end) = c    (s+1:end);
    alpha1(s+2:end) = alpha(s+1:end);
    
    dep    = dep1;
    c      = c1;
    rho    = rho1;
    alpha  = alpha1;
    Layers = Layers + 1;
    Ns     = [Ns(1:s); Ns(s); Ns(s+1:end)];
    
    R     = zeros(sum(Ns+1), 1);
    R(sum(Ns-1)+2*s) = 0.5 / pi;

end