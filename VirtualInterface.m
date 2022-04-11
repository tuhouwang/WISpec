function [dep, c, rho, alpha, Layers, Ns, R] = VirtualInterface(dep, c, rho, alpha, zs, Layers, Ns)

    for i = 1 : Layers
         if(zs > dep{i}(1) && zs < dep{i}(end) )  
             break;
         end
    end
    
    for m = 1 : length(dep{i})
         if(zs > dep{i}(m) && zs < dep{i}(m+1) )  
             break;
         end
    end
    
    dep1   = cell(Layers+1,1);
    rho1   = cell(Layers+1,1);
    c1     = cell(Layers+1,1); 
    alpha1 = cell(Layers+1,1);
    
    if(i == 1)
        dep1  (i) = {[dep{i}(1:m);zs]};
        rho1  (i) = {interp1(dep{i},rho{i},  dep1{i},'linear','extrap')};
        c1    (i) = {interp1(dep{i},c{i},    dep1{i},'linear','extrap')};
        alpha1(i) = {interp1(dep{i},alpha{i},dep1{i},'linear','extrap')};
        dep1  (i+1) = {[zs;dep{i}(m+1:end)]};
        rho1  (i+1) = {interp1(dep{i},rho{i},  dep1{i+1},'linear','extrap')};
        c1    (i+1) = {interp1(dep{i},c{i},    dep1{i+1},'linear','extrap')};
        alpha1(i+1) = {interp1(dep{i},alpha{i},dep1{i+1},'linear','extrap')};
    else
        dep1  (1:i-1) = dep  (1:i-1);
        rho1  (1:i-1) = rho  (1:i-1);
        c1    (1:i-1) = c    (1:i-1);
        alpha1(1:i-1) = alpha(1:i-1);
        
        dep1  (i) = {[dep{i}(1:m);zs]};
        rho1  (i) = {interp1(dep{i},rho{i},  dep1{i},'linear','extrap')};
        c1    (i) = {interp1(dep{i},c{i},    dep1{i},'linear','extrap')};
        alpha1(i) = {interp1(dep{i},alpha{i},dep1{i},'linear','extrap')};
        dep1  (i+1) = {[zs;dep{i}(m+1:end)]};
        rho1  (i+1) = {interp1(dep{i},rho{i},  dep1{i+1},'linear','extrap')};
        c1    (i+1) = {interp1(dep{i},c{i},    dep1{i+1},'linear','extrap')};
        alpha1(i+1) = {interp1(dep{i},alpha{i},dep1{i+1},'linear','extrap')};
        
    end
    
    dep1  (i+2:end) = dep  (i+1:end);
    rho1  (i+2:end) = rho  (i+1:end);
    c1    (i+2:end) = c    (i+1:end);
    alpha1(i+2:end) = alpha(i+1:end);

    
    dep    = dep1;
    c      = c1;
    rho    = rho1;
    alpha  = alpha1;
    Layers = Layers + 1;
    Ns     = [Ns(1:i);Ns(i);Ns(i+1:end)];
    
    R     = zeros(sum(Ns+1), 1);
    R(sum(Ns-1)+2*i) = - 0.5 / pi;

end