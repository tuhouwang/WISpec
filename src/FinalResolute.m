function [z, rhoi] = FinalResolute(dep, dz, rho, Layers)

    z    = 0 : dz : dep{end}(end);
    zt   = [];
    rhoi = [];
    zi   = cell(Layers, 1);
    temp = cell(Layers, 1);
    for m = 1 : Layers
        zi  (m) = {dep{m}(1) : dz : dep{m}(end)};
        temp(m) = {interp1(dep{m},rho{m},zi{m},'linear','extrap')};       
    end
    
    for m = 1 : Layers
        if(isempty(rhoi) == 1)
           zt   = zi{1};
           rhoi = temp{m};
        elseif(zt(end) == zi{m}(1))   
           zt   = [zt,     zi{m}(2:end)];  
           rhoi = [rhoi, temp{m}(2:end)];
        else
           zt   = [zt,     zi{m}];
           rhoi = [rhoi, temp{m}];        
        end  
    end
    
    rhoi = interp1(zt, rhoi, z, 'linear', 'extrap');

end