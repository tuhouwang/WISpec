function [c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Nl)

    for i = 1 : Layers
        x = cos( (0 : Nl(i)) * pi / Nl(i) )';
        z = ( ( dep{i}(1) + dep{i}(end) ) / ( dep{i}(end) - dep{i}(1) ) - x )...
                                            * ( dep{i}(end) - dep{i}(1) ) / 2.0;

        c(i)     = {interp1(dep{i}, c{i},     z, 'linear', 'extrap')};
        rho(i)   = {interp1(dep{i}, rho{i},   z, 'linear', 'extrap')};
        alpha(i) = {interp1(dep{i}, alpha{i}, z, 'linear', 'extrap')};
    end
    
end
