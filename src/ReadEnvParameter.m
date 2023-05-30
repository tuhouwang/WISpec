function [casename, Src, Layers, Nl, kmax, M, freq, zs, dz, rmax, dr, ...
     tlmin, tlmax, dep, c, rho, alpha, Lb, ch, rhoh, alphah] = ReadEnvParameter(env_file)

    fid           = fopen(env_file);
    casename      = fgetl(fid);
    Src           = fscanf(fid, '%s', 1);
    Layers        = fscanf(fid, '%d', 1);
    Nl            = fscanf(fid, '%d', Layers);   
    kmax          = fscanf(fid, '%f', 1);
    M             = fscanf(fid, '%d', 1);
    freq          = fscanf(fid, '%f', 1);
    zs            = fscanf(fid, '%f', 1);
    dz            = fscanf(fid, '%f', 1);
    rmax          = fscanf(fid, '%f', 1);
    dr            = fscanf(fid, '%f', 1);
    tlmin         = fscanf(fid, '%f', 1);
    tlmax         = fscanf(fid, '%f', 1);
    interface     = fscanf(fid, '%f', Layers);
    nprofile      = fscanf(fid, '%d', Layers);
    
    dep   = cell(Layers,1);
    rho   = cell(Layers,1);
    c     = cell(Layers,1); 
    alpha = cell(Layers,1);
    
    for i = 1 : Layers
        if (i < Layers && interface(i) > 0.0 && ...
            interface(i) < interface(i+1) && nprofile(i) >= 2)
            Profile     = fscanf(fid, '%f %f', [4, nprofile(i)]);
            dep(i)      = {Profile(1, 1:nprofile(i))'};
            c(i)        = {Profile(2, 1:nprofile(i))'};
            rho(i)      = {Profile(3, 1:nprofile(i))'};
            alpha(i)    = {Profile(4, 1:nprofile(i))'};
        elseif(interface(i) > 0.0 && nprofile(i) >= 2)
            Profile     = fscanf(fid, '%f %f', [4, nprofile(i)]);
            dep(i)      = {Profile(1, 1:nprofile(i))'};
            c(i)        = {Profile(2, 1:nprofile(i))'};
            rho(i)      = {Profile(3, 1:nprofile(i))'};
            alpha(i)    = {Profile(4, 1:nprofile(i))'};
        elseif(i == 1)
            Profile     = fscanf(fid, '%f %f', [4, nprofile(i)]);
            dep(i)      = {Profile(1)'};
            c(i)        = {Profile(2)'};
            rho(i)      = {Profile(3)'};
            alpha(i)    = {Profile(4)'};
        else
            error('Error! h must greater than 0 and less than H !');
        end
    end
    
    Lb = fscanf(fid, '%s', 1);
    if (Lb ~= 'V' && Lb ~= 'R' && Lb ~= 'A')
        disp('Error! The lower boundary must be vaccum, rigid or halfspace!');
    end
    
    if (Lb == 'A')
        halfspace = fscanf(fid, '%f', 3);
        ch     = halfspace(1);
        rhoh   = halfspace(2);
        alphah = halfspace(3);
    else
        ch     = 0.0;
        rhoh   = 0.0;
        alphah = 0.0;
    end
    
    % Check the input underwater sound profile
    for i = 1 : Layers
        if (dep{i}(end) ~= interface(i))
            error('Error! input sound profile is unsuitable !');
        end
    end 

    if (zs <= 0  || zs >= interface(end))
        error('zs and zr must be greater than 0 and less than H !');
    end
    
    if (tlmin >= tlmax)
        error('tlmin must less than tlmax !');
    end   
    
    if (Src ~= 'P' && Src ~= 'L')
        disp('Error! The source must be point (P) or line (L)!');
    end

    fclose(fid);

end
