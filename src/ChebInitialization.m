function [r, k, kh, w] = ChebInitialization(Layers, freq, rmax, dr, c, alpha, ch, alphah)

    w  = 2 * pi * freq;
    r  = dr : dr : rmax;
    
    k = cell(Layers,1);  
    for i = 1 : Layers
        k(i) = {w ./ c{i} .* (1.0 + 1i * alpha{i} / (40.0 * pi * log10(exp(1.0))))};
    end
    kh = w / ch * (1.0 + 1i * alphah / (40.0 * pi * log10(exp(1.0))));
    
end
