function Vec = ChebDepthSolution(Lb, Ns, Layers, dep, k, rho, kh, rhoh, kr, R)

    Dsave = cell(Layers, 1);
    U     = zeros(sum(Ns+1), sum(Ns+1));
    Vec   = cell(Layers, 1);
    
    n = 1;
    for i = 1 : Layers
        D    = ChebDerivationMatrix(Ns(i)+1);
                                          
        A = 4.0 / ( dep{i}(end) - dep{i}(1) ) ^ 2 * ...
            ConvolutionMatrix(ChebTransFFT(Ns(i), rho{i})) *  D * ...
            ConvolutionMatrix(ChebTransFFT(Ns(i), 1.0 ./ rho{i})) * D + ...
            ConvolutionMatrix( ChebTransFFT(Ns(i), k{i} .^ 2 - kr ^ 2));                                 
                                          
        Dsave(i) = {D};       
        
        U(n:n+Ns(i)-2, n:n+Ns(i)-2) = A(1:Ns(i)-1, 1:Ns(i)-1);
        U(n:n+Ns(i)-2, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = A(1:Ns(i)-1,Ns(i):Ns(i)+1);
        n = n + Ns(i) - 1;
    end  
    
    % boundary condition
    n = 1;
    for i = 1 : Layers - 1
        %displacement potential function 
        Pu = -2 / (dep{i}(end) - dep{i}(1)) * ((-1.0).^(0 : Ns(i))) * Dsave{i};
        Pd =  2 / (dep{i+1}(end) - dep{i+1}(1)) * ones(1, Ns(i+1)+1) * Dsave{i+1};
        
        %sound pressure is continuous
        U(sum(Ns-1)+2*i-1,n:n+Ns(i)-2)   =  rho{i}(end) * (-1.0).^(0:Ns(i)-2);
        U(sum(Ns-1)+2*i-1,sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = rho{i}(end) * (-1.0).^(Ns(i)-1:Ns(i));
        %second interface boundary
        U(sum(Ns-1)+2*i, n:n+Ns(i)-2) = Pu(1:Ns(i)-1);
        U(sum(Ns-1)+2*i, sum(Ns-1)+2*i-1:sum(Ns-1)+2*i) = Pu(Ns(i):Ns(i)+1);
        
        n = n + Ns(i) - 1;
        %sound pressure is continuous
        U(sum(Ns-1)+2*i-1,n:n+Ns(i+1)-2  )   =  -rho{i+1}(1);
        U(sum(Ns-1)+2*i-1,sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2)   =  -rho{i+1}(1);
        %second interface boundary
        U(sum(Ns-1)+2*i,  n:n+Ns(i+1)-2) = Pd(1 : Ns(i+1)-1);
        U(sum(Ns-1)+2*i,  sum(Ns-1)+2*i+1:sum(Ns-1)+2*i+2) = Pd(Ns(i+1): Ns(i+1)+1);

    end  

    %upper boundary, pressure-free boundary
    U(end-1, 1          :    Ns(1)-1) = 1.0;
    U(end-1, sum(Ns-1)+1:sum(Ns-1)+2) = 1.0;
    
    %lower boundary
    %perfectly free / rigid
    Low = (-1.0) .^ (0 : Ns(end));
    if(Lb == 'R')        
        Low = Low * D;%这里的D是最后一个D，恰好没被覆盖
    elseif(Lb == 'A')
        Low = rho{end}(end) * sqrt(kr ^ 2 - kh ^ 2) * Low - ... 
              Low * D * 2.0 / (dep{end}(end) - dep{end}(1)) * rhoh;
    end
    
    U(end, sum(Ns-1)-Ns(end)+2:sum(Ns-1)) = Low(1:Ns(end)-1);
    U(end, end-1:end) = Low(Ns(end):Ns(end)+1);
   
    vec  = U \ R;

    n = 1;
    for i = 1 : Layers
        Vec(i) = {[vec(n:n+Ns(i)-2); vec(sum(Ns-1)+2*i-1:sum(Ns-1)+2*i)]};
        n = n + Ns(i) - 1;
    end

end
