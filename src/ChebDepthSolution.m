function Vec = ChebDepthSolution(Lb, Nl, Layers, dep, k, rho, kh, rhoh, kr, R)

    Dsave = cell(Layers, 1);
    U     = zeros(sum(Nl+1), sum(Nl+1));
    Vec   = cell(Layers, 1);
    
    n = 1;
    for i = 1 : Layers
        D    = ChebDerivationMatrix(Nl(i)+1);
                                          
        A = 4.0 / ( dep{i}(end) - dep{i}(1) ) ^ 2 * ...
            ConvolutionMatrix(ChebTransFFT(Nl(i), rho{i})) *  D * ...
            ConvolutionMatrix(ChebTransFFT(Nl(i), 1.0 ./ rho{i})) * D + ...
            ConvolutionMatrix( ChebTransFFT(Nl(i), k{i} .^ 2 - kr ^ 2));                                 
                                          
        Dsave(i) = {D};       
        
        U(n:n+Nl(i)-2, n:n+Nl(i)-2) = A(1:Nl(i)-1, 1:Nl(i)-1);
        U(n:n+Nl(i)-2, sum(Nl-1)+2*i-1:sum(Nl-1)+2*i) = A(1:Nl(i)-1,Nl(i):Nl(i)+1);
        n = n + Nl(i) - 1;
    end  
    
    % boundary condition.
    n = 1;
    for i = 1 : Layers - 1
        % displacement potential function.
        Pu = -2 / (dep{i}(end) - dep{i}(1)) * ((-1.0).^(0 : Nl(i)))  * Dsave{i};
        Pd =  2 / (dep{i+1}(end) - dep{i+1}(1)) * ones(1, Nl(i+1)+1) * Dsave{i+1};
        
        % sound pressure is continuous.
        U(sum(Nl-1)+2*i-1,n:n+Nl(i)-2)   =  rho{i}(end) * (-1.0).^(0:Nl(i)-2);
        U(sum(Nl-1)+2*i-1,sum(Nl-1)+2*i-1:sum(Nl-1)+2*i) = rho{i}(end) * (-1.0).^(Nl(i)-1:Nl(i));
        % second interface boundary.
        U(sum(Nl-1)+2*i, n:n+Nl(i)-2) = Pu(1:Nl(i)-1);
        U(sum(Nl-1)+2*i, sum(Nl-1)+2*i-1:sum(Nl-1)+2*i) = Pu(Nl(i):Nl(i)+1);
        
        n = n + Nl(i) - 1;
        % sound pressure is continuous.
        U(sum(Nl-1)+2*i-1,n:n+Nl(i+1)-2  )   =  -rho{i+1}(1);
        U(sum(Nl-1)+2*i-1,sum(Nl-1)+2*i+1:sum(Nl-1)+2*i+2)   =  -rho{i+1}(1);
        % second interface boundary.
        U(sum(Nl-1)+2*i,  n:n+Nl(i+1)-2) = Pd(1:Nl(i+1)-1);
        U(sum(Nl-1)+2*i,  sum(Nl-1)+2*i+1:sum(Nl-1)+2*i+2) = Pd(Nl(i+1):Nl(i+1)+1);

    end  

    % upper boundary, pressure-free boundary.
    U(end-1, 1          :    Nl(1)-1) = 1.0;
    U(end-1, sum(Nl-1)+1:sum(Nl-1)+2) = 1.0;
    
    % lower boundary.
    % perfectly free / rigid.
    Low = (-1.0) .^ (0 : Nl(end));
    if(Lb == 'R')        
        Low = Low * D; % Last D£¬coincidentally not covered.
    elseif(Lb == 'A')
        Low = rho{end}(end) * sqrt(kr ^ 2 - kh ^ 2) * Low - ... 
              Low * D * 2.0 / (dep{end}(end) - dep{end}(1)) * rhoh;
    end
    
    U(end, sum(Nl-1)-Nl(end)+2:sum(Nl-1)) = Low(1:Nl(end)-1);
    U(end, end-1:end) = Low(Nl(end):Nl(end)+1);
   
    vec  = U \ R;

    n = 1;
    for i = 1 : Layers
        Vec(i) = {[vec(n:n+Nl(i)-2); vec(sum(Nl-1)+2*i-1:sum(Nl-1)+2*i)]};
        n = n + Nl(i) - 1;
    end

end