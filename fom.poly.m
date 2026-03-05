function [Z] = form_poly(X,Y,Z,order)
    [M,N] = size(X);
    % Find and remove all missing data points
    xv  = reshape(X,M*N,1);
    yv  = reshape(Y,M*N,1);
    zv  = reshape(Z,M*N,1);
    % Auxiliary variable...
    x = [xv,yv];
    
    % Fit polinomial
    p1 = mypolyfitn(x,zv,order);

    % Evaluate untitled surface on measurement domain
    zp1v = mypolyvaln(p1,x);
    
    Z    = reshape(zp1v,M,N);
