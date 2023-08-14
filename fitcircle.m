function [x,r] = fitcircle(A)

    [n,m] = size(A);
    y = [A', ones(m,1)]\sum(A.*A)';
    x = .5*y(1:n);
    r = sqrt(y(n+1) + x'*x);

end