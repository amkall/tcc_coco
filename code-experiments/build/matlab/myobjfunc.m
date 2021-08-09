function F = myobjfunc(x, dimension)
    fit = 0;
    for j = 1:dimension
        fit = fit + (j+1)*(x(j)*x(j));
    end
%x1=x(1);
%x2=x(2);
%x3=x(3);
%F=(x3 + 2)*x2*x1^2;

    F = fit;
end