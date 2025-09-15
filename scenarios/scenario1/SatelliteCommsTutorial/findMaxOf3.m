function y = findMaxOf3(a,b,c)
    y = a;
    if y < b
        y = b;
    end

    if y < c
        y = c;
    end
end