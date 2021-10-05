function b1 = disturbance_b(b)
    n = length(b);
    for i=1:n
        b1(i)=(1+rand*0.01)*b(i);
    end
    b1=b1';
end