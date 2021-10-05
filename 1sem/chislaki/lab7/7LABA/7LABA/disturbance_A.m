function A1 = disturbance_A(A,b)
    n = length(b);
    for i=1:n
        for j=1:n
            A1(i,j)=(1+rand*0.01)*A(i,j);
        end
    end
end