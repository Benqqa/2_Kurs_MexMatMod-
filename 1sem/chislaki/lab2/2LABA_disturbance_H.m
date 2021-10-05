function H1 = disturbance_H(H,b)

    n = length(b);
    for i=1:n
        for j=1:n   
            H1(i,j)=(1+rand*0.01)*H(i,j);
        end 
    end
    
end