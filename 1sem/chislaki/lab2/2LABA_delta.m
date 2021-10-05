function s = delta(M,b,n)
    n=length(b);
    s=0;
    for i=1:n
        s=s+M(i).*M(i);
    end
    s=sqrt(s)
end
        