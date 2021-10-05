
a=0.4
b=0.8
def fun(x):
    y=(x**3)-5*(x**8)
    return y
if(fun(a)*fun(b)<0):
    while fun((a+b)/2) !=0:
        if fun(a)*fun(b)<0:
            c=(a+b)/2
            if fun(a)*fun(c)<0:
                b=c
            else:
                a=c
print((a+b)/2)
    
