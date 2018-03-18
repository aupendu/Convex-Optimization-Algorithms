% Maximum Power Point Tracking Algorithm (Gradient Search) %
m=1;
E1 = input('Enter the minimum error in function value : ');
Beta = input('Please enter a value for Beta : ');
lk = input('Enter the Step Size : ');
syms l;
Cp = -0.5*((116/((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1))-0.4*Beta*((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1)-5)*exp(-21/((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1));
dk = -diff(Cp);
for m=1:1
    p(m) = input('Enter the Starting Value of Lamda : ');
end
Cpp(1,m) = subs(Cp,l,p(1,m));
dkk(1,m) = subs(dk,l,p(1,m));
p(1,m+1) = p(1,m) + lk*dkk(1,m);
Cpp(1,m+1) = subs(Cp,l,p(1,m+1));
dkk(1,m+1) = subs(dk,l,p(1,m+1));
disp(p(m))
a1 = Cpp(1,m+1)-Cpp(1,m);
while abs(a1)>=E1
    m=m+1;
    Cpp(1,m) = subs(Cp,l,p(1,m));
    dkk(1,m) = subs(dk,l,p(1,m));
    p(1,m+1) = p(1,m)+lk*dkk(1,m);
    Cpp(1,m+1) = subs(Cp,l,p(1,m+1));
    dkk(1,m+1) = subs(dk,l,p(1,m+1));
    disp(p(m))
    a1 = Cpp(1,m+1)-Cpp(1,m);
end
disp('No of Iteration=')
disp(m)
disp('Maximum Value of Cp=')
double(-Cpp(m))
disp(p(1,m))