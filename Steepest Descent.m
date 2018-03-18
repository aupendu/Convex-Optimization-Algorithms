% Maximum Power Point Tracking Algorithm (Steepest Descent Method) %
m=1;
E1 = input('Enter the minimum error in function value');
Beta = input('Please enter a value for Beta:');
syms l;
Cp = -0.5*((116/((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1))-0.4*Beta*((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1)-5)*exp(-21/((1/(l+0.08*Beta)-(0.035/(Beta^3+1))))^(-1));
dk = -diff(Cp);
lk = (diff(diff(Cp)))^(-1);
for m=1:1
    p(m) = input('Enter the Starting Value of Lamda');
end
Cpp(1,m) = subs(Cp,l,p(1,m));
dkk(1,m) = subs(dk,l,p(1,m));
lkk(1,m) = subs(lk,l,p(1,m));
p(1,m+1) = p(1,m) + lkk(1,m)*dkk(1,m);
Cpp(1,m+1) = subs(Cp,l,p(1,m+1));
dkk(1,m+1) = subs(dk,l,p(1,m+1));
a1 = Cpp(1,m+1)-Cpp(1,m);
while abs(a1)>=E1
    m=m+1;
    Cpp(1,m) = subs(Cp,l,p(1,m));
    dkk(1,m) = subs(dk,l,p(1,m));
    lkk(1,m) = subs(lk,l,p(1,m));
  if dkk(1,m)>0
        if lkk(1,m)>0
    p(1,m+1) = p(1,m) + lkk(1,m)*dkk(1,m);
        end
        if lkk(1,m)<0
    p(1,m+1) = p(1,m) - lkk(1,m)*dkk(1,m);
        end
  end
  if dkk(1,m)<0
        if lkk(1,m)>0
    p(1,m+1) = p(1,m) + lkk(1,m)*dkk(1,m);
        end
        if lkk(1,m)<0
    p(1,m+1) = p(1,m) - lkk(1,m)*dkk(1,m);
        end
  end
    Cpp(1,m+1) = subs(Cp,l,p(1,m+1));
    dkk(1,m+1) = subs(dk,l,p(1,m+1));
    disp(p(m))
    a1 = Cpp(1,m+1)-Cpp(1,m);
end
disp('No of Iteration=')
disp(m)
disp('Maximum Value of Cp=')
double(-Cpp(m))
disp(p(m))