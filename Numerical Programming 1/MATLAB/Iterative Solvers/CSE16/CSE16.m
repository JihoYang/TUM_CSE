
	Academic License

cd('C:\Users\Folkmar Bornemann\Desktop\SOR')
n = 10000; rng(843); A = sprandsym(n,10/n,1/100,2);
clc
condest(A)

ans =

           162.76159479457

tic, x_chol = A\b; toc
{Undefined function or variable 'b'.
} 
clc
b = ones(n,1);
tic, x_chol = A\b; toc
Elapsed time is 0.692230 seconds.
n = 20000; rng(843); A = sprandsym(n,10/n,1/100,2);
b = ones(n,1);
condest(A)
{Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('lu')" style="font-weight:bold">lu</a>
Out of memory. Type HELP MEMORY for your options.

Error in <a href="matlab:matlab.internal.language.introspective.errorDocCallback('condest', 'C:\Program Files (x86)\MATLAB\R2015a\toolbox\matlab\matfun\condest.m', 48)" style="font-weight:bold">condest</a> (<a href="matlab: opentoline('C:\Program Files (x86)\MATLAB\R2015a\toolbox\matlab\matfun\condest.m',48,0)">line 48</a>)
   [L,U,~,~] = lu(A,'vector');
} 
clc
clear all
n = 20000; rng(843); A = sprandsym(n,10/n,1/100,2);
{Operation terminated by user during <a href="matlab:matlab.internal.language.introspective.errorDocCallback('sprandsym', 'C:\Program Files (x86)\MATLAB\R2015a\toolbox\matlab\sparfun\sprandsym.m', 187)" style="font-weight:bold">sprandsym</a> (<a href="matlab: opentoline('C:\Program Files (x86)\MATLAB\R2015a\toolbox\matlab\sparfun\sprandsym.m',187,0)">line 187</a>)

} 
clear all
n = 10000; rng(843); A = sprandsym(n,10/n,1/100,2);
n = 10000; rng(843); A = sprandsym(n,10/n,1/100,2); b = ones(n,1);
condest(A)

ans =

           162.76159479457

format short g
tic, x_chol = A\b; toc
Elapsed time is 0.696332 seconds.
tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1); toc
Elapsed time is 0.104253 seconds.
tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1.5); toc
Elapsed time is 0.069597 seconds.
tic, [x_cg,flag,relres,iter] = pcg(A,b,1e-10,1000); toc
Elapsed time is 0.076480 seconds.
iter

iter =

    67

clc
n = 10000; rng(843); A = sprandsym(n,10/n,1/100,2); b = ones(n,1);
spy(A), shg
tic, x_chol = A\b; toc
Elapsed time is 0.691305 seconds.
condest(A)

ans =

       162.76

tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1); toc
Elapsed time is 0.097569 seconds.
iter

iter =

   106

norm(x_sor-x_chol)/norm(x_chol)

ans =

   7.3606e-11

tic, [x_sor,relres,iter] = SOR(A,b,1e-10,0.5); toc
Elapsed time is 0.204186 seconds.
iter

iter =

   220

tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1.5); toc
Elapsed time is 0.071937 seconds.
iter

iter =

    69

tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1.7); toc
Elapsed time is 0.092067 seconds.
iter

iter =

    96

tic, [x_sor,relres,iter] = SOR(A,b,1e-10,1.3); toc
Elapsed time is 0.075708 seconds.
iter

iter =

    80

edit SOR
tic, [x_cg,flag,relres,iter] = pcg(A,b,1e-10,1000); toc
Elapsed time is 0.054714 seconds.
iter

iter =

    67

norm(x_cg-x_chol)/norm(x_chol)

ans =

   3.4523e-11

n = 10000; rng(843); A = sprandsym(n,10/n,1/10000,2); b = ones(n,1);
tic, [x_cg,flag,relres,iter] = pcg(A,b,1e-10,1000); toc
Elapsed time is 0.301278 seconds.
iter

iter =

   724

