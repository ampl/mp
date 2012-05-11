% /****************************************************************
% Copyright (C) 1997-1998 Lucent Technologies
% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and
% its documentation for any purpose and without fee is hereby
% granted, provided that the above copyright notice appear in all
% copies and that both that the copyright notice and this
% permission notice and warranty disclaimer appear in supporting
% documentation, and that the name of Lucent or any of its entities
% not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior
% permission.
%
% LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
% IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
% SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
% IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
% ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.
% ****************************************************************/

% enewt.m

% Much simplified form of nlppd.m (of David Gay, Michael Overton and
% Margaret Wright), omitting the reduced system, checks for negative
% curvature, the watchdog (and watchcat), and separate step lengths
% in x and z, and using an Armijo linesearch to minimize the norm of
% the residual, rather than using a merit function.
%
%  Primal Dual method for solving the nonlinear program
%             min  f(x)               n variables, all nonnegative
%             s.t.  c(x) = 0          m nonlinear equality constraints
%             and   d(x) >= 0.        p nonlinear inequality constraints
%  The optimality conditions are
%                   A(x)'y  + B(x)'z = g(x),  z >= 0.     z is p by 1
%  where g = grad f, A = (grad c)', B = (grad d)',
%  together with the complementarity condition
%                  d(x) .* z = 0
%
%  We shall apply Newton's method to
%		A(x)'y + B(y)'z = g(x)
%		c(x) = 0
%  and		D(x) Z = mu I
%
%  where D=diag(d), Z=diag(z)
%
%  The Newton system reduces to
%  [-W(x,y)    A'(x)   B(x)'  ] [deltax]       [g(x) - A(x)'y - B(x)'z]   [rhs1]
%  [ A(x)      0       0      ] [deltay]   =   [-c(x)]                  = [rhs2]
%  [ Z B(x)    0       D(x)   ] [deltaz]       [diag(mu I - D(x)Z) ]      [rhs3]
%
% Assume init has been run.

maxit1 = iter + maxit;

[f,c,d] = evalf(x);
[g,A,B] = evalg(x);
W = evalw(y,z);
if min(d) < 0,
	fprintf('Negative initial d!\n');
	break;
	end;
rhs1 = g - A'*y - B'*z;
rhs2 = -c;
dtz = d'*z;

if iter == 0,
	nf = 1;
	mu = (norm(rhs1) + dtz) / p;
	z = mu ./ d;
	rhs1 = g - A'*y - B'*z;
	dtz = d'*z;
	end;

rhs3 = mu - d.*z;
rhs = [rhs1; rhs2; rhs3];

K = W + B'*diag(z./d)*B; % inequalities are incorporated in Hessian K

rc = '?';
solve_result = -1;
hiter = iter;
lambda = 0;
lssteps = 0;
starting = 1;

%
% main loop
%
while iter <=  maxit1,	%{
    if display,
	if iter == hiter,
	  hiter = hiter + 10;
          fprintf('iter  ls   lambda    p_infeas     d_infeas    mu       dtz       f\n')
          end;
	if starting,
		fprintf('                  ');
		starting = 0;
	else
		fprintf(' %3d %3d  %8.1e', iter, lssteps, lambda);
		end;
	fprintf(' %11.3e  %11.3e  %8.1e %8.1e %11.3e\n', ...
		norm(rhs2), norm(rhs1), mu, dtz, f);
	end;
    rhsnorm = norm(rhs);
    if rhsnorm < tol & mu < tol,
	rc = '*';
	solve_result = 0;
	break;
	end;
    if iter >= maxit1
	rc = 'I';
	solve_result = 400;
	break;
	end;
    iter = iter + 1;

%
% Set up block matrix to solve for p-d steps.
%
	block = [ -W              A'           B';
	           A              zeros(m,m)   zeros(m,p);
	           diag(z)*B      zeros(p,m)   diag(d)   ];
	delta = block\rhs;
	deltax = delta(1:n);
	deltay = delta(n+1:n+m);
	deltaz = delta(n+m+1:n+m+p);
%
% Estimate step to boundary of nonlinear inequalities.
%
    taumod = 1 - min([1-tau,mu,100*mu*mu]);
    xratio = max([-B*deltax;-deltaz] ./ [d; z]);
    lambda = 1;
    if xratio > 0,
	lambda = min([1, tau/xratio]);
	end;

%
% Line search to reduce norm(rhs)...
%
    oldnorm = norm([rhs1; rhs2]);
    for lssteps = 1:maxlssteps; %{
	xnew = x + lambda*deltax;
	ynew = y + lambda*deltay;
	znew = z + lambda*deltaz;

	[f,c,d] = evalf(xnew);
	nf = nf + 1;
	if d > 0,
		[g,A,B] = evalg(xnew);
		rhs1 = g - A'*ynew - B'*znew;
		rhs2 = -c;
		rhs3 = mu - d.*znew;
		if oldnorm - norm([rhs1; rhs2]) > lambda*lstol*oldnorm,
			break;
			end;
		end;
	lambda = lambda/2;
	end; %} of line search loop

   if lssteps == maxlssteps,
	fprintf('\n**** Too many line search steps.  Sayonara!\n');
	rc = 'S';
	solve_result = 500;
	break;
	end;

    x = xnew;	% c,d,A,B already contain the new values
    y = ynew;
    z = znew;
    W = evalw(y,z);
%
% Update mu
%

    dtz = d'*z;

    if lssteps < 3,
	rhs = [rhs1; rhs2; rhs3 ];
	if norm(rhs) < 10*mu,
		if mu < 1.0e-4,
 	  		mu = 10*mu^2;
	   	else
	 		mu = mu*0.1;
	   	end;
	else
		mu = 0.9*mu;
	end;
     end;

    rhs3 = mu - d.*znew;
    rhs = [rhs1; rhs2; rhs3 ];
    end; %} of main loop

fprintf('pname   rc  iter    nf      mu          rhs1        rhs2        rhs3        f\n');
fprintf('%-8s %s  %-8d%-8d%-12.3g%-12.3g%-12.3g%-12.3g%-15g\n',...
	pname,rc,iter,nf,mu,norm(rhs1),norm(rhs2),norm(rhs3),f);
