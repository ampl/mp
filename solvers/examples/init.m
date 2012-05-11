% /****************************************************************
% Copyright (C) 1997 Lucent Technologies
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

% Initialize for enewt.m
% pname must have been set to the name (or stub) of
% the problem's .nl file.

global bl bu cl clb cu cub clu ceq m mp0 n p pname xlc xuc

if ischar(pname) == 0,
	error('No ''stub'' assigned to pname.');
	end;

[x,bl,bu,v,cl,cu] = amplfunc(pname);

n = length(x);
mp0 = length(cl);
xlc = [];	% i: constraints added for x(i) >= bl(i)
xuc = [];	% i: constraints added for x(i) <= bu(i)
for i = 1:n,
	if bu(i) < 1e20,  xuc = [xuc i]; end
	if bl(i) > -1e20, xlc = [xlc i]; end
	end

% Prepare to turn AMPL constraints cl(i) <= ca(i) <= cu(i)
% into c(x) == 0, d(x) >= 0...

y = [];		% multipliers for c(x) == 0
z = [];		% multipliers for d(x) >= 0

clb = [];	% i: ca(i) >= cl(i)
cub = [];	% i: ca(i) <= cu(i)
clu = [];	% i: cl(i) <= ca(i) <= cu(i)
ceq = [];	% i: ca(i) == cl(i) == cu(i)
for i = 1:mp0,
	if cu(i) < 1e20,
		if cl(i) > -1e20,
			if cl(i) >= cu(i),
				ceq = [ceq i];
				y = [y; -v(i)];
			else
				clu = [clu i];
			end
		else
			cub = [cub i];
			end
	else
		clb = [clb i];
		end
	end
for i = clb,
	t = v(i);
	if t <= 0, t = 1; end
	z = [z; t];
	end
for i = cub,
	t = -v(i);
	if t <= 0, t = 1; end
	z = [z; t];
	end
for i = clu,
	t = 0.5*v(i);
	if t <= 0, t = 0.5; end
	z = [z; -t; t];
	end
for i = [xlc xuc], z = [z; 1]; end
m = length(y);
p = length(z);

display = 1;
iter = 0;
maxit = 50;
lstol = 0.1;		% require at least lstol of linear norm reduction
maxlssteps = 20;	% maximum number of linesearch steps
tau = 0.99;		% factor for step to boundary
tol = 1e-8;

[f,c,d] = evalf(x);
fprintf('n = %d, m = %d, p = %d; dmin = %.3g\n', n, m, p, min(d));
