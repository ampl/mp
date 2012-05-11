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

function [f,c,d] = evalf(x)

global bl bu cl clb cu cub clu ceq m mp0 p xlc xuc

[f,cc] = amplfunc(x,0);
c = zeros(m,1);
d = zeros(p,1);
v = zeros(mp0,1);
j = 0;
for i = ceq, j = j + 1; c(j) = cc(i) - cl(i); end
j = 0;
for i = clb, j = j + 1; d(j) = cc(i) - cl(i); end
for i = cub, j = j + 1; d(j) = cu(i) - cc(i); end
for i = clu,
	j = j + 2;
	d(j-1) = cc(i) - cl(i);
	d(j) = cu(i) - cc(i);
	end
for i = xlc,
	j = j + 1;
	d(j) = x(i) - bl(i);
	end
for i = xuc,
	j = j + 1;
	d(j) = bu(i) - x(i);
	end
