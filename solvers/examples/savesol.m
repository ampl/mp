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

if iter > 0,
	k = 0;
	for i = ceq,
		k = k + 1;
		v(i) = y(k);
		end;
	k = 0;
	for i = clb,
		k = k + 1;
		v(i) = z(k);
		end;
	for i = cub,
		k = k + 1;
		v(i) = -z(k);
		end;
	for i = clu,
		k = k + 2;
		v(i) = z(k-1) - z(k);
		end;
	amplfunc(sprintf( ...
	'nlppd on problem %s: "%s"\nf = %.10g\np_infeas = %.3g\nd_infeas = %.3g\nmu = %.3g',...
		 pname, rc, f, norm(rhs2), norm(rhs1), mu), x, v, solve_result);
else
	fprintf('No solution to save yet! iter = 0\n');
	end;
