function [U_3,Q]=boundary_BC3(coord,surf1,surf3,surf4,surf5,surf6)

% =========================================================================
%
%  Paper-style boundary conditions for the tapered bar benchmark:
%    - prescribed constant axial velocity on the inlet face z = 0
%    - prescribed zero axial velocity on the shell faces
%    - x/y components remain free on the shell
%
%  This mirrors the wording in paper.tex for Example 2. The inlet and shell
%  conditions are formally incompatible on their shared edge; this helper
%  gives priority to the inlet value there.
%
% ======================================================================
%

Q = true(size(coord));
Q(3, surf1(:)) = 0;
Q(3, surf3(:)) = 0;
Q(3, surf4(:)) = 0;
Q(3, surf5(:)) = 0;
Q(3, surf6(:)) = 0;

U_3 = zeros(1, size(coord, 2));
Q_s1 = false(1, size(coord, 2));
Q_s1(surf1(:)) = true;
U_3(Q_s1) = 1;

end
