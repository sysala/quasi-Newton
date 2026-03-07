function [U, it, crit_hist, itcg, timed_out] = newton(U_ini, WEIGHT, K_fix, B, f, heter_int, alpha, beta, linear_system_solver, max_runtime_seconds)
%--------------------------------------------------------------------------
% The Newton method for solution of the system
%              find U:   F(U)=f
% Optional arguments:
%   linear_system_solver   - solver object from LINEAR_SOLVERS
%   max_runtime_seconds   - optional timeout in seconds
%
% Output data:
%   U               - approximation of the solution, size(U)=(n_uknown,1)
%   it              - number of Newton's iteration
%   crit_hist       - evaluation of the stopping criterion
%   itcg            - cumulative linear iterations (0 for direct solve)
%   timed_out       - true when timeout is hit
%--------------------------------------------------------------------------

if nargin < 9
    linear_system_solver = [];
end
if nargin < 10 || isempty(max_runtime_seconds)
    max_runtime_seconds = inf;
end
use_iterative_solver = ~isempty(linear_system_solver);

it_max = 100;
tol = 1e-9;
crit_hist = zeros(1, it_max);
itcg = 0;
timed_out = false;
n_int = length(WEIGHT);
E = zeros(2, n_int);  % values of the gradient at integration points
U = U_ini;             % initial approximation of the solution

% Auxiliary arrays
AUX = reshape(1:2 * n_int, 2, n_int);
AUX1 = AUX(:, heter_int);
iD = repmat(AUX1, 2, 1);
jD = kron(AUX1, ones(2, 1));
WEIGHT1 = WEIGHT(heter_int);

it = 0;
method_timer = tic;
while true
    if toc(method_timer) >= max_runtime_seconds
        timed_out = true;
        crit_hist = crit_hist(1:it);
        fprintf('     Newton solver timed out after %.2f s\n', max_runtime_seconds)
        break
    end

    it = it + 1;

    % constitutive operator and its derivative
    E(:) = B * U; % strain at integration points
    [S, DS1] = CONSTITUTIVE_PROBLEM.constitutive_problem(E, heter_int, alpha, beta);
    % tangential stiffness matrix
    vD = repmat(WEIGHT1, 4, 1) .* DS1;
    D_p = sparse(iD(:), jD(:), vD(:), 2 * n_int, 2 * n_int);
    A = K_fix + B' * D_p * B;

    % vector of internal forces
    F = B' * reshape(repmat(WEIGHT, 2, 1) .* S, 2 * n_int, 1);
    b = f - F;

    % stopping criterion
    criterion = norm(b) / norm(f);
    crit_hist(it) = criterion;
    if criterion < tol
        crit_hist = crit_hist(1:it);
        fprintf(' Newton solver converges  ');
        fprintf(' number of iteration=%d  ', it);
        fprintf(' cumulative cg iterations=%d  ', itcg);
        fprintf(' stopping criterion=%e  \n', criterion);
        break
    end

    % test on number of iteration
    if it == it_max
        fprintf('     Newton solver converges slowly: stopping criterion=%e  ', criterion)
        fprintf('\n');
        break
    end

    remaining_time = max_runtime_seconds - toc(method_timer);
    if remaining_time <= 0
        timed_out = true;
        crit_hist = crit_hist(1:it);
        fprintf('     Newton solver timed out after %.2f s\n', max_runtime_seconds)
        break
    end

    if use_iterative_solver
        linear_system_solver.max_solve_time_seconds = remaining_time;
        linear_system_solver.setup_preconditioner(A);
        linear_system_solver.A_orthogonalize(A);
        [dU, iter] = linear_system_solver.solve(A, b);
        itcg = itcg + iter;
        if linear_system_solver.last_solve_timed_out
            timed_out = true;
            crit_hist = crit_hist(1:it);
            fprintf('     Newton solver timed out after %.2f s\n', max_runtime_seconds)
            break
        end
        linear_system_solver.expand_deflation_basis(dU);
    else
        dU = A \ b;
    end

    U = U + dU;
end

end
