function [U, it, crit_hist, omega_hist, itcg, timed_out] = newton_quasi2(U_ini,WEIGHT,K_elast,B,f,Q,mu_0,mu_infty,lambda,p,linear_system_solver,max_runtime_seconds)


    %--------------------------------------------------------------------------
    % The quasi-Newton method with Preconditioner 5 for solution of the system
    %              find U:   F(U)=f
    % Linearized systems are solved by DCG with incomplete Cholesky
    % preconditioner. It is used the fact that the stiffness matrix is
    % constant.
    %
    % Input data:
    %   U_ini - initial choice of U
    %   WEIGHT- weight coefficients of integration points, size(WEIGHT)=(1,n_int)
    %   K_elast - the elastic stiffness matrix, size(K_V)=(3*n_n,3*n_n)
    %   B     - the strain-displacement matrix, size(B)=(6*n_int,3*n_n)
    %   f     - vector of external forces, size(f)=(3,n_n)
    %   Q     - logical array indicating the nodes where the homogeneous
    %           Dirichlet boundary condition is considered, size(Q)=(3,n_n)
    %   mu_0,mu_infty,lambda,p - material parameters
    %
    % Output data:
    %   U - approximation of the solution, size(U)=(3,n_n)
    %   it - number of Newton's iteration
    %   crit_hist  - evolution of the stopping criterion
    %   omega_hist - evolution of the damped coefficients
    %
    %--------------------------------------------------------------------------

    %
    % Auxiliary arrays and initialization
    %
    n_int = length(WEIGHT); % number of integration points
    n_n = size(U_ini, 2); % number of nodes
    dU = zeros(3, n_n); % Newton's increment (in displacement)
    F = zeros(3, n_n); % vector of internal forces
    E = zeros(6, n_int); % values of the strain tensor at integration points
    U = U_ini; % initial approximation of the solution

    % elastic stiffness matrix (related to Precondtioner 2)
    A_N = K_elast(Q, Q);
    A_N = (A_N + A_N') / 2;

    %
    it_max = 200;
    crit_hist = zeros(1, it_max);
    omega_hist = zeros(1, it_max);
    linear_system_solver.setup_preconditioner(A_N);
    %
    % Quasi-Newton's solver (Preconditioner 2)
    %
    if nargin < 12 || isempty(max_runtime_seconds)
        max_runtime_seconds = inf;
    end
    method_timer = tic;

    it = 0; % iteration number
    itcg = 0;
    timed_out = false;

    while true
        if toc(method_timer) >= max_runtime_seconds
            timed_out = true;
            crit_hist = crit_hist(1:it);
            omega_hist = omega_hist(1:max(it - 1, 0));
            fprintf('     Quasi-Newton method 2 timed out after %.2f s', max_runtime_seconds)
            fprintf('\n');
            break
        end

        it = it + 1;

        % solution of the constitutive problem
        E(:) = B * U(:); % strain at integration points
        [S, m, M] = CONSTITUTIVE_PROBLEM.constitutive_problem_quasi2(E,mu_0,mu_infty,lambda,p);
        % solution of the constitutive problem

        % vector of internal forces
        F(:) = B' * reshape(repmat(WEIGHT, 6, 1) .* S, 6 * n_int, 1);
        b_N = f(Q) - F(Q);

        % stopping criterion
        criterion = norm(b_N) / norm(f(Q));
        crit_hist(it) = criterion;

        if criterion < 1e-10
            crit_hist = crit_hist(1:it);
            omega_hist = omega_hist(1:it - 1);
            fprintf(' Quasi-Newton method 2 converges ');
            fprintf(' number of iteration=%d  ', it);
            fprintf(' cumulative cg iterations=%d  ', itcg);
            fprintf(' stopping criterion=%e  ', criterion);
            fprintf('\n');
            break
        end

        % test on number of iteration
        if it == it_max
            omega_hist = omega_hist(1:it - 1);
            fprintf('     Quasi-Newton method 2 converges slowly: stopping criterion=%e  ', criterion)
            fprintf('\n');
            break
        end

        remaining_time = max_runtime_seconds - toc(method_timer);
        if remaining_time <= 0
            timed_out = true;
            crit_hist = crit_hist(1:it);
            omega_hist = omega_hist(1:max(it - 1, 0));
            fprintf('     Quasi-Newton method 2 timed out after %.2f s', max_runtime_seconds)
            fprintf('\n');
            break
        end
        linear_system_solver.max_solve_time_seconds = remaining_time;
        linear_system_solver.A_orthogonalize(A_N);
        [dU(Q), iter] = linear_system_solver.solve(A_N, b_N);
        itcg = itcg + iter;
        if linear_system_solver.last_solve_timed_out
            timed_out = true;
            crit_hist = crit_hist(1:it);
            omega_hist = omega_hist(1:max(it - 1, 0));
            fprintf('     Quasi-Newton method 2 timed out after %.2f s', max_runtime_seconds)
            fprintf('\n');
            break
        end
        linear_system_solver.expand_deflation_basis(dU(Q));

        %fprintf('%d ', iter);

        % update of the unknown vector
        omega = 2 / (M + m);
        omega_hist(it) = omega;
        U = U + omega * dU;

    end % true

end % function
