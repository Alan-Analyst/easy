# easy

function x = gaussSeidelMatrixForm(A, b, x0, maxIter, tol)
    % Matrix form: x^(k+1) = M * x^(k) + c

    D = diag(diag(A));
    L = tril(A, -1);
    U = triu(A, 1);

    M = -(D + L) \ U;  % Equivalent to inv(D + L) * (-U)
    c = (D + L) \ b;

    x = x0(:);  % Ensure column vector
    for k = 1:maxIter
        x_new = M * x + c;
        if norm(x_new - x, inf) < tol
            fprintf('Converged in %d iterations.\n', k);
            x = x_new;
            return;
        end
        x = x_new;
    end

    warning('Maximum iterations reached without convergence.');
end

function Jacobi_matrix(A, x0, b)
D = diag(diag(A));
U = triu(A, 1);
L = tril(A, -1);

M = inv(D)*(U + L);
C = inv(D)*b;

x = x0(:);

for i=1:2
    x_new = M*x + C

end

function Jacobi(b, maxIter, x0, tol)
A = input('Enter your Matrix');
[rows, cols] = size(A);

if rows ~= cols
    disp('Matrix must be square')
end

% Check for diagonal dominanece 
IsDiagonallyDominant = true;

for i=1:rows
    diagElem = abs(A(i,i));
    offDiagElem = sum(abs(A(i, :))) - diagElem;
    if diagElem < offDiagElem
        IsDiagonallyDominant = false;
        break;
    end
end

n = length(b);
x = x0;

if IsDiagonallyDominant == true
    disp('Matrix is diagonally dominant.')

    for k = 1:maxIter
        x_new = zeros(n,1);
        for i = 1:n
            sumExceptI = A(i, :) * x - A(i,i) * x(i);
            x_new(i) = (b(i) - sumExceptI) / A(i,i)
        end

          % Check for convergence
        if norm(x_new - x, inf) < tol
            fprintf('Converged in %d iterations.\n', k);
            x = x_new;
            return;
        end
        x = x_new;

    end

        

else
    disp('Matrix must be diagonally dominant.')
end

end

function x = gauss(A, b, x0, maxIter, tol)

    n = length(b);
    x = x0(:);  % Ensure column vector

    for k = 1:maxIter
        x_old = x;
        for i = 1:n
            % Calculate sum using updated x values where available
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) / A(i,i);
        end

        % Convergence check
        if norm(x - x_old, inf) < tol
            fprintf('Converged in %d iterations.\n', k);
            return;
        end
    end

    warning('Maximum iterations reached without convergence.');
end

function table = simpleForwardDiff(x, y)
    % SIMPLEFORWARDDIFF creates a forward difference table
    % x - list of x values (equally spaced)
    % y - list of f(x) values

    n = length(y);
    table = zeros(n, n);      % Make an n x n table
    table(:,1) = y(:);        % First column = original y values

    % Fill in the difference columns
    for col = 2:n
        for row = 1:n-col+1
            table(row, col) = table(row+1, col-1) - table(row, col-1);
        end
    end
end

function NR(f, g, x0, y0, iter)
syms x
syms y


fx = diff(f, x);
fy = diff(f, y);
gx = diff(g, x);
gy = diff(g, y);

J = fx*gy-fy*gx;

if subs(J, {x, y}, {x0, y0}) == 0
    error('division by zero')
else 
    x1 = x0 - (subs(g*fy-f*gy, {x, y}, {x0, y0})/subs(J, {x, y}, {x0, y0}));
    y1 = y0 - (subs(g*fx-f*gx, {x, y}, {x0, y0}))/subs(J, {x, y}, {x0, y0});
    solution_vector = zeros(iter, 2);  % Preallocate for speed


    for i=1:iter

        x0 = x1;
        y0 = y1;     
        x1 = x0 - (subs(g*fy-f*gy, {x, y}, {x0, y0})/subs(J, {x, y}, {x0, y0}));
        y1 = y0 - (subs(g*fx-f*gx, {x, y}, {x0, y0}))/subs(J, {x, y}, {x0, y0});
        solution_vector(i, :) = [double(x1), double(y1)];



    end
    
    disp(solution_vector)

end

