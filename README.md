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

