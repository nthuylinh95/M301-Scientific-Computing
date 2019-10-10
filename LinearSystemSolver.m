% Nguyen Thuy Linh
% MATH 301 MIDTERM
% Program for solving a linear system of equations of the form Ax=b
% LU factorization without pivots,like in the example from class. 
% Forward substitution and then back substitution.
% Example matrices provided.

function x = LinearSystemSolver(A, b)
% A = [1 2 -3 4; 2 2 -2 3; 0 1 1 0; 1 -1 1 -2] % 4x4
% b = [12; 10; -1; -4] % 4x4
% A = [2 1 3; 2 6 8; 6 8 18]; % 3x3
% b = [1; 3; 5]; % 3x3
% A = [1 1; -3 1]; % 2x2
% b = [6; 2]; % 2x2
[m,m] = size(A);
L=eye(m);
% We first perform LU factorization using the code from previous homework.
% We create a loop that allows us to get entries (2,1), (3,1), (3,2) etc.
for M = 2:m % Start from second row.
    for N = 1:(M-1)
        E=eye(m); % Creating an identity matrix E.
        E(M,N)=-A(M,N)/A(N,N); % Making zeros by dividing pivot columns by pivot entries.
        A = E*A; % Recording the U matrix, which is the result of performed operations.
        L = E*L; % Records a matrix whose inverse will give us the L we are looking for.
    end
end
L=inv(L);
U=A;
% We obtained L and U, such that A=LU. Now we want to solve Ax=b. We first 
% rewrite this equation as LUx=b. Multiplying both sides from the left by
% L^(-1), we obtain Ux=L^(-1)b. We let c=Ux so that we get c=L^(-1)b. 
% This means Lc=b, and we will solve this first (because we can) by using 
% forward substitution. Below, I define a subfunction for forward
% substitution which gives us vector c.
c = forwSub(L,b);
% After obtaining vector c, we can use it, and the matrix U, to calculate
% X, as in equation Ux=c. We can do this by back substitution, which is a
% subfunction I defined below.
x = backSub(U,c);
% Thus, we obtain the x vector.
% Subfunctions forwSub and backSub:
    function c = forwSub(L,b)
        % Forward substitution, as in equation Lc=b, solving for c.
        % We are creating the loop that goes through all entries: from each
        % row i from top to bottom, going through all columns j.
        for i = 1:m
            % We need to somehow define the vector c, and there is a few ways to do it.
            % One of them is setting the c(i)=b(i). In this way we first store 
            % c(1)=b(1) as we go through all the columns in the first row.
            % Then in the second, etc.
            c(i) = b(i); % We don't need to divide by L(i,i)=1, because L is the 
            % lower triangle.
            for j=i+1:m
                % We solve for b(j) on the current columns. If we first had
                % c(1)=b(1), then next we are calculating b(2), which we
                % later equal our c(2) to, and so on. This way we can keep
                % track of c. The formula for the forward substitution
                % comes from class notes (also, possible to find by hand).
                b(j)=(b(j)-L(j,i)*c(i))/L(i,i); % We don't really need L(i,i)=1
            end
        end
        % We need a transpose of the vector c to make it a column vector.
        disp('Our vector c is:')
        c=c'
    end
    function x = backSub(U,c)
        % Function for back substitution, as in equation Ux=c, solving for
        % x.
        % We start the loop from backwards, from m, m-1,...,1. Similarly 
        % like in the forward substitution, we set x(j)=c(j)/U(j,j). But, in
        % this case, unlike in forward substitution, there is a division by
        % the jth component of the our solution by the diagonal element of
        % the matrix, U(j,j). The code is more or less the same, we just go
        % backwards.
        for j=m:-1:1
            x(j)=c(j)/U(j,j);
            for k=1:j-1
            c(k)=c(k)-U(k,j)*x(j); % As we can see, we are counting backwards.
            % The formula is also from class notes and hand work, but is
            % analogical to the forward substitution. Only in this case, we
            % use U. So if we were looking at x(4) of a 4x4 matrix in the
            % first round of the loop, we would be calculating x(3) from
            % the next loop.
            end
        end
        disp('Our vector x is:')
        x=x' % Again, we want a column vector so we transpose x.
    end
disp('Thus, we have the solutions of x(1),...,x(m) in the form of the vector x (as in Ax=b):')
x
end
