% Nguyen Thuy Linh
% MATH 301 MIDTERM
% This is a function that performs LU factorization with pivots.
% We can input some matrix A. Some suggestions of matrix A that could not work
% with the old LU factorization without pivots are given.
% I provided step-by-step explanations for all actions.

function [P, L, U] = LUfactorizationWithPivots(A)
% A = [1 1 1; 2 2 3; 1 2 3];
% A = [10 -7 0; -3 2 6; 5 -1 5];
% A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
[m,n] = size(A);
% We make identity matrices.
P = eye(m);
L = eye(m);
% First, we have to make sure that A is an n x n matrix.
if m ~= n
    disp('Matrix A is not a square matrix')
    return
end
% Now, we perform LU factorization with pivoting, where PA = LU.
for M=1:m
    % We are going through all entries column-wise (e.g. A(1,1), A(2,1)
    % etc. and then the next columns) to find the biggest entry in the
    % matrix A.
    % I found from help max that [Y,I]=max(X) returns the indices of the maximum values
    % in vector I, and we are interested in the I value (the Ith row, where the max is).
    [Y,I] = max(abs(A(M:m,M)));
    % I indicates in which row we have the highest value.
    I = I+M-1;
    % Now we can proceed to switching rows I and M.
    save = A(M,:); % We are saving the whole Mth row in the variable save.
    % Since we performed the act of exchanging rows, we have to record
    % the pivot row swap in our P matrix (we do same operations, but on matrix P).
    savep = P(M,:); % We are saving the whole Mth row in P as well.
    A(M,:) = A(I,:); % We are replacing the Mth row with the Ith row so that the maximum values are in the higher row (it will be ordered this way).
    P(M,:) = P(I,:); % Recording the replacement in P.
    A(I,:) = save; % In the place of Ith row, we copy the original Mth row that we saved in the variable save.
    P(I,:) = savep; % Recording the same thing in P.
    % We swapped rows in A, recorded the act in matrix P. Now, we need
    % to perform LU factorization just like in the submitted homework, where we are 
    % trying to make zeros (first, A(2,1), then A(3,1) and A(3,2)), which means we 
    % have to start from the second row. We make an appropriate
    % condition:
    if M >= 2
        % Before proceeding to LU factorization, we make the same
        % permutations we did on the A matrix to the L matrix.
        savel = L(M,1:M-1); % We save the entries in the Mth row where we need to get zeros (e.g. ig M=3, we get entries (3,1), (3,2) etc.).
        L(M,1:M-1) = L(I,1:M-1); % We swap those entries with respective entries in the row I.
        L(I,1:M-1) = savel; % We copy the saved entries from row M to row I.
    end
    % Now we can perform LU factorization without problems we encountered in the code
    % in the previous homework (because we added pivoting), using the same code with 
    % some adjeustments, since I used 2:m before and here we started 1:m. The code is:
    for N = M+1:m
        E = eye(m); % Creating some identity matrix.
        % Making zeros by dividing pivot columns by pivot entries:
        E(N,M) = -A(N,M)/A(M,M); % Making zeros in entries (2,1), (3,1), (3,2) etc.
        A = E*A; % Recording the U matrix, which is the result of performed operations.
        L = E*L; % Records a matrix whose inverse will give us the L we are looking for.
    end
end
disp('Our P, L, U are as follows:')
P
L = inv(L) % To find L we are taking the inverse of the L we got by performing LU factorization.
U = A
end    