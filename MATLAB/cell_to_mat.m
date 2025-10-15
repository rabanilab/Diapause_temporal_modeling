function a = cell_to_mat(c)
% ----------------------------------------------------------
% convert a cell c into a matrix a
% perserving dimensionality of c
% ----------------------------------------------------------

[n,m] = size(c);

for (i=1:n)
    for (j=1:m)
        if (isempty(c{i,j}))
            a(i,j) = NaN;
        else
            a(i,j) = c{i,j};
        end
    end
end
