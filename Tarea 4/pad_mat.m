function padded_mat = pad_mat(mat, n, Ny)
    [~, columns] = size(mat);
    if (n - columns) >= 0
        max_col = columns;
    else
        max_col = n;
    end

    padded_mat = zeros(Ny, n);
    padded_mat(1:Ny, 1:max_col) = mat(1:Ny, 1:max_col);
end