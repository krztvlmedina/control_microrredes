function padded_mat = pad_mat(mat, n, Ny)
    [~, columns] = size(mat);    
    padded_mat = zeros(Ny, n);
    padded_mat(1:Ny, 1:columns) = mat(1:Ny, :);
end