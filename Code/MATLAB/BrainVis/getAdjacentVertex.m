clear all
addpath('./NII')

% Load atlas
mask = load_nii('desikan.nii');
mask.img = double(mask.img);
data = mask.img;

[nx, ny, nz] = size(data);

d = [0 0 1;
    0 1 0;
    0 1 1;
    1 0 0;
    1 0 1;
    1 1 0;
    1 1 1];

A = zeros(70, 70);

for i = 1:nx
    i
    for j = 1:ny
        for k = 1:nz
            if (data(i, j, k) ~= 0)
                for iD = 1:7
                    if (data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3)) ~= 0)
                        if (data(i, j, k) ~= data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3))) && ...
                                (((data(i, j, k) <= 35) && (data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3)) <= 35)) || ...
                                ((data(i, j, k) > 35) && (data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3)) > 35)))
                            A(data(i, j, k), data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3))) = 1;
                            A(data(i + d(iD, 1), j + d(iD, 2), k + d(iD, 3)), data(i, j, k)) = 1;
                        end
                    end
                end
            end
        end
    end
end

csvwrite('../../../Data/adjacent_list.csv', A)