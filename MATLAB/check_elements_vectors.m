function [ her ] = check_elements_vectors(vec1, vec2 )
%CHECK_ELEMENTS_VECTORS Summary of this function goes here
%   Detailed explanation goes here

count = size(vec1);
her = 0;
for j=1: count(2)
    if vec1(j) < vec2(j)
        her = j
    end
end