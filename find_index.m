function idx = find_index(vector, a)

% calculates absolute difference between every element of the vector and a
diff = abs(vector - a);

% then finds the index of the element in vector with the minimal difference
[~, idx] = min(diff);

end