function isok = isfinitearray(A)
% verify that each element of an array is finite
% 20181209 Kurt Feigl
ncount = numel(A);
nfinite = sum(isfinite(reshape(A,ncount,1)));
if ncount == nfinite
    isok = 1;
else
    isok = 0;
end
return


end

