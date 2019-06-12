function isok = isfinitearray(A)
% verify that each element of an array is finite
% 20181209 Kurt Feigl
ncount = numel(A);
nfinite = sum(isfinite(reshape(A,ncount,1)));
if ncount == nfinite
    isok = true;
else
    isok = false;
end
return


end

