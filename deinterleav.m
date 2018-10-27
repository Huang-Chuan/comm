%------------------------------------------------------------------

function output1 = deinterleav( sig, Nrows, Ncols)

	i = repmat([0 : Ncols - 1]', Nrows,1);
    j = reshape(repmat([1:Nrows],Ncols,1),Nrows*Ncols,1);

    output1 = sig(Nrows * i + j);
end