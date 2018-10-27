%------------------------------------------------------------------

function output = interleav( sig, Nrows, Ncols) 
	i = repmat([0 : Nrows - 1]', Ncols,1);
    j = reshape(repmat([1:Ncols],Nrows,1),Nrows*Ncols,1);

    output = sig(Ncols * i + j);
end