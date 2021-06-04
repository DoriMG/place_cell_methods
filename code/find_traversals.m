function tras = find_traversals(loc, thresh)
if ~exist('thresh','var')
    thresh = 30;
end
  % Find the traversals
    tra = 1;
    tras = [tra];
    for i = 1:length(loc)-1
        if loc(i+1)-loc(i) < -thresh
            tra = tra+1;
        end
        tras = [tras, tra];
    end
    
end