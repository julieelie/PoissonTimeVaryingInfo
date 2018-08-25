function [NewPath, NotLastPath] = updatePath(Path,MaxValues,W)
NotLastPath=1;
if nargin<=2
    W=length(Path);
end
if Path(W)==MaxValues(W) % we reach the end of a cycle
    if (sum(Path(1:W) == MaxValues(1:W)) == W) % Are we at the last possible path?
        NotLastPath=0;
        NewPath = [];
        return
    else
        J=0;
        while true
            if Path(W-J)==MaxValues(W-J)
                J = J+1;
            else
                break
            end
        end
        K = W - J;
    end
else
    K=W;
end

NewPath = Path;
NewPath(K) = NewPath(K) + 1;
if K<length(Path)
    NewPath(K+1 : end) = ones(length(NewPath) - K,1);
end
end