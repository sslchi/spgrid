function [ff,dff] = heirlinear(l,x,type)



if nargin < 3
    type = 'disall';
end

switch type
    case {'disall'}
        
        
end


end



function [ff,dff] = mother_linear(l,x)


if ( (x > 1) || (x < -1) )
    ff = 0; dff = 0;
else
    ff = 1 - abs(x);
    if ( x > 0 )
        dff = 1;
    else 
        dff = -1;
    end
end




end