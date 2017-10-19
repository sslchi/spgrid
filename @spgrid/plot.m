function varargout = plot(sp, varargin)
% PLOT Plot the sparse grids for 2D and 3D.

if (~isa(sp,'spgrid'))
    error('spgrid:plot:InPut',...
        'SP must be a SPGRID object for 2D or 3D.\n')
elseif ( (sp.dim > 3) || (sp.dim < 2) )
     error('spgrid:plot:InPut',...
        'SP must be a SPGRID object for 2D  or 3D.\n')
end

if ( nargin == 1 )
    
    if (sp.dim == 2)
        h = plot(sp.x(:,1),sp.x(:,2),'.');
    elseif (sp.dim == 3)
        h = plot3(sp.x(:,1),sp.x(:,2),sp.x(:,3),'.');
    end
    
else
    
    if (sp.dim == 2)
        h = plot(sp.x(:,1),sp.x(:,2),varargin{:});
    elseif (sp.dim == 3)
        h = plot3(sp.x(:,1),sp.x(:,2),sp.x(:,3),varargin{:});
    end
    
end          
        

if ( nargout > 0 )
    varargout = {h};
end

end