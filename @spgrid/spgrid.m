classdef spgrid
    %SPGRID  A class of sparse grids.
    %  SP = SPGRID(D, L, TYPE) construct an object of SPGRID.
    %  D is the dimension, L is the level and Type is the type of the sparse
    %  grids.
    %
    %
    %
    %
    
    %% For disinner & disall:   
    %                      D <= |levelset|_1 <= L + D
    %                --->  |levelset  - 1|_1 < = L    }
    %                      |levelset|_min = 1         }
    %                --->  |levelset|_inf = L + 1    
    %
    % the disadjoint points in each level is like this:
    %                                                  
    %                                           _____levelset_____              
    %   L                                       disinner    disall 
    %           o                       o                      0              
    %   0                   o                       1          1  
    %   1             o           o                 2          2  
    %   2          o     o     o     o              3          3  
    %           o  o  o  o  o  o  o  o  o
    %           1  6  4  7  3  8  5  9  2       Nodes in Hierarchical order
    %           1  2  3  4  5  6  7  8  9       Nodes in oordinary order
    %
    %% For quad & Smolyak : 
    %                      D <= |levelset|_1 <= L + D
    %                --->  |levelset - 1|_1 < = L   }
    %                      |levelset - 1|_min = 0   }
    %                --->  |levelset - 1|_inf = L   
    %
    % the disadjoint points in each level is like this:
    %
    %   L                                       levelset                 
    %
    %   1       o                       o          2               
    %   0                   o                      1               
    %   2             o           o                3               
    %   3          o     o     o     o             4              
    %           o  o  o  o  o  o  o  o  o
    %           2  6  4  7  1  8  5  9  3       Nodes in Hierarchical order
    %           1  2  3  4  5  6  7  8  9       Nodes in oordinary order

    
    
    properties
        dim;            % dimension
        level;          % level 
        type;           % type of grids
        x;              % Nodes in Hierarchical order
        I;              % index of Nodes
        w;              % quadrature weights
        v;              % barycentric weights for interpolation
    end
    
    methods
        function sp = spgrid(dim, level, type, varargin)
            
            % Parse input
            valid_type = {'disall', 'disinner', 'quad', 'smolyak'};
            if ( ~ismember(type, valid_type) )
                error('spgrid:InPut',...
                    'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
            end
            
            if  ( ismember(type, {'disinner','disall'}) )
                [x, w, v] = chebpts(2^(level+1) + 1); 
            elseif ( ismember(type, {'quad', 'smolyak'}) )
               if ( level == 0 ) 
                   [x, w, v] = chebpts(1); 
               else
                   [x, w, v] = chebpts(2^level + 1); 
               end
                
            end
            hid = spgrid.hier(level, type);
            x = x(hid); w = w(hid); v = v(hid);
            
            I = spgrid.constructor(dim, level, type);
            x = x(I);
            
            sp.dim = dim; 
            sp.level = level;
            sp.type = type;
            sp.x  = x;
            sp.I = I;
            sp.w = w;
            sp.v = v;
        end
        
        
    end
    
    methods ( Access = public, Static = false )
        varargout = plot(sp,varargin);
    end
    methods ( Access = public, Static = true )
        
        [type,number] = parsein(spvar, par)    % Parse input
        indset = starbar(d, l, type )           % indset of level l for d dimension
        indset = inner2all(indset)              % turn inner indset to all indset
        indset = levelset(d, l, varargin);      % level set
        numsp = numsp(indset, type, level)      % number of sparse grids
        [ss,ll] = hier(l, type)                 % index of grids in Hierarchical order.
        [x,y] =  constructor(d, l, type)        % main constructor of spgrid
        [ff, dff] = hierlag(l, type)            % High order Heirarchial basis
        [ff,dff] = heirlinear(l, x, type)       % piece wise linear Hierarchical basis
        
    end
    
    
    
end

