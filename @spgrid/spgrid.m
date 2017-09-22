classdef spgrid < handle
    %SPGRID A class of sparse grids.
    %  SP = SPGRID( D,L,TYPE ) construct an object of SPGRID.
    %  TYPE must be one of {'disinner','disall','quad','smolyak'}, and 'quad' is
    %  the same as 'smolyak'. The default TYPE is 'disall'. D is the dimension
    %  of variable, L is the level of sparse grids, whicth is interprated as
    %  following:
    %
    %  disall & disinner:  D <= |level|_1 <= L + D.
    %                --->  indset := level
    %                --->  |indset  - 1|_1 < = L
    %                --->  |indset|_inf = L + 1    |indset|_min = 1
    %
    %  quad & Smolyak :    D <= |level|_1 <= L + D
    %                --->  |level-1|_1 < = L
    %                --->  indset := level-1
    %                --->  |indset|_inf = L   |indset|_min = 0;
    %
    %
    %  for Smolyak, the disadjoint points in each level is like this
    %  ------------------------------------------------------------------
    %  level                                    l := level - 1.
    %   2       o                       o       1   include
    %   1                   o                   0
    %   3             o           o             2
    %   4          o     o     o     o          3
    %           o  o  o  o  o  o  o  o  o
    %           1  2  3  4  5  6  7  8  9
    %  -----------------------------------------------------------------
    %
    %  for disinner, the disadjoint points in each level is like this
    %  -----------------------------------------------------------------
    %  level                                      l := level
    %   0       o                       o         0   exclude
    %   1                   o                     1
    %   2             o           o               2
    %   3          o     o     o     o            3
    %           o  o  o  o  o  o  o  o  o
    %           1  2  3  4  5  6  7  8  9
    %------------------------------------------------------------------
    %
    %  NOTE: This is differnt from spgrid.indx. 
    %  
    %  All the defalut values of inputs are setted in parsespin.
    %
    % Example:
    %
    %  sp = spgrid(2, 5, 'disall');
    %  plot(sp, 'r.')
    %
    % See also spgrid.indx, spgrid.parsespin.
    
    % Checked: 14-Sep-2017.
    % $Last revised: 14-Sep-2017$
    % Copyright (c) Guanjie Wang, wangguanjie0@126.com
    
    properties
        d;              % dimension  [done]
        l;              % level [done]
        type;           % type of grids, {'disall','disinner','quad','smolyak'} [done]
        x;              % Nodes [done]
        I;              % index of Nodes for full tensor [done]
        w;              % quadrature weights
        v;              % barycentric weights for interpolation
    end
    
    methods
        
        function sp = spgrid(d,l,varargin)
            % SPGRID Counstruct an obejct of  SPGRID class.
            
            % parse whether all varargin is char.
            [varargin, isallchar] = spgrid.parsespin(varargin,'allchar');
            if ( ~isallchar )
                warning('SPGRID:InPut',...
                    ['Inputs after the second input',...
                    'should be a CHAR or NUMERIC, otherwise, it will be ignored.']);
            end
            
            % Parse the type of grids
            type  = spgrid.parsespin( varargin, 'type');
            valid_type = {'disall', 'disinner', 'quad', 'smolyak'};
            if ( ~ismember(type, valid_type) )
                error('spgrid:InPut',...
                    'TYPE must be one of {''disinner'',''disall'',''quad'',''smolyak''}')
            end
            
            %% compute sp.I and sp.x
            % for special case
            if ( ismember(type, {'quad','smolyak','disinner'}) && (l == 0) )
                sp.I = ones(1,d);
            elseif ismember(type, {'disinner','disall'})
                l = l + 1;
            end
            sp.I = spgrid.constructor( d, l, type );
            x = chebpts(2^l+1); sp.x = x(sp.I);
            sp.d = d; sp.l = l; sp.type = type;
        end
        
        
        
    end
    
    methods ( Access = public, Static = false )
        varargout = plotsp(sp,varargin);
    end
    methods ( Access = public, Static = true )
        
        [type,number] = parsespin(spvar, par)   % Parse input
        indset = starbar(d, l, type )           % indset of level l for d dimension
        indset = inner2all(indset)              % turn inner indset to all indset
        indset = indx(d, l, varargin);          % total indset
        numsp = numsp(indset, type, level)      % number of sparse grids
        [ss,ll] = hier(l, type)                 % index of grids in Hierarchical order.
        [x,y] =  constructor(d, l, type)        % main constructor of spgrid
        [ff, dff] = hierlag(l, type)            % High order Heirarchial basis
        [ff,dff] = heirlinear(l, x, type)       % piece wise linear Hierarchical basis
        
    end
    
end

