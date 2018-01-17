classdef ABHI_Node < handle

% Node class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        node_coord
        node_dof
        node_number    
    end
    
% Property definitions
%
%   node_coord  == a 3x1 vector containing the x, y, and z coordinates of the
%                  node
%   node_dof    == a 6x1 vector of the degrees of freedom numbers associated 
%                  with the node
%   node_number == a scalar corresponding to the node number in the structure
%                  as specified by the user in MASTAN
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = ABHI_Node(node_coord,node_number)
            self.node_coord = node_coord;
            self.node_number = node_number;
            self.node_dof = self.AssignDOF;
        end
        
        %% Get Node Coordinates
        %  Returns "node_coord" for access outside the node class
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
        
        %% Get Node Dofs
        %Return 6 DOFS associated with node for access outside the node
        %class
        function node_dof = GetNodeDof(self)
            node_dof = self.node_dof;
        end      
    end
    
    % Private methods go here
    methods (Access = private)
        %% Assign DOF numbers 
        % Assigns DOF numbers corresponding to every node number
        function node_dof=AssignDOF(self)
            node_dof=6*(self.node_number-1) + (1:6)';
            self.node_dof = node_dof;
        end 
    end
end
