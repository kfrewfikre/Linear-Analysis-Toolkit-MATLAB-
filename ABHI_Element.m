classdef ABHI_Element < handle

% Element class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
       element_nodes
       E
       v
       A
       Ayy
       Azz
       Iyy
       Izz
       J
       Length
       webdir
       Transformation
       k_local
       k_global
       elem_dof
       FEF_local      
    end

%% Property Definitions
%
%   element_nodes  == a 2x1 array of node objects that define the nodes at the
%                     ends of the element
%   Length         == the length of the element
%   Transformation == the 12x12 element transformation matrix
%   k_local        == the 12x12 element stiffness matrix in local coordinates
%   k_global       == the 12x12 element stiffness matrix in global coordinates
%   elem_dof       == a 12x1 vector containing the numbers of the degrees
%                     of freedom for the element in order of delta_x,delta_y,
%                     delta_z, theta_x, theta_y, theta_z for the first end 
%                     followed by delta_x, delta_y, delta_z, theta_x, theta_y,
%                     and theta_z for the second end
%   FEF_local      == a 12x1 vector of the fixed end forces on the element 
%                     in local coordinates where FEF_local(i) = the fixed
%                     end force on the element at degree of freedom i
%                     

    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = ABHI_Element(element_nodes,E,v,A,Ayy,Azz,Iyy,Izz,J,webdir)
        
            self.element_nodes=element_nodes;
            self.E=E;
            self.v=v;
            self.A=A;
            self.Ayy=Ayy;
            self.Azz=Azz;
            self.Iyy=Iyy;
            self.Izz=Izz;
            self.J = J;
            self.webdir=webdir;
            
            %Methods called from the constructor
            self.ComputeLength();
            self.ComputeTransformationMatrix();
            self.k_global = self.ComputeElasticStiffnessMatrix();
            self.elem_dof = self.RetrieveDOF;
        end
        
                
        %% Get global stiffness matrix
        % Function to query global stiffness matrix from outside the element class
        function k_global = GetGlobalStiffnessMatrix(self)
            k_global = self.k_global;
        end
        
        %% Get element Dofs
        %Function to query 12x1 element degree of freedom number vector
        %from outside the element class
        function elem_dof = GetElemDof(self)
            elem_dof = self.elem_dof;
        end
        
        %% Compute Forces
        % Computes the local element forces on the element based on the
        % global deltas and fixed end forces for that element
        function [elem_forces] = ComputeForces(self,delta_glob)
            
            %convert the the global deltas to local deltas by
            %multiplying by the tranformation matrix
            delta_loc = self.Transformation*delta_glob;
            
            %Get local element forces
            elem_forces = self.k_local*delta_loc;
            elem_forces = elem_forces + self.FEF_local;
        end
        
        %% Fixed End Forces
        % Function to calculate Fixed End Forces on element
        % Assigns fixed end forces at local and global degrees of freedom
        function FEF_global=FixedEndForces(self,w) 
            UDL=w;
            I=find(UDL);
            FEF_local=zeros(12,1);
            L=self.Length;
            FEF_local(2)=-UDL(2)*L/2;
            FEF_local(3)=-UDL(3)*L/2;
            FEF_local(5)=UDL(3)*L^2/12;
            FEF_local(6)=-UDL(2)*L^2/12;
            FEF_local(8)=-UDL(2)*L/2;
            FEF_local(9)=-UDL(3)*L/2;
            FEF_local(11)=-UDL(3)*L^2/12;
            FEF_local(12)=UDL(2)*L^2/12;
            self.FEF_local = FEF_local;
            
            T=self.Transformation;
            FEF_global=T'*FEF_local;          
        end
    end
    
        % Private methods go here
        
        methods (Access = private)
        %% Compute Transformation Matrix
        % Compute the element's geometric transformation matrix
        function ComputeTransformationMatrix(self)
                 Transformation=zeros(12,12);      % Initializes the transformation matrix with 
                                                   % zeros 
                 
                 Node1=self.element_nodes(1);      % calls the coordinates of node 'i'
                 Node2=self.element_nodes(2);      % calls the coordinates of node 'j'
                 N1=Node1.GetNodeCoord();          % extracts the coordinates of node 'i'
                 N2=Node2.GetNodeCoord();          % extracts the coordinates of node 'j'
                 N=N2-N1;
                 Dir_X=N/norm(N);                  % finds the direction cosines for member 
                                                   % axis direction
                 Dir_Y=self.webdir;                % extracts the direction cosine for web 
                                                   % direction
                 M=cross(Dir_X,Dir_Y);             % Performs vector cross product
                 Dir_Z=M/norm(M);                  % finds the direction cosine in the Z' axis
                
                 c=1;
                 for i=1:3:12                      % Nested for-loop to generate 12X12 
                                                   % Tranformation Matrix
                  for j=i:i+2                          
                     Transformation(i,j)=Dir_X(c);
                     Transformation(i+1,j)=Dir_Y(c);
                     Transformation(i+2,j)=Dir_Z(c);
                     c=c+1;
                  end
                  c=1;
                 end
                 self.Transformation = Transformation;

        end
        
        %% Retrieve DOF
        % Function to assign degrees of freedom to element
        function element_dof=RetrieveDOF(self)
            node1=self.element_nodes(1);  
            node2=self.element_nodes(2);
            element_dof=[node1.GetNodeDof();node2.GetNodeDof()];%assign degree of freedom numbers                                                                           
            self.elem_dof = element_dof;                                       
        end
        
    
        %% Compute Length
        % Compute the element's length
        function ComputeLength(self)
            
            node1 = self.element_nodes(1);         %get node objects
            node2 = self.element_nodes(2);
            
            node1coord = node1.GetNodeCoord();	   %get node coordinate vectors
            node2coord = node2.GetNodeCoord();
            
            d = node2coord - node1coord;           %find euclidean distance between nodes
            Length = norm(d); 
            self.Length = Length;
            
        end
        
        %% Compute Elastic Stiffness Matrix
        % Compute the element's elastic stiffness matrix in local and global coordinates
        function k_global = ComputeElasticStiffnessMatrix(self)
            %Define method variables
            E = self.E;
            A = self.A;
            v = self.v;
            Iyy = self.Iyy;
            Izz = self.Izz;
            J = self.J;
            k = zeros(12);
            L = self.Length;
            Ayy = self.Ayy;
            Azz = self.Azz;
            G = E/(2*(1+v));
            phiy = (12*E*Izz)/(G*Ayy*L^2);
            phiz = (12*E*Iyy)/(G*Azz*L^2);
            T = self.Transformation;
            
            %Stiffness Coefficients for typical 12 DOF frame element
            %including shear
            S1 = (E*A)/L;
            S2 = (12*E*Izz)/((L^3)*(1+phiy));
            S3 = (6*E*Izz)/(L^2*(1+phiy));
            S4 = (12*E*Iyy)/(L^3*(1+phiz));
            S5 = (6*E*Iyy)/(L^2*(1+phiz));
            S6 = ((4+phiy)*E*Izz)/(L*(1+phiy));
            S7 = ((4+phiz)*E*Iyy)/(L*(1+phiz));
            S8 = ((2-phiz)*E*Iyy)/(L*(1+phiz));
            S9 = ((2-phiy)*E*Izz)/(L*(1+phiy));
            S10 = ((G*J)/L);
            
            %Assemble top half of stiffness matrix
            k(1,1) = S1;
            k(1,7) = -S1;
            k(2,2) = S2;
            k(2,6) = S3;
            k(2,8) = -S2;
            k(2,12) = S3;
            k(3,3) =  S4;
            k(3,5) = -S5;
            k(3,9) = -S4;
            k(3,11) = -S5;
            k(4,4) = S10;
            k(4,10) = -S10;
            k(5,5) = S7;
            k(5,9) = S5;
            k(5,11) = S8;
            k(6,6) = S6;
            k(6,8) = -S3;
            k(6,12) = S9;
            k(7,7) = S1;
            k(8,8) = S2;
            k(8,12) = -S3;
            k(9,9) = S4;
            k(9,11) = S5;
            k(10,10) = S10;
            k(11,11) = S7;
            k(12,12) = S6;
            
            %fill out bottom half of stiffness matrix
            
            %create matrix with just diagonal terms
            kD = k.*eye(12); 
            
            % create upper triangular
            kU = k - kD; 
            
            % create lower triangular 
            kL= kU'; 
            
            % sum upper, diagonal, and lower for full stiffness matrix
            k = kL + kD + kU;
            k_local = k;
            
            %Find global stiffness matrix
            k_global = T'*k_local*T;
            
            self.k_local = k_local;
            self.k_global = k_global;
        end      
    end
end
