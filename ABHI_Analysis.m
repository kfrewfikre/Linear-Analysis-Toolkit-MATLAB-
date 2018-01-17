classdef ABHI_Analysis < handle
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Analysis class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        nnodes
        nele
        coord
        ends
        w
        concen
        webdir
        E
        v
        A
        Ayy
        Azz
        Iyy
        Izz
        J
        fixity
       
        eleobj
        nodeobj
        free_dof
        supp_dof
        disp_dof
        K
        Kff
        Kfn
        Ksf
        Ksn
        Knf
        Knn
        Pff
        delf
        deln
        DEFL
        REACT
        ELE_FOR
        FEF
    end
    
%% Property Definitions
% Note: variables defined in ud_3d1el.m are not repeated here
%
%   eleobj   == a nele x 1 array of element objects where eleobj(i)
%               is an element object corresponding the ith element
%   nodeobj  == a nnodes x 1 array of node objects where nodeobj(i) is the 
%               node object corresponding to the ith node
%   free_dof == a column vector that stores the numbers of the free degrees
%               of freedom in numerical order
%   supp_dof == a column vector that stores the numbers of the supported 
%               degrees of freedom in numerical order
%   disp_dof == a column vector that stores the numbers of the displaced
%               degrees of freedom in numerical order
%   K        == the ndof x ndof global stiffness matrix for the structure
%   Kff      == the structure stiffness submatrix that has rows corresponding to 
%               the free degrees of freedom and columns corresponding to the 
%               free degrees of freedom
%   Kfn      == the structure stiffness submatrix that has rows corresponding 
%               to the free degrees of freedom and columns corresponding to the 
%               displaced degrees of freedom
%   Ksf      == the structure stiffness submatrix that has rows corresponding 
%               to the supported degrees of freedom and columns corresponding 
%               to the free degrees of freedom
%   Ksn      == the structure stiffness submatrix that has rows corresponding to the 
%               supported degrees of freedom and columns corresponding to
%               the displaced degrees of freedom
%   Knf      == the structure stiffness submatrix that has rows corresponding to the
%               displaced degrees of freedom and columns corresponding to
%               the free degrees of freedom
%   Knn      == the structure stiffness submatrix that has rows
%               corresponding to the displaced degrees of freedom and
%               columns corresponding to the displaced degrees of freedom
%   Pff      == a column vector that contains the concentrated loads minus
%               the fixed end forces that are applied to the free degrees of 
%               freedom of the structure listed in numerical order
%   delf     == a column vector that contains the displacements calculated
%               at the free degrees of freedom listed in numerical order
%   deln     == a column vector that contains the displacements at the 
%               displacement specified supports in numerical order
     
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = ABHI_Analysis(nnodes,nele,coord,ends,w,concen,webdir,...
                E,v,A,Ayy,Azz,Iyy,Izz,J,fixity)
            
            self.nnodes=nnodes;
            self.nele=nele;
            self.coord=coord;
            self.ends=ends;
            self.w=w;
            self.concen=concen;
            self.webdir=webdir;
            self.E=E;
            self.v=v;
            self.A=A;
            self.Ayy=Ayy;
            self.Azz=Azz;
            self.Iyy=Iyy;
            self.Izz=Izz;
            self.J=J;
            self.fixity = fixity;
            
        end
        
        %% Run analysis
        % Performs the series of function calls in the Analysis, Element,
        % and Node classes to perform the analysis
        function RunAnalysis(self)
            [self.free_dof, self.supp_dof, self.disp_dof] = self.classifyDOF();
            self.nodeobj=self.CreateNodes();
            self.eleobj=self.CreateElement();
            K=self.CreateStiffnessMatrix();
            self.CreateLoadVectors();
            self.K = K;
            [Kfn, Knn, Ksn, Kff, Knf, Ksf] = self.CreateStiffnessSubMatrices();
            [DEFL, REACT] = self.ComputeDisplacementsReactions();
            [ELE_FOR] = self.RecoverElementForces();                   
        end
        
        %% Get Mastan2 Returns
        % Function to return the necessary outputs to Mastan2 to be used by
        % the Mastan postprocessor
        function [DEFL, REACT, ELE_FOR, AFLAG] = GetMastan2Returns(self)
           DEFL = self.DEFL;
           REACT = self.REACT;
           ELE_FOR = self.ELE_FOR;
           AFLAG = self.CheckKffMatrix();
           self.ComputeError()
       end
    end     
    
    % Private methods go here
    methods (Access = private)
        
        %% Create Nodes Function
        % Creates an array of node objects for all of the coordinates given
        % by Mastan in the matrix named "coord"
        function nodes=CreateNodes(self)
            nodes=ABHI_Node(self.coord(1,:),1);
            for i=2:self.nnodes
                nodes=[nodes;ABHI_Node(self.coord(i,:),i)];
            end
        end
        
        %% Create Element Function
        % Creates an array of element objects for all the elements given
        % provided by mastan 
        function element=CreateElement(self)
            %Creates first entry in the element object array
            elem_ends=self.ends(1,:);
            node_ele=[self.nodeobj(elem_ends(1));self.nodeobj(elem_ends(2))];
            element=ABHI_Element(node_ele,self.E(1),self.v(1),self.A(1),...
                self.Ayy(1),self.Azz(1),self.Iyy(1),self.Izz(1),self.J(1)...
                ,self.webdir(1,:)');
            
            %Generates entries 2 through nele in the element object array
            for i=2:self.nele
                elem_ends=self.ends(i,:);
                node_ele=[self.nodeobj(elem_ends(1));self.nodeobj(elem_ends(2))];
                element=[element;ABHI_Element(node_ele,self.E(i),self.v(i)...
                    ,self.A(i),self.Ayy(i),self.Azz(i),self.Iyy(i),...
                    self.Izz(i),self.J(i),self.webdir(i,:)')];
            end
        end
        
        %% Create Stiffness Matrix Function
        % Create the entire global stiffness matrix for the structure
        function K=CreateStiffnessMatrix(self)
             
            nDOF=self.nnodes*6;
            %initialize a matrix with columns and rows equal to nDOF
            K=zeros(nDOF,nDOF);
            
            %loop through all elements
            for k=1:self.nele
                %Query the kth elements k_global matrix
                Kele=self.eleobj(k).GetGlobalStiffnessMatrix();
                %Query degrees of freedom associated with element
                element_dof=self.eleobj(k).GetElemDof();
               
                %Set the structure global stiffness matrix at the locations
                %corresponding to the element degrees of freedom
                K(element_dof, element_dof) = K(element_dof,element_dof) + Kele;
            end
            K = sparse(K);
        end
         
        %% Create Load Vectors Function
        % Create Fixed End Force and Concentrated Force Vectors in Global
        % Coordinates. Also stores Pff which is the concentrated forces
        % minus the fixed end forces acting on the free degrees of freedom
        function CreateLoadVectors(self)
             nDOF=6*self.nnodes;
             FEF=zeros(nDOF,1);
             FixedEndElements = zeros(12,self.nele);
             for i = 1:self.nele
                 FixedEndElements(:,i)=self.eleobj(i).FixedEndForces(self.w(i,:));
             end
                %loop through all elements
                  for j=1:self.nele
                      fef_elem_global=FixedEndElements(:,j); %Col Vector
                      elem_dof=self.eleobj(j).GetElemDof();
                      FEF(elem_dof) = FEF(elem_dof) + fef_elem_global;
                  end
            self.FEF = FEF;
            %Using linear indexing, create a nDOF x 1 vector from nnodes
            %matrix called 'concen' that lists the concentrated loads on
            %all degrees of freedom
            concen_trans=(self.concen)';
            concen_loads=zeros(nDOF,1);
            for i=1:nDOF
                concen_loads(i)=concen_trans(i);
            end
            concen_trans=self.concen';
            concen_loads=zeros(nDOF,1);
            for i=1:nDOF
                concen_loads(i)=concen_trans(i);
            end
            %Subtract fixed end forces from concentrated loads vector
            Load=concen_loads-FEF;
            %Create load vector for free degrees of freedom
            L=length(self.free_dof);
            Pff=zeros(L,1);
            for i=1:L;
                Pff(i)=Load(self.free_dof(i));
            end
            self.Pff=Pff;
        end
        
        
            
           %% Classify DOF function
           % create vectors containing the numbers of the free, displaced,
           % and support degrees of freedom
           
            function [free_dof, supp_dof, disp_dof] = classifyDOF(self)
                fixity = self.fixity;
                fix_t = fixity';
                free_dof = find(isnan(fix_t));
                supp_dof = find(fix_t == 0);
                disp_dof = find(fix_t~=0 & ~isnan(fix_t));
                self.free_dof = free_dof;
                self.supp_dof = supp_dof;
                self.disp_dof = disp_dof;
            end
            %% Create Stiffness submatrices function
            % creates Kfn, Knn, Ksn, Kff, Knf, and Ksf submatrices based
            % from the global K matrix
            function [Kfn, Knn, Ksn,Kff,Knf,Ksf] = CreateStiffnessSubMatrices(self)
              [free_dof,supp_dof,disp_dof]=self.classifyDOF();
              L1=length(free_dof);
              L2=length(supp_dof);
              L3=length(disp_dof);
              Kfn=zeros(L1,L3);
              Knn=zeros(L3,L3);
              Ksn=zeros(L2,L3);
              Kff=zeros(L1,L1);
              Knf=zeros(L3,L1);
              Ksf=zeros(L2,L1);
              K = self.K;
              
              %Creates submatrices using the classified degree of freedom
              %vectors as indices
              Kfn = K(free_dof,disp_dof);
              Knn = K(disp_dof, disp_dof);
              Ksn = K(supp_dof, disp_dof);
              Kff = K(free_dof,free_dof);
              Knf = K(disp_dof, free_dof);
              Ksf = K(supp_dof, free_dof);
                
              %Set property variables
              self.Kff=Kff;
              self.Kfn=Kfn;
              self.Ksf=Ksf;
              self.Ksn=Ksn;
              self.Knf=Knf;
              self.Knn=Knn;
                
            end
            %% Check Kff Matrix Function
            % This function checks the condition number of the Kff matrix
            % and reports the number of significant digits lost in the
            % analysis
            function flag=CheckKffMatrix(self)
                k=condest(self.Kff);
                p=16;
                s=p-log10(k);
             
                if s<3
                    flag=0;
                else
                flag=1;
                end
                fprintf('The condition number of the Kff matrix is %f\n',k)
                fprintf('The number of significant digits lost in the analysis is %f\n',log10(k)) 
            end
            %% Compute Displacements Reactions Function
            % Computes the global displacements at free degrees of freedom
            % and formats all displacements for output to MASTAN.
            % coputes the global forces at supported and displaced degrees
            % of freedom and formats all the reactions for output to MASTAN
            
            function [DEFL,REACT]=ComputeDisplacementsReactions(self)
                fix_t=self.fixity';
                Ldisp=length(self.disp_dof);
                deln=zeros(Ldisp,1);
                %Create vector of specified displacements deln
                for i=1:Ldisp
                    deln(i)=fix_t(self.disp_dof(i));
                end
                %Compute displacements at free dofs and reactions at supports
                if ~isempty(deln) || ~isempty(self.Kfn)
                    delf=self.Kff\(self.Pff - (self.Kfn*deln));
                    Ps=(self.Ksf*delf)+ self.Ksn*deln + self.FEF(self.supp_dof);
                else
                    delf=self.Kff\self.Pff;
                    Ps=(self.Ksf*delf) + self.FEF(self.supp_dof);
                end
                self.delf = delf;
                self.deln = deln;
                %Compute forces at displaced dofs Pn
                if ~isempty(deln) || ~isempty(self.Knn)
                    Pn = self.Knf*delf+self.Knn*deln + self.FEF(self.disp_dof);
                end
                %Format Reaction Matrix for Mastan
                REACT_t = zeros(6,self.nnodes);
                REACT_t(self.supp_dof) = Ps;
                %Incorporate forces at displaced dofs in REACT matrix
                if ~isempty(deln) || ~isempty(self.Knn)
                    REACT_t(self.disp_dof) = Pn;
                end
                REACT = REACT_t';
                self.REACT = REACT;
                DEFL_t = zeros(6,self.nnodes);
                DEFL_t(self.free_dof) = self.delf;
               
                %check to see if there are any prescribed deflections and add
                %them in to the DEFL matrix if so
                if ~isempty(self.disp_dof)
                    DEFL_t(self.disp_dof) = self.deln ;
                end
                DEFL = DEFL_t';
                self.DEFL = DEFL;
            end
            
           %% Function Recover Element Forces
           % Iterates trough each element of the structure, passes down
           % deflections and returns the element forces in local
           % coordinates. Called from the runanalysis method
           function [ELE_FOR] = RecoverElementForces(self)
                DEFL = self.DEFL;
                DEFL_t = DEFL';
                %Initialize element forces vector
                ELE_FOR = zeros(self.nele,12);
                defl = zeros(12,1);
                for i = 1:self.nele
                    defl(:,1) = DEFL_t(self.eleobj(i).GetElemDof());
                    ELE_FOR(i,:) = self.eleobj(i).ComputeForces(defl);
                    self.ELE_FOR = ELE_FOR;
                end
           end
           
           %% Function Compute error
           % Computes the difference between the applied loads Pf and the
           % back calculated Pf loads
            function ComputeError(self)
                Error = self.Pff - (self.Kff*self.delf + self.Kfn*self.deln);
                disp('Error is')
                disp(Error)
            end         
    end
end

