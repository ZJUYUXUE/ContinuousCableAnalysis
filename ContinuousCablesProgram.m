%% pre-stress adjustment of cable-strut structures using noisy data
function ContinuousCablesProgram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic information:
    % Function: Compute the load responses of continuous cable systems
    % Author: Xue yu
    % Institution: Space Structure Research Centre, Zhejiang University
    % Email: 11812048@zju.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input parameters
FileName = 'Levy1';
% Example : Levy dome with continuous loop cables
NodeTable = xlsread(FileName,'Node');
ElementTable = xlsread(FileName,'Element');
LoadTable = xlsread(FileName,'Load');
ContinuousCableSequence = xlsread(FileName,'Continuous');
% The structural information is input by the excel file.
% Sheet "Node":
%   column 1 :      node numbering
%   columns 2-4:    nodal coordinates in x-, y-, and z- directions(Unit: m)
%   column 5:       "1" means unfixed node, and "0" means fixed node
% Sheet "Element":
%   column 1:       segment numnering
%   column 2£º      segment type£»"0" means cable and "1" means strut
%   column 3£º      segment axial stiffness (Unit: kPa)
%   columns 4-5£º   connecting nodes
%   column 6:       segment initial forces(Unit:kN)
%   column 7:       segment lengths;
% Sheet "Load":
%   column 1:       external load on corresponding DOF(Unit:kN)
% Sheet "Continuous":
%   column 1-2:     continuous segments
AnalysisType = 0;
% 1:definite load path
% 0:indefinite load path

ElementNumber = size(ElementTable,1);
% Number of elements
NodeNumber = size(NodeTable,1);
% Number of Nodes
UnfixedNodeNumber = sum(NodeTable(:,5));
OutputResult(ElementTable,NodeTable);
% Check the initial geometry
ElementStiffness = ElementTable(:,3);
% Segment stiffness vector

C = zeros(ElementNumber,NodeNumber);
for i = 1 : 1 :ElementNumber
        C(i,ElementTable(i,4)) = 1;
        C(i,ElementTable(i,5)) = -1;
end
% Node connectivity matrix

L = zeros(ElementNumber,ElementNumber);
dX = zeros(ElementNumber,ElementNumber);
dY = zeros(ElementNumber,ElementNumber);
dZ = zeros(ElementNumber,ElementNumber);
for i = 1 :1 : ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    L(i,i) = LEN(NodeTable(n,2),NodeTable(n,3),NodeTable(n,4),NodeTable(m,2),NodeTable(m,3),NodeTable(m,4));
    dX(i,i) = NodeTable(n,2)-NodeTable(m,2);
    dY(i,i) = NodeTable(n,3)-NodeTable(m,3);
    dZ(i,i) = NodeTable(n,4)-NodeTable(m,4);
end
D0 = [C'* dX/L;
    C'* dY /L;
    C'* dZ /L];
ZERO = diag([NodeTable(:,5);
    NodeTable(:,5);
    NodeTable(:,5)]);
ZERO (all(ZERO == 0, 2),:) = [];
D = ZERO * D0;
% Equilibrium matrix

K_e = diag(diag(ElementStiffness)/L);
% Segmenet stiffness vector

SegmentPrestress = ElementTable(:,6);
InitialElongations = diag(K_e)^(-1) * SegmentPrestress;
FrictionCoefficient = 0.1;

SlidingLimit = GetSlidingLimit(NodeTable,ElementTable,ContinuousCableSequence,FrictionCoefficient);
% Computing the sliding limits of the continuous segment force ratios

%% computing structural responses when the load is continuously changed
if (AnalysisType == 1)
    ContinuousMatrix = eye(size(K_e,1));
    ForceRelationMatrix = eye(size(K_e,1));
    % Initialize the segment connectivity matrix and the segment force relation matrix
    SegmentForce(:,1) = SegmentPrestress;
    % Initalize the segment forces
    for k = 1 : 1 : size(ContinuousCableSequence,1) + 2
        [ContinuousMatrix,ForceRelationMatrix] = GetMaxtrix(SegmentPrestress,SegmentForce(:,k),ContinuousCableSequence,SlidingLimit,ContinuousMatrix,ForceRelationMatrix);
        % Update the segment connectivity matrix and the segment force relation matrix
        K_E = D * ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * D';
        % Linear stiffness matrix
        RedistributionForce = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * InitialElongations;
        % Redistributed forces caused by the change of nodal slide state
        Q = L^(-1) * diag(RedistributionForce);
        G = C'* Q * C;
        K_G0 = zeros(3 * NodeNumber,3 * NodeNumber);
        K_G0(1:NodeNumber,1:NodeNumber) = G;
        K_G0(NodeNumber + 1:2 * NodeNumber,NodeNumber + 1:2 * NodeNumber) = G;
        K_G0(2 * NodeNumber + 1:3 * NodeNumber,2 * NodeNumber + 1:3 * NodeNumber) = G;
        K_G = ZERO * K_G0* ZERO';
        % Geometrical stiffness matrix
        K_F = D * (RedistributionForce - SegmentPrestress);
        % Unbalanced force caused by the redistribution of the intial forces
        NodalDisplacment(:,k) = (K_E+K_G)^(-1) * (LoadTable - K_F);
        % Nodal displacment solved in current iteration
        SegmentForce(:,k+1) = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * (D' * NodalDisplacment(:,k) + InitialElongations);
        % Segment force solved in current iteration
        if (~IsSlide(SegmentForce(:,k+1),ContinuousCableSequence,SlidingLimit))
            % judge whether the segment force ratio exceeds the sliding limits
            Displacement = NodalDisplacment(:,k);
            FinalForce = SegmentForce(:,k+1);
            break;
        end
    end
PlotResponse(NodalDisplacment(:,k),SegmentForce(:,k+1),ElementNumber,UnfixedNodeNumber);
end
%% Indifinte load path
if (AnalysisType == 0)
    ContinuousCableNumber = size(ContinuousCableSequence,1);
    MaxmimalForce = nan(ElementNumber,1);
    MinimalForce = nan(ElementNumber,1);
    MaximalDeformation = nan(3 * UnfixedNodeNumber,1);
    MinimalDeformation = nan(3 * UnfixedNodeNumber,1);
    for m = 1 : 1 : 2^(ContinuousCableNumber)
        BinIndex = dec2bin(m,ContinuousCableNumber);
        ContinuousMatrix = eye(size(K_e,1));
        ForceRelationMatrix = eye(size(K_e,1));
        [ContinuousMatrix,ForceRelationMatrix] = GetUncertainMaxtrix(ContinuousCableSequence,SlidingLimit,ContinuousMatrix,ForceRelationMatrix,BinIndex);
        % change segment connectivity matrix and segment force relation matrix by binary encoding
        K_E = D * ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * D';
        RedistributionForce = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * InitialElongations;
        Q = L^(-1) * diag(RedistributionForce);
        G = C'* Q * C;
        K_G0 = zeros(3 * NodeNumber,3 * NodeNumber);
        K_G0(1:NodeNumber,1:NodeNumber) = G;
        K_G0(NodeNumber + 1:2 * NodeNumber,NodeNumber + 1:2 * NodeNumber) = G;
        K_G0(2 * NodeNumber + 1:3 * NodeNumber,2 * NodeNumber + 1:3 * NodeNumber) = G;
        K_G = ZERO * K_G0* ZERO';
        K_F = D * (RedistributionForce - SegmentPrestress);
        NodalDisplacment = (K_E+K_G)^(-1) * (LoadTable - K_F);
        FinalForce = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * (D' * NodalDisplacment + InitialElongations);
        MaxmimalForce = max(MaxmimalForce,FinalForce);
        MinimalForce = min(MinimalForce,FinalForce);
        MaximalDeformation = max(MaximalDeformation,NodalDisplacment);
        MinimalDeformation = min(MinimalDeformation,NodalDisplacment);
    end
    PlotResponseRange(MaxmimalForce,MinimalForce,MaximalDeformation,MinimalDeformation,ElementNumber,UnfixedNodeNumber);
end
end

function SlidingLimit = GetSlidingLimit(NodeTable,ElementTable,ContinuousCableSequence,FrictionCoefficient)
ContinuousCableNumber = size(ContinuousCableSequence,1);
SlidingLimit = zeros(ContinuousCableNumber,1);
for i = 1 : ContinuousCableNumber
    Element1 = ContinuousCableSequence(i,1);
    Element2 = ContinuousCableSequence(i,2);
    Node2 = intersect(ElementTable(Element1,4:5),ElementTable(Element2,4:5));
    Node1 = setdiff(ElementTable(Element1,4:5),Node2);
    Node3 = setdiff(ElementTable(Element2,4:5),Node2);
    Angle = SolveThreeNodeAngle(NodeTable,Node1,Node2,Node3);
    SlidingLimit(i) = exp(FrictionCoefficient*Angle);
end
end

function [ContinuousMatrix,ForceRelationMatrix] = GetMaxtrix(InitialSegmentForce,SegmentForce,ContinuousCableSequence,SlidingLimit,ContinuousMatrix,ForceRelationMatrix)
ContinuousCableNumber = size(ContinuousCableSequence,1);
ForceRatio = zeros(ContinuousCableNumber,1);
SlideTime = ones(ContinuousCableNumber,1);
DeltaSegmentForce = SegmentForce - InitialSegmentForce;
for i = 1 : ContinuousCableNumber
    Element1 = ContinuousCableSequence(i,1);
    Element2 = ContinuousCableSequence(i,2);
    ForceRatio(i) = SegmentForce(Element1)/SegmentForce(Element2);
    if(ForceRatio(i) > SlidingLimit(i))
    SlideTime(i) = (InitialSegmentForce(Element1) - SlidingLimit(i) * InitialSegmentForce(Element2))/(SlidingLimit(i) * DeltaSegmentForce(Element2) - DeltaSegmentForce(Element1));
    end
    if(ForceRatio(i) < 1/SlidingLimit(i))
    SlideTime(i) = (InitialSegmentForce(Element1) - 1/SlidingLimit(i) * InitialSegmentForce(Element2))/(1/SlidingLimit(i) * DeltaSegmentForce(Element2) - DeltaSegmentForce(Element1));
    end
end
if(min(SlideTime)<0.999)
    [~,SlidingElement] = min(SlideTime);
    Element1 = ContinuousCableSequence(SlidingElement,1);
    Element2 = ContinuousCableSequence(SlidingElement,2);
if(ForceRatio(SlidingElement) > 1)
    RowElement1 = find(ContinuousMatrix(:,Element1));
    RowElement2 = find(ContinuousMatrix(:,Element2));
    ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
    ContinuousMatrix(Element2,:) = 0;
    ForceRelationMatrix(RowElement1,:) =  ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:)/SlidingLimit(SlidingElement);
    ForceRelationMatrix(RowElement2,:) = 0;
end
if(ForceRatio(SlidingElement) < 1)
    RowElement1 = find(ContinuousMatrix(:,Element1));
    RowElement2 = find(ContinuousMatrix(:,Element2));
    ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
    ContinuousMatrix(Element2,:) = 0;
    ForceRelationMatrix(RowElement1,:) = ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:) * SlidingLimit(SlidingElement);
    ForceRelationMatrix(RowElement2,:) = 0;
end
end
end

function result = IsSlide(SegmentForce,ContinuousCableSequence,SlidingLimit)
result = 0;
ContinuousCableNumber = size(ContinuousCableSequence,1);
ForceRatio = zeros(ContinuousCableNumber,1);
for i = 1 : ContinuousCableNumber
    Element1 = ContinuousCableSequence(i,1);
    Element2 = ContinuousCableSequence(i,2);
    ForceRatio(i) = SegmentForce(Element1)/SegmentForce(Element2);
    if(ForceRatio(i) > SlidingLimit(i) || ForceRatio(i) < 1/SlidingLimit(i))
        result = 1;
    end
end
end

function [ContinuousMatrix,ForceRelationMatrix] = GetUncertainMaxtrix(ContinuousCableSequence,SlidingLimit,ContinuousMatrix,ForceRelationMatrix,BinIndex)
ContinuousCableNumber = size(ContinuousCableSequence,1);
for i = 1 : ContinuousCableNumber
    Element1 = ContinuousCableSequence(i,1);
    Element2 = ContinuousCableSequence(i,2);
    Index(i) = str2num(BinIndex(i));
    if(Index(i) == 1)
        RowElement1 = find(ContinuousMatrix(:,Element1));
        RowElement2 = find(ContinuousMatrix(:,Element2));
        ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
        ContinuousMatrix(Element2,:) = 0;
        ForceRelationMatrix(RowElement1,:) =  ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:)/SlidingLimit(i);
        ForceRelationMatrix(Element2,:) = 0;
    end
    if(Index(i) == 0)
        RowElement1 = find(ContinuousMatrix(:,Element1));
        RowElement2 = find(ContinuousMatrix(:,Element2));
        ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
        ContinuousMatrix(Element2,:) = 0;
        ForceRelationMatrix(RowElement1,:) = ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:) * SlidingLimit(i);
        ForceRelationMatrix(Element2,:) = 0;
    end
end
end

function OutputResult(ElementTable,NodeTable)
ElementNumber = size(ElementTable,1);
figure(1); clf;
for i = 1 :1 : ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    x = [NodeTable(n,2),NodeTable(m,2)];
    y = [NodeTable(n,3),NodeTable(m,3)];
    z = [NodeTable(n,4),NodeTable(m,4)];
    plot3(x,y,z,'-black','LineWidth',4);
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);
axis on
end

function result = LEN(a,b,c,d,e,f)
result = ((a-d)^2 + (b-e)^2+(c-f)^2)^0.5;
end

function result = SolveThreeNodeAngle(NodeTable,Node1,Node2,Node3)
    a = NodeTable(Node2,2:4) - NodeTable(Node1,2:4);
    b = NodeTable(Node2,2:4) - NodeTable(Node3,2:4);
    result = acos(-a*b'/(norm(a)*norm(b)));
end

function result = StiffnessOpterator(A,B)
    [mA,nA] = size(A);
    result = zeros(mA,1);
    for i = 1 : 1 : mA
        for j = 1 : 1 : nA
            if (result(i) * A(i,j) ~=0)
                result(i) =  (1/result(i)+A(i,j)/B(j))^(-1);
            end
            if (result(i) == 0 && A(i,j)~=0)
                result(i) = B(j)/A(i,j);
            end
            if (result(i) ~= 0 && A(i,j)==0)
                result(i) = result(i);
            end
            if (result(i) == 0 && A(i,j)==0)
                result(i) = 0;
            end
        end
    end
end

function PlotResponse(Displacement,FinalForce,ElementNumber,UnfixedNodeNumber)
figure;
set(gcf,'position',[597,500,800,300])
plot(1:UnfixedNodeNumber,Displacement(2*UnfixedNodeNumber+1:3*UnfixedNodeNumber),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5);
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Node number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Displacment(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([1,18]);

figure;
set(gcf,'position',[597,500,800,300])
plot(1:ElementNumber,FinalForce,'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5)
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Force(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,65]);
end

function PlotResponseRange(MaxmimalForce,MinimalForce,MaximalDeformation,MinimalDeformation,ElementNumber,UnfixedNodeNumber)
figure;
set(gcf,'position',[597,500,800,300])
fill([1:ElementNumber,ElementNumber:-1:1],[MaxmimalForce;rot90(MinimalForce,2)]','r','FaceColor','[1 0.8 0.8]','EdgeColor','[1 0.8 0.8]')
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Force(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,65]);

figure;
set(gcf,'position',[597,500,800,300])
fill([1:UnfixedNodeNumber,UnfixedNodeNumber:-1:1],[MaximalDeformation(2*UnfixedNodeNumber+1:3*UnfixedNodeNumber);rot90(MinimalDeformation(2*UnfixedNodeNumber+1:3*UnfixedNodeNumber),2)]','r','FaceColor','[1 0.8 0.8]','EdgeColor','[1 0.8 0.8]','LineWidth',1)
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Node number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Displacment(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([1,18]);
end