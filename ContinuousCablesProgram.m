%% Mechanical analysis of frictional continuous cable system considering the infuluence of load path
function ContinuousCablesProgram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic information:
    % Function: Compute the load responses of frictional continuous cable systems
    % Author: Xue Yu
    % Institution: Space Structure Research Centre, Zhejiang University
    % Email: 11812048@zju.edu.cn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input parameters
FileName = 'Levy1';

AnalysisType = 1;
% 1:definite load path
% 0:indefinite load path

% FrictionCoefficient = 0;

% for mm = 1 : 1 : 4
%     AA =  [0,0.01,0.1,100];
FrictionCoefficient = 0.1;
NodeTable = xlsread(FileName,'Node');
ElementTable = xlsread(FileName,'Element');
LoadTable = xlsread(FileName,'Load');
ContinuousCableSequence = xlsread(FileName,'Continuous');
OutputResult(ElementTable,NodeTable)
% Input structural information via excel file
% two example(a Levy dome and a Geiger dome) are provided
% Example : Levy dome with continuous loop cables
% Sheet "Node":
%   column 1 :      node numbering
%   columns 2-4:    nodal coordinates in x-, y-, and z- directions(Unit: m)
%   column 5:       "1" means unfixed node, and "0" means fixed node
% Sheet "Element":
%   column 1:       segment numnering
%   column 2：      segment type；"0" means cable and "1" means strut
%   column 3：      segment axial stiffness (Unit: kN)
%   columns 4-5：   connecting nodes
%   column 6:       segment initial forces(Unit:kN)
%   column 7:       segment lengths;
% Sheet "Load":
%   column i:       external load on corresponding DOF in the ith load stage(Unit:kN)
% Sheet "Continuous":
%   column 1-2:     continuous segments

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
L0 = zeros(ElementNumber,ElementNumber);
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
    L0(i,i) = L(i,i)/(1+ElementTable(i,6)/ElementStiffness(i));
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
K_e = diag(diag(ElementStiffness)/L0);
% Segmenet stiffness vector
SegmentPrestress = ElementTable(:,6);
InitialElongations = L0 * diag(ElementStiffness)^(-1) * SegmentPrestress;
InitialElongations1 = InitialElongations;

SlidingLimit = GetSlidingLimit(NodeTable,ElementTable,ContinuousCableSequence,FrictionCoefficient);
% Compute sliding limits of the continuous segment force ratios
NodeTable1 = NodeTable;

%% compute structural response ranges when the load path is definite

if (AnalysisType == 1)
    LoadStageNumber = size(LoadTable,2);
    SegmentForce(:,1) = SegmentPrestress;
    DisplacementInCurrentStage = zeros(3*UnfixedNodeNumber,1);
    for c = 1 : 1 : LoadStageNumber
        NodeTable(1:UnfixedNodeNumber,2) = NodeTable1(1:UnfixedNodeNumber,2) + DisplacementInCurrentStage(1:UnfixedNodeNumber);
        NodeTable(1:UnfixedNodeNumber,3) = NodeTable1(1:UnfixedNodeNumber,3) + DisplacementInCurrentStage(UnfixedNodeNumber + 1 : 2 * UnfixedNodeNumber);
        NodeTable(1:UnfixedNodeNumber,4) = NodeTable1(1:UnfixedNodeNumber,4) + DisplacementInCurrentStage(2 * UnfixedNodeNumber + 1 : 3 * UnfixedNodeNumber);
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
            L0(i,i) = L(i,i)/(1+ElementTable(i,6)/ElementStiffness(i));
        end
        D0 = [C'* dX/L;
            C'* dY /L;
            C'* dZ /L];
        ZERO = diag([NodeTable(:,5);
            NodeTable(:,5);
            NodeTable(:,5)]);
        ZERO (all(ZERO == 0, 2),:) = [];
        D = ZERO * D0;
        if (c == 1)
            UnbalancedForce(:,c) = D * ElementTable(:,6);
        end
        if (c > 1)
            UnbalancedForce(:,c) = D * ElementTable(:,6) - sum(LoadTable(:,1:c-1),2);
        end
        K_e = diag(diag(ElementStiffness)/L0);
        InitialElongations = L0 * diag(ElementStiffness)^(-1) * SegmentPrestress;
        ResidualError(c) = norm(UnbalancedForce(:,c))/(3*UnfixedNodeNumber);
        LoadinCurrentStage = - UnbalancedForce(:,c) + LoadTable(:,c);
        ContinuousMatrix = eye(size(K_e,1));
        ForceRelationMatrix = eye(size(K_e,1));
        SegmentForce = zeros(ElementNumber,size(ContinuousCableSequence,1) + 2);
        for k = 1 : 1 : size(ContinuousCableSequence,1) + 2
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
            NodalDisplacment(:,k) = (K_E+K_G)^(-1) * (LoadinCurrentStage - K_F);
            SegmentForce(:,k+1) = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * (D' * NodalDisplacment(:,k) + InitialElongations);
            IsSlide(SegmentForce(:,k+1),ContinuousCableSequence,SlidingLimit)
            if (IsSlide(SegmentForce(:,k+1),ContinuousCableSequence,SlidingLimit))
                [ContinuousMatrix,ForceRelationMatrix] = GetMaxtrix(SegmentPrestress,SegmentForce(:,k+1),ContinuousCableSequence,SlidingLimit,ContinuousMatrix,ForceRelationMatrix);
            end
        end
        DisplacementInCurrentStage = DisplacementInCurrentStage + NodalDisplacment(:,k);
        ElementTable(:,6) = SegmentForce(:,k+1);
        force(:,c) = SegmentForce(:,k+1);
        dis(:,c) = DisplacementInCurrentStage;
        disn(c) = norm(NodalDisplacment(:,k));
        SegmentPrestress = ElementTable(:,6);
        SlidingLength(:,c) = D'* NodalDisplacment(:,k) - diag(K_e)^(-1) * ElementTable(:,6) + InitialElongations;
        %     InitialElongations = diag(K_e)^(-1) * SegmentPrestress;
        TotalSlidingLength(:,c) = sum(SlidingLength,2);
        AllSlidingLimit(:,c) = SlidingLimit;
    end
    % figure;plot(ResidualError);
    % residual errors
    PlotResponse(DisplacementInCurrentStage,force(:,c),ElementNumber,UnfixedNodeNumber);
%     OutputResult(ElementTable,NodeTable);
%     PlotSlidingStage(TotalSlidingLength,force,ElementTable,ContinuousCableSequence,AllSlidingLimit);
%     plot sliding-induced length of continuous segments
end

%% compute structural response ranges when the load path is indefinite
if (AnalysisType == 0)
    FinalLoad = LoadTable(:,LoadStageNumber);
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
        % Redistributed force density matrix
        G = C'* Q * C;
        K_G0 = zeros(3 * NodeNumber,3 * NodeNumber);
        K_G0(1:NodeNumber,1:NodeNumber) = G;
        K_G0(NodeNumber + 1:2 * NodeNumber,NodeNumber + 1:2 * NodeNumber) = G;
        K_G0(2 * NodeNumber + 1:3 * NodeNumber,2 * NodeNumber + 1:3 * NodeNumber) = G;
        K_G = ZERO * K_G0* ZERO';
        K_F = D * (RedistributionForce - SegmentPrestress);
        NodalDisplacment = (K_E+K_G)^(-1) * (FinalLoad - K_F);
        FinalForce = ForceRelationMatrix' * diag(StiffnessOpterator(ForceRelationMatrix,K_e)) * ContinuousMatrix * (D' * NodalDisplacment + InitialElongations);
        MaxmimalForce = max(MaxmimalForce,FinalForce);
        MinimalForce = min(MinimalForce,FinalForce);
        MaximalDeformation = max(MaximalDeformation,NodalDisplacment);
        MinimalDeformation = min(MinimalDeformation,NodalDisplacment);
    end
    PlotResponseRange(MaxmimalForce,MinimalForce,MaximalDeformation,MinimalDeformation,ElementNumber,UnfixedNodeNumber);
end
end

function PlotSlidingStage(TotalSlidingLength,force,ElementTable,ContinuousCableSequence,AllSlidingLimit)

for i = 1 : size(ContinuousCableSequence,1)
    Element1 = ContinuousCableSequence(i,1);
    Element2 = ContinuousCableSequence(i,2);
    ForceElement1 = force(Element1,:);
    ForceElement2 = force(Element2,:);
    ForceRatio = diag(ForceElement2)^(-1) * ForceElement1';
    figure;
    yyaxis left;
    plot(TotalSlidingLength(Element1,1:10:size(force,2)),'--^','MarkerSize',8,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',1),hold on
    ylabel('Sliding length(cm)', 'Color','black','fontsize',18,'FontName','Times New Roman');
    set(gca,'ycolor','black','FontSize',16,'FontName','Times New Roman');
    yyaxis right;
    plot(AllSlidingLimit(i,:),'LineWidth',1.5),hold on;
    plot((AllSlidingLimit(i,:)).^(-1),'LineWidth',1.5),hold on
    plot(ForceRatio(1:size(force,2)),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1);hold on;
    ylabel('Force ratio', 'Color','black','fontsize',18,'FontName','Times New Roman');
    set(gca,'ycolor','black','FontSize',16,'FontName','Times New Roman');
    xlabel('Load(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
    legend('Sliding length','Sliding limit range','Force ratio');
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
    if(ForceRatio(i) > 1.000001 * SlidingLimit(i))
        SlideTime(i) = (InitialSegmentForce(Element1) - SlidingLimit(i) * InitialSegmentForce(Element2))/(SlidingLimit(i) * DeltaSegmentForce(Element2) - DeltaSegmentForce(Element1));
    end
    if(ForceRatio(i) < 0.999999/SlidingLimit(i))
        SlideTime(i) = (InitialSegmentForce(Element1) - 1/SlidingLimit(i) * InitialSegmentForce(Element2))/(1/SlidingLimit(i) * DeltaSegmentForce(Element2) - DeltaSegmentForce(Element1));
    end
end
    [~,SlidingElement] = min(SlideTime);
    disp(min(SlideTime))
    if (min(SlideTime)<1)
        Element1 = ContinuousCableSequence(SlidingElement,1);
        Element2 = ContinuousCableSequence(SlidingElement,2);
        if(ForceRatio(SlidingElement) > 1)
            RowElement1 = find(ContinuousMatrix(:,Element1));
            RowElement2 = find(ContinuousMatrix(:,Element2));
            ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
            ContinuousMatrix(RowElement2,:) = 0;
            ForceRelationMatrix(RowElement1,:) =  ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:)/SlidingLimit(SlidingElement);
            ForceRelationMatrix(RowElement2,:) = 0;
        end
        if(ForceRatio(SlidingElement) < 1)
            RowElement1 = find(ContinuousMatrix(:,Element1));
            RowElement2 = find(ContinuousMatrix(:,Element2));
            ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
            ContinuousMatrix(RowElement2,:) = 0;
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
    if(ForceRatio(i) > 1.0001 * SlidingLimit(i) || ForceRatio(i) < 0.9999/SlidingLimit(i))
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
        ContinuousMatrix(RowElement2,:) = 0;
        ForceRelationMatrix(RowElement1,:) =  ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:)/SlidingLimit(i);
        ForceRelationMatrix(RowElement2,:) = 0;
    end
    if(Index(i) == 0)
        RowElement1 = find(ContinuousMatrix(:,Element1));
        RowElement2 = find(ContinuousMatrix(:,Element2));
        ContinuousMatrix(RowElement1,:) = ContinuousMatrix(RowElement1,:) + ContinuousMatrix(RowElement2,:);
        ContinuousMatrix(RowElement2,:) = 0;
        ForceRelationMatrix(RowElement1,:) = ForceRelationMatrix(RowElement1,:) + ForceRelationMatrix(RowElement1,Element1)/ForceRelationMatrix(RowElement2,Element2) * ForceRelationMatrix(Element2,:) * SlidingLimit(i);
        ForceRelationMatrix(RowElement2,:) = 0;
    end
end
end

function OutputResult(ElementTable,NodeTable)
ElementNumber = size(ElementTable,1);
NodeNumber = size(NodeTable,1);
figure(1); hold on;
for i = 1 :1 : ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    x = [NodeTable(n,2),NodeTable(m,2)];
    y = [NodeTable(n,3),NodeTable(m,3)];
    z = [NodeTable(n,4),NodeTable(m,4)];
    if(ElementTable(i,2) == 1)
        plot3(x,y,z,'-r','LineWidth',3);
    end
    if(ElementTable(i,2) == 0)
        plot3(x,y,z,'-b','LineWidth',1);
    end
    hold on;
end
set(gca,'DataAspectRatio',[1 1 1]);
axis on
for i = 1 : 1 : NodeNumber
    text(NodeTable(i,2),NodeTable(i,3),NodeTable(i,4),num2str(i));
end
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
plot(1:UnfixedNodeNumber,Displacement(2 * UnfixedNodeNumber+1:3 * UnfixedNodeNumber),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5);
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Node number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Displacment(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([1,UnfixedNodeNumber]);

figure;
set(gcf,'position',[597,500,800,300])
plot(1:ElementNumber,FinalForce,'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5)
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Force(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,ElementNumber]);
end

function PlotResponseRange(MaxmimalForce,MinimalForce,MaximalDeformation,MinimalDeformation,ElementNumber,UnfixedNodeNumber)
figure;
set(gcf,'position',[597,500,800,300])
fill([1:ElementNumber,ElementNumber:-1:1],[MaxmimalForce;rot90(MinimalForce,2)]','r','FaceColor','[1 0.8 0.8]','EdgeColor','[1 0.8 0.8]')
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Force(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,ElementNumber]);

figure;
set(gcf,'position',[597,500,800,300])
fill([1:UnfixedNodeNumber,UnfixedNodeNumber:-1:1],[MaximalDeformation(2*UnfixedNodeNumber+1:3*UnfixedNodeNumber);rot90(MinimalDeformation(2*UnfixedNodeNumber+1:3*UnfixedNodeNumber),2)]','r','FaceColor','[1 0.8 0.8]','EdgeColor','[1 0.8 0.8]','LineWidth',1)
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Node number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Displacment(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([1,UnfixedNodeNumber]);
end
