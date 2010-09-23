function [handles,order] = plotX(tr,varargin)
% set defaults
dispBranchLabels = false;
dispLeafLabels = false;
dispTerminalLabels = false;
renderType = 'phylogram';
orientation = 'left';
rotation = 0;

tr = struct(tr);
tr.numBranches = size(tr.tree,1);

if nargin>1 && islogical(varargin{1})
    activeBranches = varargin{1};  
    argStart = 2;
else     
    activeBranches = true(tr.numBranches,1);
    argStart = 1;
end

if nargin - argStart > 0
        if rem(nargin - argStart,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
              'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'type','orientation','rotation',...
              'branchlabels','leaflabels','terminallabels'};   
    for j = argStart:2:nargin-argStart          
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); %#ok
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                  'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                  'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1 % type
                    oktypes={'phylogram','cladogram','radial','dendrogram'};
                    l = strmatch(lower(pval),oktypes); %#ok
                    if isempty(l)
                         error('Bioinfo:UnknownTypeName',...
                               'Unknown option for %s.',upper(okargs{k}));
                    else
                        if l==4 
                            l=1; 
                        end
                        renderType = oktypes{l};
                    end 
                case 2 % orientation
                    oktypes={'left','right','top','bottom'};
                    l = strmatch(lower(pval),oktypes); %#ok
                    if isempty(l)
                         error('Bioinfo:UnknownOrientation',...
                               'Unknown option for %s.',upper(okargs{k}));
                    else
                        orientation = oktypes{l};
                    end 
                case 3 % rotation
                    if isreal(pval(1))
                        rotation = double(pval(1));
                    else
                        error('Bioinfo:NotValidType',...
                              'ROTATION must be numeric and real');
                    end
                case 4 % branch labels
                    dispBranchLabels = pval == true;
                case 5 % leaf labels
                    dispLeafLabels = pval == true;
                case 6 % terminal labels
                    dispTerminalLabels = pval == true;
            end
        end
    end
end    
 
tr = doBasicCalculations(tr,activeBranches,renderType);

nodeIndex   = 1:tr.numLabels;
leafIndex   = 1:tr.numLeaves;
branchIndex = tr.numLeaves+1:tr.numLabels;
    

% check empty names
for ind = nodeIndex
    if isempty(tr.names{ind}) 
        if ind > tr.numLeaves
            tr.names{ind} = ['Branch ' num2str(ind-tr.numLeaves)]; 
        else
            tr.names{ind} = ['Leaf ' num2str(ind)]; 
        end
    end
end

% rendering graphic objects
fig = gcf;%figure('Renderer','ZBuffer');
h.axes =gca;% axes; hold on;
hold on;
box off;
%axis off;

switch orientation
    case 'bottom'
set(gca,'XGrid','on');
set(gca,'YColor','w');
    case 'left'
set(gca,'YGrid','on');
set(gca,'XColor','w');
end

sepUnit = max(tr.x)*[-1/20 21/20];


% setting the axes
switch renderType
    case {'phylogram','cladogram'}
        switch orientation
            case 'left'
                set(h.axes,'YTick',1:numel(tr.terminalNodes),...
                           'YtickLabel','','YAxisLocation','Right')
                if dispTerminalLabels
                    %set(h.axes,'Position',[.05 .10 .7 .85])
                else
                    %set(h.axes,'Position',[.05 .10 .9 .85])
                end
                xlim(sepUnit); 
                ylim([0 numel(tr.terminalNodes)+1]);
            case 'right'
                set(h.axes,'YTick',1:numel(tr.terminalNodes),'Xdir','reverse',...
                           'YtickLabel','','YAxisLocation','Left')
                if dispTerminalLabels
                    %set(h.axes,'Position',[.25 .10 .7 .85])
                else
                    %set(h.axes,'Position',[.05 .10 .9 .85])
                end
                xlim(sepUnit); 
                ylim([0 numel(tr.terminalNodes)+1]);
            case 'top'
                set(h.axes,'XTick',1:numel(tr.terminalNodes),...
                           'XtickLabel','','XAxisLocation','Top')
                if dispTerminalLabels
                    %set(h.axes,'Position',[.10 .05 .85 .7])
                else
                    %set(h.axes,'Position',[.10 .05 .85 .9])
                end
                ylim(sepUnit); 
                xlim([0 numel(tr.terminalNodes)+1]);
            case 'bottom'
                set(gca,'XDir','reverse');
                set(h.axes,'XTick',1:numel(tr.terminalNodes),'Ydir','reverse',...
                           'XtickLabel','','XAxisLocation','Bottom')
                if dispTerminalLabels
                    %set(h.axes,'Position',[.10 .25 .85 .7])
                else
                    %set(h.axes,'Position',[.10 .05 .85 .9])
                end
                ylim(sepUnit); 
                xlim([0 numel(tr.terminalNodes)+1]);
        end
    case 'radial'
        set(h.axes,'XTick',[],'YTick',[])
        set(h.axes,'Position',[.05 .05 .9 .9])
        dispTerminalLabels = false;
        axis equal
end

% drawing lines
switch renderType
    case 'phylogram'
        X = tr.x([nodeIndex;repmat([tr.par(1:tr.numLabels-1) tr.numLabels],2,1)]);
        Y = tr.y([repmat(nodeIndex,2,1);[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        switch orientation
            case {'left','right'}
                h.BranchLines = plot(X,Y,'-k');
                delete(h.BranchLines(~tr.activeNodes))
                h.BranchLines = h.BranchLines(tr.activeNodes);
            case {'top','bottom'}
                h.BranchLines = plot(Y,X,'-k');
                delete(h.BranchLines(~tr.activeNodes))
                h.BranchLines = h.BranchLines(tr.activeNodes);                
        end
    case 'cladogram'
        X = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        Y = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
        switch orientation
            case {'left','right'}
                h.BranchLines = plot(X,Y,'-k');
                delete(h.BranchLines(~tr.activeNodes))
                h.BranchLines = h.BranchLines(tr.activeNodes);
            case {'top','bottom'}
                h.BranchLines = plot(Y,X,'-k');
                delete(h.BranchLines(~tr.activeNodes))
                h.BranchLines = h.BranchLines(tr.activeNodes);                
        end
    case 'radial'
       R = tr.x; 
       A = tr.y / numel(tr.terminalNodes)*2*pi+rotation*pi/180;
       tr.x = R .* sin(A);
       tr.y = R .* cos(A);
       X = tr.x([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]);
       Y = tr.y([nodeIndex;[tr.par(1:tr.numLabels-1) tr.numLabels]]); 
       h.BranchLines = plot(X,Y,'-k');
       delete(h.BranchLines(~tr.activeNodes))
       h.BranchLines = h.BranchLines(tr.activeNodes);
end
               
% drawing nodes
switch renderType
    case {'phylogram','cladogram','radial'}
        switch orientation
            case {'left','right'}
                %h.BranchDots = 
                for i=branchIndex(tr.activeNodes(branchIndex))
                    if ~strncmp(tr.names(i),'Branch',6)
                    res=sscanf(tr.names{i},'%d/%f');
                    %plot(tr.x(i),tr.y(i),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[1-res(2) 1-res(2) 1-res(2)]);
                    %text(tr.x(i)+0.01,tr.y(i),sprintf('%f',res(2)));
                    end
                end
                %plot(tr.x(branchIndex(tr.activeNodes(branchIndex))),...
                %    tr.y(branchIndex(tr.activeNodes(branchIndex))),'o',...
                %    'MarkerSize',5,'MarkerEdgeColor','k',...
                %    'MarkerFaceColor','b');
                h.LeafDots = plot(tr.x(leafIndex(tr.activeNodes(leafIndex))),...
                    tr.y(leafIndex(tr.activeNodes(leafIndex))),'square',...
                    'MarkerSize',4,'MarkerEdgeColor','k',...
                    'MarkerFaceColor','w');
            case {'top','bottom'}
                %h.BranchDots = plot(tr.y(branchIndex(tr.activeNodes(branchIndex))),...
                %    tr.x(branchIndex(tr.activeNodes(branchIndex))),'o',...
                %    'MarkerSize',5,'MarkerEdgeColor','k',...
                %    'MarkerFaceColor','b');
                h.LeafDots = plot(tr.y(leafIndex(tr.activeNodes(leafIndex))),...
                    tr.x(leafIndex(tr.activeNodes(leafIndex))),'square',...
                    'MarkerSize',4,'MarkerEdgeColor','k',...
                    'MarkerFaceColor','w');
        end
    case 'radial0'
        h.BranchDots = plot(tr.x(branchIndex(tr.activeNodes(branchIndex))),...
            tr.y(branchIndex(tr.activeNodes(branchIndex))),'o',...
            'MarkerSize',5,'MarkerEdgeColor','k',...
            'MarkerFaceColor','b');
        h.LeafDots = plot(tr.x(leafIndex(tr.activeNodes(leafIndex))),...
            tr.y(leafIndex(tr.activeNodes(leafIndex))),'square',...
            'MarkerSize',4,'MarkerEdgeColor','k',...
            'MarkerFaceColor','w');
end

% resize figure if needed
switch renderType
    case {'phylogram','cladogram'}
        switch orientation
            case {'left','right'}  
                %correctFigureSize(fig, 15 * numel(tr.terminalNodes),0);
                fontRatio = max(get(fig,'Position').*[0 0 0 1])/numel(tr.terminalNodes);
            case {'top','bottom'}
                %correctFigureSize(fig, 0, 15 * numel(tr.terminalNodes));
                fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
        end
    case 'radial'
        temp = 10/pi*numel(tr.terminalNodes);
        %correctFigureSize(fig,temp,temp);
        fontRatio = max(get(fig,'Position').*[0 0 1 0])/numel(tr.terminalNodes);
end

set(h.axes,'Fontsize',min(9,ceil(fontRatio/1.5)));  

% set branch node labels
X = tr.x(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));
Y = tr.y(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels)));
switch renderType
    case {'phylogram','cladogram'}
        switch orientation
            case {'left'}  
               h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
               set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
               set(h.branchNodeLabels,'vertical','bottom')
               set(h.branchNodeLabels,'horizontal','right')
               set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));  
            case {'right'} 
               h.branchNodeLabels = text(X+sepUnit(1)/2,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
               set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
               set(h.branchNodeLabels,'vertical','bottom')
               set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));  
            case {'top'} 
               h.branchNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
               set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
               set(h.branchNodeLabels,'vertical','bottom','Rotation',30)
               set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));  
            case {'bottom'} 
               h.branchNodeLabels = text(Y,X+sepUnit(1)/2,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
               set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
               set(h.branchNodeLabels,'vertical','bottom','Rotation',30)
               set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio/2)));  
        end
    case 'radial'
         h.branchNodeLabels = text(X,Y,tr.names(branchIndex(tr.activeNodes(tr.numLeaves+1:tr.numLabels))));
         set(h.branchNodeLabels,'color',[0 0 .8],'clipping','on')
         set(h.branchNodeLabels,'vertical','bottom')
         set(h.branchNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2)));  
         for ind = 1:numel(h.branchNodeLabels)
             if X(ind)<0
                set(h.branchNodeLabels(ind),'horizontal','right')
                set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')+[sepUnit(1)/2 0 0])
             else
                set(h.branchNodeLabels(ind),'horizontal','left')
                set(h.branchNodeLabels(ind),'Position',get(h.branchNodeLabels(ind),'Position')-[sepUnit(1)/2 0 0])
             end
         end
end

% set leaf nodes labels
X = tr.x(leafIndex(tr.activeNodes(1:tr.numLeaves)));
Y = tr.y(leafIndex(tr.activeNodes(1:tr.numLeaves)));
switch renderType
    case {'phylogram','cladogram'}
        switch orientation
            case {'left'}  
               h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
               order=tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves)));
               set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
               set(h.leafNodeLabels,'horizontal','left')
               set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2))); 
            case {'right'} 
               h.leafNodeLabels = text(X-sepUnit(1)/2,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
               set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
               set(h.leafNodeLabels,'horizontal','right')
               set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2))); 
            case {'top'} 
               h.leafNodeLabels = text(Y,X-sepUnit(1)/2,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
               set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
               set(h.leafNodeLabels,'horizontal','left','Rotation',60)
               set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2))); 
            case {'bottom'} 
               h.leafNodeLabels = text(Y,X-sepUnit(1),tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
               set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
               set(h.leafNodeLabels,'horizontal','right','Rotation',60)
               set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio/2))); 
        end
    case 'radial'
        h.leafNodeLabels = text(X,Y,tr.names(leafIndex(tr.activeNodes(1:tr.numLeaves))));
        set(h.leafNodeLabels,'color',[.5 .5 .5],'clipping','on')
        set(h.leafNodeLabels,'Fontsize',min(8,ceil(fontRatio*1.2))); 
       % textHeight = mean(cell2mat(get(h.leafNodeLabels,'Extent')))*[0 0 0 1]';
        for ind = 1:numel(h.leafNodeLabels)
            if X(ind)<0
                set(h.leafNodeLabels(ind),'horizontal','right')
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')+[sepUnit(1) 0 0])
             else
                set(h.leafNodeLabels(ind),'horizontal','left')
                set(h.leafNodeLabels(ind),'Position',get(h.leafNodeLabels(ind),'Position')-[sepUnit(1) 0 0])
            end
%             a=atan(Y(ind)/X(ind))*180/pi;
%             if a > 0  a = max(0,a-60)/2; else
%                       a = min(0,a+60)/2; end
%             set(h.leafNodeLabels(ind),'Rotation',a)
        end
            
end

% correct axis limits given the extent of labels
if dispBranchLabels
    E = cell2mat(get(h.branchNodeLabels,'Extent'));
    if strcmp(get(gca,'XDir'),'reverse')
        E(:,1) = E(:,1) - E(:,3);
    end
    if strcmp(get(gca,'YDir'),'reverse')
        E(:,2) = E(:,2) - E(:,4);
    end
    E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
    mins = min(E(:,[1 2]));
    maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
    axis([mins(1) maxs(1) mins(2) maxs(2)])
end

if dispLeafLabels
    E = cell2mat(get(h.leafNodeLabels,'Extent'));
    if strcmp(get(gca,'XDir'),'reverse')
        E(:,1) = E(:,1) - E(:,3);
    end
    if strcmp(get(gca,'YDir'),'reverse')
        E(:,2) = E(:,2) - E(:,4);
    end
    E=[E;[xlim*[1;0] ylim*[1;0] diff(xlim) diff(ylim)]];
    mins = min(E(:,[1 2]));
    maxs = max([sum(E(:,[1 3]),2) sum(E(:,[2 4]),2)]);
    axis([mins(1) maxs(1) mins(2) maxs(2)])
end

% set terminal nodes labels
switch renderType
    case {'phylogram','cladogram'}
        X = tr.x(tr.terminalNodes) * 0; 
        Y = tr.y(tr.terminalNodes);
        switch orientation
            case {'left'}  
               X = X + max(xlim) - sepUnit(1)/2;
               h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes)); 
            case {'right'}
               X = X + max(xlim) - sepUnit(1)/2;
               h.terminalNodeLabels = text(X,Y,tr.names(tr.terminalNodes));
               set(h.terminalNodeLabels,'Horizontal','right')
            case {'top'} 
               X = X + max(ylim) - sepUnit(1)/2;
               h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
               set(h.terminalNodeLabels,'Rotation',90)
            case {'bottom'} 
               X = X + max(ylim) - sepUnit(1)/2; 
               h.terminalNodeLabels = text(Y,X,tr.names(tr.terminalNodes));
               set(h.terminalNodeLabels,'Rotation',270)
         end
    case 'radial'
         h.terminalNodeLabels = text(0,0,' '); 
end

if dispTerminalLabels
    set(h.terminalNodeLabels,'Fontsize',min(9,ceil(fontRatio/1.5)));  
end

if ~dispBranchLabels 
    set(h.branchNodeLabels,'visible','off'); 
end
if ~dispLeafLabels 
    set(h.leafNodeLabels,'visible','off'); 
end
if ~dispTerminalLabels 
    set(h.terminalNodeLabels,'visible','off'); 
end

%box on

% store handles
set(fig,'UserData',h)   
if nargout
    handles = h;
end

%Revert the Xaxis
% m=max(tr.x(leafIndex(tr.activeNodes(leafIndex))));
% ticks=get(h.axes,'XTick');
%  tickslabel=str2num(get(h.axes,'XTickLabel'));
%  di=m-max(ticks);
%  nticks=ticks+di;
%  set(h.axes,'XTick',nticks);
%  tickslabel=tickslabel(end:-1:1);
%  set(h.axes,'XTickLabel',num2str(tickslabel));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = doBasicCalculations(tr,activeBranches,renderType)

% helper function to compute and find some features of the tree
tr.numLeaves = tr.numBranches + 1;
tr.numLabels = tr.numBranches + tr.numLeaves; 

% remove uderscores from names
for ind = 1:tr.numLabels
    tr.names{ind}(tr.names{ind}=='_')=' ';
end
    
% obtain parents for every node
tr.par(tr.tree(:)) = tr.numLeaves + [1:tr.numBranches 1:tr.numBranches];

% find active nodes
tr.activeNodes = true(tr.numLabels,1);
for ind =tr.numBranches:-1:1
    tr.activeNodes(tr.tree(ind,:)) = tr.activeNodes(tr.numLeaves+ind) & activeBranches(ind);
end

% propagate last leaf
tr.lastleaf = 1:tr.numLabels;
for ind = tr.numBranches:-1:1
    if ~tr.activeNodes(tr.tree(ind,1))
        tr.lastleaf(tr.tree(ind,:))=tr.lastleaf(ind+tr.numLeaves);
    end
end

tr.activeBranches = tr.activeNodes(tr.numLeaves+1:tr.numLabels)&activeBranches;
tr.activeLeaves = tr.activeNodes(1:tr.numLeaves);

% find x coordinates of branches
tr.x = tr.dist;
for i=1:length(tr.x)
    if tr.x(i)>50,tr.x(i)=0;end
end
for ind = tr.numBranches:-1:1
    tr.x(tr.tree(ind,:)) = tr.x(tr.tree(ind,:)) + tr.x(ind+tr.numLeaves);
end

% find y coordinates of branches 
tr.terminalNodes = tr.lastleaf([true,diff(tr.lastleaf(1:tr.numLeaves))~=0]);
tr.y=zeros(tr.numLabels,1);
tr.y(tr.terminalNodes)=1:length(tr.terminalNodes);
switch renderType
    case 'phylogram'
        for ind = 1:tr.numBranches
            if tr.activeBranches(ind)
                tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
            end
        end
    case {'cladogram','radial'}
        for ind = 1:tr.numBranches
            if tr.activeBranches(ind)
                if tr.x(tr.tree(ind,1))/tr.x(tr.tree(ind,2))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,1));
                elseif tr.x(tr.tree(ind,2))/tr.x(tr.tree(ind,1))>3
                    tr.y(ind+tr.numLeaves) = tr.y(tr.tree(ind,2));
                else
                    tr.y(ind+tr.numLeaves) = mean(tr.y(tr.tree(ind,:)));
                end
            end
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function correctFigureSize(fig,recommendedHeight,recommendedWidth)
% helper function to increase initial figure size depending on the screen &
% tree sizes
screenSize = diff(reshape(get(0,'ScreenSize'),2,2),[],2)-[0;100];
            % 100 gives extra space for the figure header and win toolbar
position = get(fig,'Position');
if recommendedHeight > position(4)
    if recommendedHeight < sum(position([2 4]))
        position(2) = sum(position([2 4])) - recommendedHeight;
        position(4) = recommendedHeight;
    elseif recommendedHeight < screenSize(2)
        position(2) = 30; 
        position(4) = recommendedHeight;
    else 
        position(2) = 30; 
        position(4) = screenSize(2);
    end
end
if recommendedWidth > position(3)
    if recommendedWidth < sum(position([1 3]))
        position(1) = sum(position([1 3])) - recommendedWidth;
        position(3) = recommendedWidth;
    elseif recommendedWidth < screenSize(1)
        position(1) = 0; 
        position(3) = recommendedHeight;
    else 
        position(1) = 0; 
        position(3) = screenSize(1);
    end
end    
set(fig,'Position',position)
