function [h,ax,BigAx,patches,pax,cc] = plotmatrix_lower_david(varargin)
%PLOTMATRIX_LOWER Scatter plot matrix.
%   PLOTMATRIX(X,Y) scatter plots the columns of X against the columns
%   of Y.  If X is P-by-M and Y is P-by-N, PLOTMATRIX will produce a
%   N-by-M matrix of axes. PLOTMATRIX(Y) is the same as PLOTMATRIX(Y,Y)
%   except that the diagonal will be replaced by HIST(Y(:,i)). 
%
%   PLOTMATRIX(...,'LineSpec') uses the given line specification in the
%   string 'LineSpec'; '.' is the default (see PLOT for possibilities). 
%
%   PLOTMATRIX(...,Z) plots the optimum in each of the scatter plots.
%   For the histogram it is represented by a line. The number of columns X,Y of
%   the data needs to be equal to length of the vector inputed as Z
%
%   PLOTMATRIX(AX,...) uses AX as the BigAx instead of GCA.
%
%   [H,AX,BigAx,P,PAx,cc] = PLOTMATRIX(...) returns a matrix of handles
%   to the objects created in H, a matrix of handles to the individual
%   subaxes in AX, a handle to big (invisible) axes that frame the
%   subaxes in BigAx, a matrix of handles for the histogram plots in
%   P, and a matrix of handles for invisible axes that control the
%   histogram axes scales in PAx.  BigAx is left as the CurrentAxes so
%   that a subsequent TITLE, XLABEL, or YLABEL will be centered with
%   respect to the matrix of axes. CC is colorbar handle
%
%   Example:
%       x = randn(50,3); y = x*[-1 2 1;2 0 1;1 -2 3;]';
%       plotmatrix(y)

%   Clay M. Thompson 10-3-94
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.19.4.6 $  $Date: 2005/06/21 19:37:51 $
%
% modifications
% Bekaert David 09/2012     set the scatter plots in color
% Bekaert David 09/2012     Fix issue with the axis limits of the
%                           histograms.
% Bekaert David 09/2012     Allow the histograms to have identical y-axis.
%                           To modify this change the variable fix_ylim_hist.
% Bekaert David 09/2012     Allow extra input argument to plot the optimum
%                           of each dataset
% Bekaert David 09/2013     Plot the scatter of components as a matrix 
% Bekaert David 02/2014     Include a matrix to be ploted
% Bekaert David 02/2014     Include freezeColors option to maintain
%                           different colormaps
% Bekaert David 04/2014     Fix the scater plot maximum colors and plto a
%                           colorbar and out put its handle.

% Parse possible Axes input
[cax,args,nargs] = axescheck(varargin{:});
%error(nargchk(1,6,nargs,'struct'));            % RUCH HACK - commented out
nin = nargs;
sym = 'sq'; % Default scatter plot symbol.
dohist = 0;
fix_ylim_hist=1;
ncolors=300;
matrix_scatter_plot = 0;
markersize = 2;
matrix_plot = 1;
%               flipud(hot)

n_colors = 20;
colormaps_used_str = 'jet';
              
optimum_flag =0;
roughness_flag =0;
max_smoothness_flag =0;
counter = 1;
if nin>1
   for k=2:nin
      if isnumeric(args{k})
          if counter ==1
                optimum_data = args{k};
                optimum_flag =1;
          elseif counter ==2
                max_smoothness_data = args{k};
                max_smoothness_flag =1;
          elseif counter==3
                roughest_data = args{k};
                roughness_flag =1;
          end
          counter = counter +1;

          if k==nin
              nin = nin-counter+1;
          end
      end
   end
end


if ischar(args{nin}),
  sym = args{nin};
  if ~strcmpi(sym,'contour') && ~strcmpi(sym,'plot_color')
      [l,c,m,msg] = colstyle(sym); %#ok
      if ~isempty(msg), error(msg); end %#ok
  end
  nin = nin - 1;
end

if nin==1, % plotmatrix(y)
  rows = size(args{1},2); cols = rows;
  x = args{1}; y = args{1};
  dohist = 1;
elseif nin==2, % plotmatrix(x,y)
  rows = size(args{2},2); cols = size(args{1},2);
  x = args{1}; y = args{2};
else
  error('MATLAB:plotmatrix:InvalidLineSpec',...
        'Invalid marker specification. Type ''help plot''.');
end
x = double(x);
y = double(y);

% the maximum number of data points which ever could get ina  single
% histogram
max_data = max([size(x,1) size(y,1)]);

% Don't plot anything if either x or y is empty
patches = [];
pax = [];
if isempty(rows) || isempty(cols),
   if nargout>0, h = []; ax = []; BigAx = []; end
   return
end

if ndims(x)>2 || ndims(y)>2,
  error(id('InvalidXYMatrices'),'X and Y must be 2-D.')
end
if size(x,1)~=size(y,1) || size(x,3)~=size(y,3),
  error(id('XYSizeMismatch'),'X and Y must have the same number of rows and pages.');
end


% Create/find BigAx and make it invisible
BigAx = newplot(cax);
fig = ancestor(BigAx,'figure');
hold_state = ishold(BigAx);
set(BigAx,'Visible','off','color','none')

if any(sym=='.'),
  units = get(BigAx,'units');
  set(BigAx,'units','pixels');
  pos = get(BigAx,'Position');
  set(BigAx,'units',units);
  markersize = max(1,min(15,round(15*min(pos(3:4))/max(1,size(x,1))/max(rows,cols))));
else
  markersize = get(0,'defaultlinemarkersize');
end

% Create and plot into axes
ax = zeros(rows,cols);
pos = get(BigAx,'Position');
width = pos(3)/cols;
height = pos(4)/rows;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
m = size(y,1);
k = size(y,3);
xlim = nan([rows cols 2]);
ylim = nan([rows cols 2]);
BigAxHV = get(BigAx,'HandleVisibility');
hh=[];




% putting the same colorbar to all the scatter plots requires to now what
% the maximum is prior ot plotting:
max_value = [];
for i=1:1:rows-1,
    for j=1:1:i,

%         n_bins=min([round(m/40),100]);
%         n_bins2=min([round(sqrt(m/4)),50]);
        n_bins=32;
        n_bins2 = n_bins;

        % dummy run on the histogram (2D)
        [h_temp,bins_temp] = hist3([reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k])],[n_bins2,n_bins2]);

        % matrix or scatterplot
        data_matrix = h_temp;
        temp = sort(reshape(data_matrix,[],1));
        temp(isnan(temp)==1)=[];
        ix = ceil(0.97*length(temp));
        max_value = max([max_value temp(ix)]);
        
    end
end



for i=1:1:rows,
  for j=1:1:i,
    axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height width*(1-space) height*(1-space)];
    findax = findobj(fig,'Type','axes','Position',axPos);
    if isempty(findax),
      ax(i,j) = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',fig);
      set(ax(i,j),'visible','on');
    else
      ax(i,j) = findax(1);
      
    end
%     n_bins=min([round(m/40),100]);
%     %n_bins=15
%     n_bins2=min([round(sqrt(m/4)),50]);
%     %n_bins2=15
%     
%     
%     n_bins=32;
%     n_bins2 = n_bins;
    
    
    if i==rows & dohist==1;
      colormap(jet);
      [nn,xx] = hist(reshape(y(:,j,:),[m k]),n_bins);
      patches(i,:) = bar(ax(i,j),xx,nn,1,'b');
      if matrix_scatter_plot~=1;
        shading flat;
      end
      %set(ax(i,j),'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
      %set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
      pax(j) = ax(i,j);  % ax handles for histograms
      axis(ax(i,j),'tight');
      box on;
      ylim(i,j,:) = get(ax(i,j),'ylim');
      
      
      if optimum_flag==1;
          % plotting the optimum values incase they are specified
          hold on;
          plot([optimum_data(j) optimum_data(j)],[0 max_data],'r-','linewidth',2);
      end
      if max_smoothness_flag ==1;
         % plotting the maximum smoothness probability solution if its specified
         hold on;
         plot([max_smoothness_data(j) max_smoothness_data(j)],[0 max_data],'m-','linewidth',2);
      end
      if roughness_flag ==1;
         % plotting the roughest solution if its specified
         hold on;
         plot([roughest_data(j) roughest_data(j)],[0 max_data],'g-','linewidth',2);
      end
 
      
      % fixing the colorbar and colormap
      freezeColors;
      
    elseif strcmpi (sym,'contour');
      [h,bins] = hist3([reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k])],[n_bins2,n_bins2]);
      contour(bins{1},bins{2},h','parent',ax(i,j))';
      axis tight;
      ylim(i,j,:) = get(ax(i,j),'ylim');
      box on;
    elseif strcmpi (sym,'plot_color')
      
      
      [h_temp,bins_temp] = hist3([reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k])],[n_bins2,n_bins2]);
      
        
      if matrix_scatter_plot ==1;
    %       imagesc([min(bins_temp{1}) max(bins_temp{1})],[min(bins_temp{2}) max(bins_temp{2})],h_temp)
%           Jet_colors = jet(max(max(h_temp)));
%           Jet_colors = [[1 1 1]; Jet_colors(70)];
    %       ix_colors = round(h_temp./max(max(h_temp))*ncolors);
    
    
    
            h_temp = h_temp';
          h_sorted = sort(h_temp);
          h_normalise = h_sorted(end-floor(length(h_sorted)*0.025));
          imagesc([min(bins_temp{1}) max(bins_temp{1})],[min(bins_temp{2}) max(bins_temp{2})],(h_temp./h_normalise).*ncolors)
          Jet_colors = [[1 1 1]; jet(ncolors)];
          colormap(Jet_colors)
     
          axis xy
          
          fprintf('matrix scatter plot \n')
      
      else

          
          % matrix or scatterplot
          if matrix_plot ==1;
              n_bins2_temp = n_bins2;
              data_matrix = h_temp;
             % imagesc(bins_temp{1},bins_temp{2},data_matrix)
             % axis xy

              % saturating the scatter plot
              data_matrix(data_matrix>=max_value)=max_value;
             
              % reshifting the zeros such they can get white asociated with
              % them.
              max_new_zero = max_value*(n_colors+1)/(n_colors);
              data_matrix(data_matrix==0)=max_new_zero;

              % plotting the result
              imagesc((bins_temp{1}),(bins_temp{2}),data_matrix')
              axis xy
              caxis([0 max_new_zero])
              colormap([eval([colormaps_used_str,'(',num2str(n_colors),')']) ; 1 1 1])
              freezeColors

              %fprintf('matrix plot showing all different from 0\n');
                         
             
%              
%              
%              figure
%               data_matrix = h_temp;
% 
%               max_new_zero = ceil(max(max(data_matrix))*n_colors/(n_colors-1));
%               data_matrix(data_matrix==0)=max_new_zero;
% 
%               imagesc((bins_temp{1}),(bins_temp{2}),data_matrix')
%               axis xy
%               colormap([eval([colormaps_used_str,'(',num2str(n_colors),')']) ; 1 1 1])
%               freezeColors
% 
%               fprintf('matrix plot showing all different from 0\n')

              
          else
              
              % copying the outer bins to extend the convex hull such no points
              % have NaN values.
              h_temp_new = [h_temp(:,1) h_temp h_temp(:,end)];
              h_temp_new_new = [h_temp_new(1,:); h_temp_new; h_temp_new(end,:)];
              bins_temp_new{1} = [bins_temp{1}(1)-(bins_temp{1}(2)-bins_temp{1}(1))  bins_temp{1} bins_temp{1}(end)+(bins_temp{1}(2)-bins_temp{1}(1))];
              bins_temp_new{2} = [bins_temp{2}(1)-(bins_temp{2}(2)-bins_temp{2}(1))  bins_temp{2} bins_temp{2}(end)+(bins_temp{2}(2)-bins_temp{2}(1))];
              bins_temp = bins_temp_new;
              h_temp = h_temp_new_new;
              n_bins2_temp = n_bins2+2;

              C = NaN(size(reshape(x(:,j,:),[m k])));
              bins_matrix{1} = repmat(bins_temp{1}',1,n_bins2_temp);
              bins_matrix{2} = repmat(bins_temp{2},n_bins2_temp,1);
              bins_vector{1} = reshape(bins_matrix{1},[],1);
              bins_vector{2} = reshape(bins_matrix{2},[],1);
              h_temp = reshape(h_temp,[],1);
              ZI = griddata(bins_vector{1},bins_vector{2},h_temp,reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k]));

              %       hh(i,j,:) = scatter3(reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k]),ZI,markersize,ZI,'filled','parent',ax(i,j))';

              x_vector = reshape(x(:,j,:),[m k]);
              y_vector = reshape(y(:,i+1,:),[m k]);
              ix = find(isnan(ZI)==1);
              ZI(ix)=[];
              x_vector(ix)=[];
              y_vector(ix)=[];
              Jet_colors = jet(ncolors);

              ix_colors = round(ZI./max(ZI)*ncolors);
              ix_colors(ix_colors==0)=1;
              ZI_colors = Jet_colors(ix_colors,:);
              clear ix_colors Jet_colors

              for kk=1:length(x_vector)
                hh(i,j,:) = plot(x_vector(kk), y_vector(kk),'.','Color',ZI_colors(kk,:),'parent',ax(i,j))';
                hold on
              end
              hold off
              set(hh(i,j,:),'markersize',markersize);
          end
    
      end      
      
      
      if optimum_flag==1
          % plotting the optimum values incase they are specified
          hold on
%           plot(optimum_data(j), optimum_data(i+1),'k.','markersize',25)
          plot(optimum_data(j), optimum_data(i+1),'o','markersize',8,'MarkerFaceColor','r','markeredgecolor','k','linewidth',2)

      end
      if max_smoothness_flag ==1
         % plotting the maximum smoothness probability solution if its specified
         hold on
%          plot(max_smoothness_data(j), max_smoothness_data(i+1),'g.','markersize',25)
         plot(max_smoothness_data(j), max_smoothness_data(i+1),'o','markersize',8,'MarkerFaceColor','m','markeredgecolor','k','linewidth',2)

      end
      if roughness_flag ==1
         % plotting the roughest solution if its specified
         hold on
%          plot(roughest_data(j), roughest_data(i+1),'b.','markersize',25)
         plot(roughest_data(j), roughest_data(i+1),'o','markersize',8,'MarkerFaceColor','g','markeredgecolor','k','linewidth',2)

      end
      
      axis tight
      view(0,90)
      ylim(i,j,:) = get(ax(i,j),'ylim');
      box on
            
      % fixing the colorbar and colormap
      freezeColors
      
      
      
    else
      hh(i,j,:) = plot(reshape(x(:,j,:),[m k]),reshape(y(:,i+1,:),[m k]),sym,'parent',ax(i,j))';
      set(hh(i,j,:),'markersize',markersize);
      axis tight
      ylim(i,j,:) = get(ax(i,j),'ylim');
      box on
    end
    %set(ax(i,j),'xlimmode','auto','ylimmode','auto','xgrid','off','ygrid','off')
    set(ax(i,j),'xgrid','off','ygrid','off')
    set(gca,'fontweight','bold','fontsize',8,'YAxisLocation','left')
    xlim(i,j,:) = get(ax(i,j),'xlim');
    
    
    % add a colorbar to the plot
    if matrix_plot ==1 && strcmpi(sym,'plot_color') && i==rows && j==i
      
        i_temp = i-1;
        j_temp = i_temp+1;
        
        if i_temp~=1
            axPos = [pos(1)+(j_temp-1)*width pos(2)+(rows-i_temp)*height+height*0.1 width*(1-space)-width*0.1 height*(1-space)+height*0.5];
        else
            axPos = [pos(1)+0.1*width+(j_temp-1)*width pos(2)+(rows-i_temp)*height+height*0.1 width*(1-space)-width*0.1 height*(1-space)/2];
        end
        findax = findobj(fig,'Type','axes','Position',axPos);
        if isempty(findax),
          ax_colorbar = axes('Position',axPos,'HandleVisibility',BigAxHV,'parent',fig);
        else
          findax(1);
        end
        set(gca,'visible','off')

        caxis([0 max_new_zero])
        colormap([eval([colormaps_used_str,'(',num2str(n_colors),')']) ; 1 1 1])
        cc = colorbar('location','northoutside');
        title(cc,'Frequency','fontsize',12);
        xlabel_pos = floor(max_new_zero/(n_colors+1)*(n_colors)/100)*100;
        set(cc,'XTick',[0 xlabel_pos],'xticklabel',num2str([0 xlabel_pos]'))

        freezeColors        
    end

  end
end

xlimmin = nanmin(xlim(:,:,1),[],1); xlimmax = nanmax(xlim(:,:,2),[],1);
ylimmin = nanmin(ylim(:,:,1),[],2); ylimmax = nanmax(ylim(:,:,2),[],2);

% Try to be smart about axes limits and labels.  Set all the limits of a
% row or column to be the same and inset the tick marks by 10 percent.
inset = .15;

% fixing all the y-axes except hist
for i=1:rows-dohist,
  set(ax(i,i),'ylim',[ylimmin(i,1) ylimmax(i,1)])
  
  if optimum_flag==1
  end  
      
  
  dy = diff(get(ax(i,i),'ylim'))*inset;
  set(ax(i,ax(i,:)~=0),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
  %set(ax(i,ax(i+dohist,:)~=0),'ylim',[ylimmin(i,1) ylimmax(i,1)])
end

% fixing all the y-axes of the hist
if dohist==1 & fix_ylim_hist==1
  % fixing the y limits of the histogram to be equal for all of them
  set(ax(rows,rows),'ylim',[ylimmin(rows,1) ylimmax(rows,1)])  
  dy = diff(get(ax(rows,rows),'ylim'))*inset;
  set(ax(rows,ax(rows,:)~=0),'ylim',[ylimmin(rows,1) ylimmax(rows,1)+dy])
end

% fixing the x-axes of all figures
dx = zeros(1,cols);
for j=1:cols,
  set(ax(j,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
  if dohist==1
      dx(j) = diff(get(ax(j,j),'xlim'))*inset;
%     dx(j)=0;
  else
    dx(j) = diff(get(ax(j,j),'xlim'))*inset;
  end
  set(ax(ax(:,j)~=0,j),'xlim',[xlimmin(1,j)-dx(j) xlimmax(1,j)+dx(j)])
end


% 
% for i=1:rows-dohist
%     for j=1:cols,
% 
%         
%         keyboard
%           set(ax(i,j),'ylim',[ylimmin(i,1) ylimmax(i,1)])
%             color_limits = get(ax(i,j),'CLim');
%             n_colors = 59;
%             
%             
%           
%           colormap(jet(59))
%           
%           if optimum_flag==1
%           end  
% 
% 
%           dy = diff(get(ax(i,i),'ylim'))*inset;
%           set(ax(i,ax(i,:)~=0),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
%     end
% end




for i=1:rows-1
    set(ax(i,1:i),'xticklabel','')
end
if fix_ylim_hist==0
    set(ax(rows,1),'yticklabel','')

end
for j=2:cols
    set(ax(j:rows,j),'yticklabel','')
    %set(ax(1:j-1,j),'yticklabel','')
end


set(BigAx,'XTick',get(ax(1,1),'xtick'),'YTick',get(ax(1,1),'ytick'), ...
          'userdata',ax,'tag','PlotMatrixBigAx')

%if dohist, % Put a histogram on the diagonal for plotmatrix(y) case
%  for i=rows:-1:1,
%    histax = axes('Position',get(ax(i,i),'Position'),'HandleVisibility',BigAxHV,'parent',fig);
%    [nn,xx] = hist(reshape(y(:,i,:),[m k]),50);
%    patches(i,:) = bar(histax,xx,nn,'hist');
%    set(histax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
%    set(histax,'xlim',[xlimmin(1,i)-dx(i) xlimmax(1,i)+dx(i)])
%    pax(i) = histax;  % ax handles for histograms
%    %ax(i,i)=histax;
%  end
%  patches = patches';
%end

% Make BigAx the CurrentAxes
set(fig,'CurrentAx',BigAx)
if ~hold_state,
   set(fig,'NextPlot','replace')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
 'String','','Visible','on')

if nargout~=0,
  h = hh;
end
 
function str=id(str)
str = ['MATLAB:plotmatrix:' str];
