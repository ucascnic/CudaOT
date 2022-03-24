clear
flies = { 'automn1.jpg' , 'automn2.jpg' , 'automn3.jpg' , 'automn4.jpg' , 'automn5.jpg' , 'automn6.jpg',...
    'automn7.jpg' , 'automn8.jpg' , 'automn9.jpg' , 'automn10.jpg' , 'automn11.jpg' , 'automn12.jpg',...
    'automn13.jpg' , 'automn14.jpg' , 'automn15.jpg' , 'automn16.jpg' , 'automn17.jpg' , 'automn18.jpg',...
    'picture1.jpg' , 'picture2.jpg','picture3.jpg' , 'picture4.jpg' , 'picture5.jpg' , 'picture6.jpg' ,'winter1.jpg' , 'winter2.jpg' , ...
    'winter3.jpg' , 'winter4.jpg' , 'winter5.jpg' , 'winter6.jpg'};


%%
load('output_distmatrix.mat','Wass_xy')

Wass_xy = (Wass_xy + Wass_xy')/2;

% Wass_xy = Wass_xy([1:6,25:30],[1:6,25:30]);
[wcoeff,score,latent,~,explained] = pca(Wass_xy,'VariableWeights','variance');

% [eigs,s] = eig(xx);
data_PCA = mdscale(Wass_xy,3);
x = data_PCA(:,1);
y = data_PCA(:,2);

% scatter3(x,y,z)
%%

%%
% scatter(x,y,'fill')

width = 1920;
hight = 1200;

scalex = 4e4;
scaley = scalex*2;
h = cell(length(flies),1);
lim = [-scalex,scalex]/10;



n = length(x);
rng(10)

%%

 

x0 = 1920:-1:1;
y0 = 1200:-1:1;

hold on
for i =  n:-1:1
    f = flies{i};
    f = imread(f);
    
    h{i} = imshow(f); 
    hold on
    
    axis on
    
    h{i}.XData = [x(i)*scalex,x(i)*scalex+width];
    h{i}.YData = [y(i)*scaley,y(i)*scaley+hight];
 
end

%% plot the distance between the min pictures
Wass_xy = Wass_xy + eye(n,n)*1e10;
s = cell(n,1);
for i = 1:n
    [~,inds] = sort(Wass_xy(i,:));
    cnt = 0;
    for j = 1:1:2
        ind = inds(j);
        minv = min(i,ind);
        maxv = i+ind-minv;
        if  isempty(find(s{minv}==maxv, 1))
            s{minv} = [s{minv} maxv];
            cnt = cnt + 1;
            xx = [x(i)*scalex+width/2,x(ind)*scalex+width/2];
            yy = [y(i)*scaley+hight/2,y(ind)*scaley+hight/2];
            line(xx,yy,'linewidth',0.5)
            hold on
            b = sprintf('%.2f',Wass_xy(i,ind));
            text(mean(xx),mean(yy)-100,b,'fontname','New Times','fontsize',8);
            
        end
        if cnt == 2
            break
        end
    end
    
    
end
hold on


% 
%%
axis on
% set(gcf,'ToolBar','none','ReSize','off');   % remove the toolbar
set(gcf,'color','w');
xlim([-10000,12000])
ylim([-10000,12000])
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
box on
% print(gcf, 'test2d.png','-r600','-dpng');