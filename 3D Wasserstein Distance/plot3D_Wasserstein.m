clear
flies = { 'automn1.jpg' , 'automn2.jpg' , 'automn3.jpg' , 'automn4.jpg' , 'automn5.jpg' , 'automn6.jpg',...
    'automn7.jpg' , 'automn8.jpg' , 'automn9.jpg' , 'automn10.jpg' , 'automn11.jpg' , 'automn12.jpg',...
    'automn13.jpg' , 'automn14.jpg' , 'automn15.jpg' , 'automn16.jpg' , 'automn17.jpg' , 'automn18.jpg',...
    'picture1.jpg' , 'picture2.jpg','picture3.jpg' , 'picture4.jpg' , 'picture5.jpg' , 'picture6.jpg' ,'winter1.jpg' , 'winter2.jpg' , ...
    'winter3.jpg' , 'winter4.jpg' , 'winter5.jpg' , 'winter6.jpg'};


%%
load('output_distmatrix.mat','Wass_xy')

Wass_xy = (Wass_xy + Wass_xy')/2;
[wcoeff,score,latent,~,explained] = pca(Wass_xy,'VariableWeights','variance');

% [eigs,s] = eig(xx);
data_PCA = Wass_xy*wcoeff;
x = data_PCA(:,1);
y = data_PCA(:,2);
z = data_PCA(:,3);
% scatter3(x,y,z)
%%

%%
% scatter(x,y,'fill')

width = 1920;
hight = 1200;

scalex = 3e5;
scaley = scalex;
h = cell(length(flies),1);
lim = [-scalex,scalex]/10;



n = length(x);
rng(10)

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
            zz = [z(i)*scalex,z(ind)*scalex];
            line(zz,xx,yy,'linewidth',0.5)
            hold on
            b = sprintf('%.2f',Wass_xy(i,ind));
            text((zz(1)+2*zz(2))/3,(xx(1)+2*xx(2))/3,mean(yy),b,'fontname','New Times','fontsize',7);
            
        end
        if cnt == 2
            break
        end
    end
    
    
end
hold on

%%

 

x0 = 1920:-1:1;
y0 = 1200:-1:1;
%%
hold on
for i = 1:30
    f = flies{i};
    f = imread(f);
    [xx,yy] = meshgrid(x0,y0);
    zz = ones(size(xx))*z(i)*scalex;
    
    xx = xx + x(i) * scalex;
    yy = yy + y(i) * scaley;
    mesh(zz,xx,yy,f);     
end
% 
%%
axis off
box off
% set(gcf,'ToolBar','none','ReSize','off');   % remove the toolbar
set(gcf,'color','w');



ylim([-9000,23000])
zlim([5000,30000])
xlim([-10000,9000])

[xx,yy]=meshgrid([-10000,9000],[-9000,23000]);
s = mesh(xx,yy,5000*ones(size(xx)));
s.FaceColor = [105,105,105]/255;
s.EdgeColor =  [105,105,105]/255;
s.FaceAlpha = 0.1;
hold on 
[yy,zz]=meshgrid([-9000,23000],[5000,30000]);
s = mesh(-10000*ones(size(yy)),yy,zz);
s.FaceColor = [105,105,105]/255;
s.EdgeColor =  [105,105,105]/255;
s.FaceAlpha = 0.1;

[xx,zz]=meshgrid([-10000,9000],[5000,30000]);
s = mesh(xx,23000*ones(size(yy)),zz);
s.FaceColor = [105,105,105]/255;
s.EdgeColor =  [105,105,105]/255;
s.FaceAlpha = 0.1;

 

view([12,-16,0])

set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
% print(gcf, 'test3.png','-r800','-dpng');