close all
clear all
run('data1.m') ;
run('comp2TEST.m') ;
run('pmc3TEST.m') ;

% data plot
figure(1)
i1 = plot([pP2C pP1C],[qP2C qP1C],'ko') ;
hold on
i2 = plot([pP2E pP1E],[-qP2E -qP1E],'ko') ;
hold on
i3 = plot([pP2o pP1o],[qP2o qP1o],'g*') ;
xlim([-30 400]);
ylim([-300 400]);
TFo = isempty(po) ;
TFe = isempty(pE) ;

if TFo==1
    legend([i1 i2],'Triaxial Compression data','Triaxial Extension data')
elseif TFe==1
    legend([i1 i3],'Triaxial Compression data','True Triaxial data')
elseif TFo==0 && TFe==0
    legend([i1 i2 i3],'Triaxial Compression data','Triaxial Extension data','True Triaxial data')
end
Xlab = xlabel('p [MPa]','Interpreter','latex','FontSize',12);
Ylab = ylabel('q [MPa]','Interpreter','latex','FontSize',12);
tit = title(name,'Interpreter','latex','FontSize',14);
set(tit,'FontSize',14);
set(Ylab,'FontSize',12);
set(Xlab,'FontSize',12);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% p-q plot
xfitP2 = linspace(-30,300) ;
xfitP1 = linspace(-400,400) ;

figure(2)

plot([pP2C pP1C],[qP2C qP1C],'ko') ;
hold on
plot([pP2E pP1E],[-qP2E -qP1E],'ko') ;
hold on
plot([pP2o pP1o],[qP2o qP1o],'g*') ;
hold on
plot(xfitP2,qC2fit(xfitP2),'r') ;
hold on
plot(xfitP2,qE2fit(xfitP2),'--r') ;
hold on
plot(xfitP1,qC1fit(xfitP1),'b') ;
hold on
plot(xfitP1,qE1fit(xfitP1),'--b') ;
hold on
grid on

xlim([-30 400]);
ylim([-300 400]);
Xlab = xlabel('p [MPa]','Interpreter','latex','FontSize',12);
Ylab = ylabel('q [MPa]','Interpreter','latex','FontSize',12);
tit = title(name,'Interpreter','latex','FontSize',14);
set(tit,'FontSize',14);
set(Ylab,'FontSize',12);
set(Xlab,'FontSize',12);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% sig2-sig1 plot
figure(3)
x = linspace(0,max(sig2)+100) ;
y = linspace(0,max(sig1)+100) ;
plot(sig2,sig1,'ko') ;
hold on
[e f] = size(sig3) ;
for i = 1:e
    hold on
    if sig3(i) == sigconf(1) 
        plot(sig2(i),sig1(i),'ko','MarkerFaceColor','white') 
        hold on
    elseif sig3(i) == sigconf(2)
        plot(sig2(i),sig1(i),'ko','MarkerFaceColor','black') 
        hold on
    elseif sig3(i) == sigconf(3)
        plot(sig2(i),sig1(i),'ks','MarkerFaceColor','white') 
        hold on
    elseif sig3(i) == sigconf(4)
        plot(sig2(i),sig1(i),'ks','MarkerFaceColor','black') 
        hold on
    elseif sig3(i) == sigconf(5)
        plot(sig2(i),sig1(i),'k^','MarkerFaceColor','white') 
        hold on
    elseif sig3(i) == sigconf(6)
        plot(sig2(i),sig1(i),'k^','MarkerFaceColor','black') 
        hold on
    elseif sig3(i) == sigconf(7)
        plot(sig2(i),sig1(i),'kd','MarkerFaceColor','white') 
        hold on
    end
    hold on
end

ye = @(x) x ;
yt = @(x) 2*x ;
yc2 = @(x) (VoP2/NcP2)+((NcP2+1)/NcP2)*x ;
yc1 = @(x) (VoP1/NcP1)+((NcP1+1)/NcP1)*x ;
% % yt = @(x) (1/(NcP2+NeP2-1))*(VoP2+(NcP2-NeP2)*x) ;
ii = plot(x,sig1P2(x,sigconf(2)),'r') ;
hold on
jj = plot(x,sig1P1(x,sigconf(2)),'b') ;
hold on
kk = plot(x,ye(x),'--k') ;
hold on
ll = plot(x,yc2(x),'k:') ;
hold on
nn = plot(x,yt(x),'k:') ;
hold on
% plot(x,yc1(x),'--k') ;
% hold on
% % plot(x,yt(x),'--g') ;
% % hold on
%
%'P2 plane','P1 plane',
xlim([0 max(sig2)+100]);
ylim([0 max(sig1)+100]);
legend([kk ll],'\sigma_{II}=\sigma_{I} line (extension)','\sigma_{II}=\sigma_{III} line (compression)')
Xlab = xlabel('$$\sigma_{II}$$ [MPa]','Interpreter','latex','FontSize',12);
Ylab = ylabel('$$\sigma_I$$ [MPa]','Interpreter','latex','FontSize',12);
tit = title(name,'Interpreter','latex','FontSize',14);
set(tit,'FontSize',14);
set(Ylab,'FontSize',12);
set(Xlab,'FontSize',12);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


figure(8)
% 3D plot
plot3(sig2,sig3,sig1,'ko','MarkerFaceColor','black') ;
hold on
% x1 = linspace(-30,200) ;
% x2 = linspace(-30,200) ;
% x3 = linspace(-30,200) ;

% % point 1111
% syms x y z
% eqn1 = (z*AP1l+x*BP1l-y*CP1l)==1 ;
% eqn2 = (z*AP2l+x*BP2l-y*CP2l)==1 ;
% eqn3 = z*AP1l+x*BP1l-y*CP1l==z*AP2l+x*BP2l-y*CP2l ;
% eqn4 = (z*AP1l+y*BP1l-x*CP1l)==1 ;
% eqn5 = (z*AP2l+y*BP2l-x*CP2l)==1 ;
% eqn6 = z*AP1l+y*BP1l-x*CP1l==z*AP2l+y*BP2l-x*CP2l ;
% S1111 = solve([eqn1 eqn2 eqn3 eqn4 eqn5 eqn6] ,[x y z]) ;
% x1111 = double(S1111.x) ;
% y1111 = double(S1111.y) ;
% z1111 = double(S1111.z) ;
% plot3(x1111,y1111,z1111,'s','MarkerFaceColor','green')
% hold on
% % point 1122x
% eqn1 = (z*AP1l+x*BP1l-y*CP1l)==1 ;
% eqn2 = (z*AP2l+x*BP2l-y*CP2l)==1 ;
% eqn3 = (z*BP2l+x*AP2l-y*CP2l)==1 ;
% eqn4 = z*AP1l+x*BP1l-y*CP1l==z*AP2l+x*BP2l-y*CP2l ;
% eqn5 = z*AP1l+x*BP1l-y*CP1l==z*BP2l+x*AP2l-y*CP2l ;
% S112x = solve([eqn1 eqn2 eqn3 eqn4 eqn5] ,[x y z]) ;
% x112x = double(S112x.x) ;
% y112x = double(S112x.y) ;
% z112x = double(S112x.z) ;
% plot3(x112x,y112x,z112x,'s','MarkerFaceColor','green')
% hold on
% % point 1122y
% eqn1 = (z*AP1l+y*BP1l-x*CP1l)==1 ;
% eqn2 = (z*AP2l+y*BP2l-x*CP2l)==1 ;
% eqn3 = (z*BP2l+y*AP2l-x*CP2l)==1 ;
% eqn4 = z*AP1l+y*BP1l-x*CP1l==z*AP2l+y*BP2l-x*CP2l ;
% eqn5 = z*AP1l+y*BP1l-x*CP1l==z*BP2l+y*AP2l-x*CP2l ;
% S112y = solve([eqn1 eqn2 eqn3 eqn4 eqn5] ,[x y z]) ;
% x112y = double(S112y.x) ;
% y112y = double(S112y.y) ;
% z112y = double(S112y.z) ;
% plot3(x112y,y112y,z112y,'s','MarkerFaceColor','green')
% hold on
% % point 3333
% eqn1 = (-z*CP2l+x*AP2l+y*BP2l)==1 ;
% eqn2 = (-z*CP1l+x*AP1l+y*BP1l)==1 ;
% eqn3 = (-z*CP1l+y*AP1l+x*BP1l)==1 ;
% eqn4 = (-z*CP2l+y*AP2l+x*BP2l)==1 ;
% eqn5 = -z*CP2l+x*AP2l+y*BP2l==-z*CP1l+x*AP1l+y*BP1l;
% eqn6 = -z*CP1l+y*AP1l+x*BP1l==-z*CP2l+y*AP2l+x*BP2l ;
% S333 = solve([eqn1 eqn2 eqn3 eqn4 eqn5 eqn6] ,[x y z]) ;
% x333 = double(S333.x) ;
% y333 = double(S333.y) ;
% z333 = double(S333.z) ;
% plot3(x333,y333,z333,'s','MarkerFaceColor','green')
% hold on
% % point 3322x
% eqn1 = (-z*CP2l+x*AP2l+y*BP2l)==1 ;
% eqn2 = (-z*CP1l+x*AP1l+y*BP1l)==1 ;
% eqn3 = (z*BP2l+x*AP2l-y*CP2l)==1 ;
% eqn4 = -z*CP2l+x*AP2l+y*BP2l==-z*CP1l+x*AP1l+y*BP1l;
% eqn5 = -z*CP2l+x*AP2l+y*BP2l==z*BP2l+x*AP2l-y*CP2l;
% S332x = solve([eqn1 eqn2 eqn3 eqn4 eqn5] ,[x y z]) ;
% x332x = double(S332x.x) ;
% y332x = double(S332x.y) ;
% z332x = double(S332x.z) ;
% plot3(x332x,y332x,z332x,'s','MarkerFaceColor','green')
% hold on
% % point 3322y
% eqn1 = (-z*CP2l+y*AP2l+x*BP2l)==1 ;
% eqn2 = (-z*CP1l+y*AP1l+x*BP1l)==1 ;
% eqn3 = (z*BP2l+y*AP2l-x*CP2l)==1 ;
% eqn4 = -z*CP2l+y*AP2l+x*BP2l==-z*CP1l+y*AP1l+x*BP1l;
% eqn5 = -z*CP2l+y*AP2l+x*BP2l==z*BP2l+y*AP2l-x*CP2l;
% S332y = solve([eqn1 eqn2 eqn3 eqn4 eqn5] ,[x y z]) ;
% x332y = double(S332y.x) ;
% y332y = double(S332y.y) ;
% z332y = double(S332y.z) ;
% plot3(x332y,y332y,z332y,'s','MarkerFaceColor','green')
% hold on

[x1,y1] = meshgrid(-100:5:max([max(sig1) max(sig2) max(sig3)])+600);
[x2,y2] = meshgrid(-100:5:max([max(sig1) max(sig2) max(sig3)])+600);
surf2 = surf(x2,y2,zs1P2x(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf1 = surf(x1,y1,zs1P1x(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf3 = surf(x2,y2,zs3P2x(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf4 = surf(x1,y1,zs3P1x(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf5 = surf(x1,y1,zs2P1x(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf6 = surf(x2,y2,zs2P2x(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf7 = surf(x2,y2,zs1P2y(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf8 = surf(x1,y1,zs1P1y(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf9 = surf(x2,y2,zs3P2y(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf10 = surf(x1,y1,zs3P1y(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf11 = surf(x1,y1,zs2P1y(x1,y1),'FaceColor',[0 0 1],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on
surf12 = surf(x2,y2,zs2P2y(x2,y2),'FaceColor',[1 0 0],'FaceAlpha',0.4,'EdgeColor','none') ;
hold on

% 
% zd112 = zs1P2x(x2,y2) - zs1P2y(x2,y2);
% zd111 = zs1P1x(x1,y1) - zs1P1y(x1,y1);
% zd122x = zs1P2x(x2,y2) - zs2P2x(x2,y2);
% zd121x = zs1P1x(x1,y1) - zs2P1x(x1,y1);
% zd122y = zs1P2y(x2,y2) - zs2P2y(x2,y2);
% zd121y = zs1P1y(x1,y1) - zs2P1y(x1,y1);
% zd332 = zs3P2x(x2,y2) - zs3P2y(x2,y2);
% zd331 = zs3P1x(x1,y1) - zs3P1y(x1,y1);
% zd322x = zs3P2x(x2,y2) - zs2P2x(x2,y2);
% zd321x = zs3P1x(x1,y1) - zs2P1x(x1,y1);
% zd322y = zs3P2y(x2,y2) - zs2P2y(x2,y2);
% zd321y = zs3P1y(x1,y1) - zs2P1y(x1,y1);
% 
% zd11x = zs1P2x(x2,y2) - zs1P1x(x1,y1);
% zd22x = zs2P2x(x2,y2) - zs2P1x(x1,y1);
% zd33x = zs3P2x(x2,y2) - zs3P1x(x1,y1);
% zd11y = zs1P2y(x2,y2) - zs1P1y(x1,y1);
% zd22y = zs2P2y(x2,y2) - zs2P1y(x1,y1);
% zd33y = zs3P2y(x2,y2) - zs3P1y(x1,y1);
% %-------------
% [NxL112,NyL112,NzL112,Bx112,By112,Bz112] = contP2(x2,y2,zs1P2x(x2,y2),zd112,x1111,y1111,z1111) ;
% NL112 = [NxL112;NyL112;NzL112] ;
% B112 = [Bx112,By112,Bz112];
% l112 = line(NxL112,NyL112,NzL112, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL111,NyL111,NzL111,Bx111,By111,Bz111] = contP1(x1,y1,zs1P1x(x1,y1),zd111,x1111,y1111,z1111) ;
% NL111 = [NxL111;NyL111;NzL111] ;
% B111 = [Bx111,By111,Bz111];
% l111 = line(NxL111, NyL111, NzL111, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL122x,NyL122x,NzL122x,Bx122x,By122x,Bz122x] = contP2(x2,y2,zs1P2x(x2,y2),zd122x,x112x,y112x,z112x) ;
% NL122x = [NxL122x;NyL122x;NzL122x] ;
% B122x = [Bx122x,By122x,Bz122x];
% l122x = line(NxL122x, NyL122x, NzL122x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL122y,NyL122y,NzL122y,Bx122y,By122y,Bz122y] = contP2(x2,y2,zs1P2y(x2,y2),zd122y,x112y,y112y,z112y) ;
% NL122y = [NxL122y;NyL122y;NzL122y] ;
% B122y = [Bx122y,By122y,Bz122y];
% l122y = line(NxL122y, NyL122y, NzL122y, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL121x,NyL121x,NzL121x,Bx121x,By121x,Bz121x] = contP1(x1,y1,zs1P1x(x1,y1),zd121x,x112x,y112x,z112x) ;
% NL121x = [NxL121x;NyL121x;NzL121x] ;
% B121x = [Bx121x,By121x,Bz121x];
% l121x = line(NxL121x, NyL121x, NzL121x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL121y,NyL121y,NzL121y,Bx121y,By121y,Bz121y] = contP1(x1,y1,zs1P1y(x1,y1),zd121y,x112y,y112y,z112y) ;
% NL121y = [NxL121y;NyL121y;NzL121y] ;
% B121y = [Bx121y,By121y,Bz121y];
% l121y = line(NxL121y, NyL121y, NzL121y, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL332,NyL332,NzL332,Bx332,By332,Bz332] = contP2(x2,y2,zs3P2x(x2,y2),zd332,x333,y333,z333) ;
% NL332 = [NxL332;NyL332;NzL332] ;
% B332 = [Bx332,By332,Bz332];
% l332 = line(NxL332, NyL332, NzL332, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL331,NyL331,NzL331,Bx331,By331,Bz331] = contP1(x1,y1,zs3P1x(x1,y1),zd331,x333,y333,z333) ;
% NL331 = [NxL331;NyL331;NzL331] ;
% B331 = [Bx331,By331,Bz331];
% l331 = line(NxL331, NyL331, NzL331, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL322x,NyL322x,NzL322x,Bx322x,By322x,Bz322x] = contP2(x2,y2,zs3P2x(x2,y2),zd322x,x332x,y332x,z332x) ;
% NL322x = [NxL322x;NyL322x;NzL322x] ;
% B322x = [Bx322x,By322x,Bz322x];
% l322x = line(NxL322x, NyL322x, NzL322x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL322y,NyL322y,NzL322y,Bx322y,By322y,Bz322y] = contP2(x2,y2,zs3P2y(x2,y2),zd322y,x332y,y332y,z332y) ;
% NL322y = [NxL322y;NyL322y;NzL322y];
% B322y = [Bx322y,By322y,Bz322y];
% l322y = line(NxL322y, NyL322y, NzL322y, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL321x,NyL321x,NzL321x,Bx321x,By321x,Bz321x] = contP1(x1,y1,zs3P1x(x1,y1),zd321x,x332x,y332x,z332x) ;
% NL321x = [NxL321x;NyL321x;NzL321x] ;
% B321x = [Bx321x,By321x,Bz321x];
% l321x = line(NxL321x, NyL321x, NzL321x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL321y,NyL321y,NzL321y,Bx321y,By321y,Bz321y] = contP1(x1,y1,zs3P1y(x1,y1),zd321y,x332y,y332y,z332y) ;
% NL321y = [NxL321y;NyL321y;NzL321y] ;
% B321y = [Bx321y,By321y,Bz321y];
% l321y = line(NxL321y, NyL321y, NzL321y, 'Color', 'k', 'LineWidth', 1);
% hold on
% %------%%%%--------%%%---------%%%-----------------------
% [NxL11x, NyL11x, NzL11x,Bx11x,By11x,Bz11x] = contPint(x2,y2,zs1P2x(x2,y2),zd11x,x1111,y1111,z1111,x112x,y112x,z112x) ;
% NL11x = [NxL11x; NyL11x; NzL11x] ;
% B11x = [Bx11x,By11x,Bz11x];
% l11x = line(NxL11x, NyL11x, NzL11x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL11y, NyL11y, NzL11y,Bx11y,By11y,Bz11y] = contPint(x2,y2,zs1P2y(x2,y2),zd11y,x1111,y1111,z1111,x112y,y112y,z112y) ;
% NL11y = [NxL11y; NyL11y; NzL11y] ;
% B11y = [Bx11y,By11y,Bz11y];
% l11y = line(NxL11y, NyL11y, NzL11y, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL22x, NyL22x, NzL22x,Bx22x,By22x,Bz22x] = contPint(x2,y2,zs2P2x(x2,y2),zd22x,x332x,y332x,z332x,x112x,y112x,z112x) ;
% NL22x = [NxL22x; NyL22x; NzL22x] ;
% B22x = [Bx22x,By22x,Bz22x];
% l22x = line(NxL22x, NyL22x, NzL22x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL22y, NyL22y, NzL22y,Bx22y,By22y,Bz22y] = contPint(x2,y2,zs2P2y(x2,y2),zd22y,x332y,y332y,z332y,x112y,y112y,z112y) ;
% NL22y = [NxL22y; NyL22y; NzL22y] ;
% B22y = [Bx22y,By22y,Bz22y];
% l22y = line(NxL22y, NyL22y, NzL22y, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL33x, NyL33x, NzL33x,Bx33x,By33x,Bz33x] = contPint(x2,y2,zs3P2x(x2,y2),zd33x,x332x,y332x,z332x,x333,y333,z333) ;
% NL33x = [NxL33x; NyL33x; NzL33x] ;
% B33x = [Bx33x,By33x,Bz33x];
% l33x = line(NxL33x, NyL33x, NzL33x, 'Color', 'k', 'LineWidth', 1);
% hold on
% [NxL33y, NyL33y, NzL33y,Bx33y,By33y,Bz33y] = contPint(x2,y2,zs3P2y(x2,y2),zd33y,x332y,y332y,z332y,x333,y333,z333) ;
% NL33y = [NxL33y; NyL33y; NzL33y] ;
% B33y = [Bx33y,By33y,Bz33y];
% l33y = line(NxL33y, NyL33y, NzL33y, 'Color', 'k', 'LineWidth', 1);
% hold on
% %-------------------------------------------------------------------
% % 
% % [X1,Y1,Z1] = ptFill(NL122x,NL112,NL11x,B122x,B112,B11x) ;
% % fill3(X1,Y1,Z1,'r','FaceAlpha',0.4);
% % hold on
% % [X2,Y2,Z2] = ptFill(NL122y,NL112,NL11y,B122y,B112,B11y) ;
% % fill3(X2,Y2,Z2,'r','FaceAlpha',0.4);
% % hold on
% % [X3,Y3,Z3] = ptFill(NL121x,NL111,NL11x,B121x,B111,B11x) ;
% % fill3(X3,Y3,Z3,'b','FaceAlpha',0.4);
% % hold on
% % [X4,Y4,Z4] = ptFill(NL121y,NL111,NL11y,B121y,B111,B11y) ;
% % fill3(X4,Y4,Z4,'b','FaceAlpha',0.4);
% % hold on
% % %-----------------
% % [X5,Y5,Z5] = ptFill(NL122x,NL322x,NL22x,B122x,B322x,B22x) ;
% % fill3(X5,Y5,Z5,'r','FaceAlpha',0.4);
% % hold on
% % [X6,Y6,Z6] = ptFill(NL122y,NL322y,NL22y,B122y,B322y,B22y) ;
% % fill3(X6,Y6,Z6,'r','FaceAlpha',0.4);
% % hold on
% % [X7,Y7,Z7] = ptFill(NL121x,NL321x,NL22x,B121x,B321x,B22x) ;
% % fill3(X7,Y7,Z7,'b','FaceAlpha',0.4);
% % hold on
% % [X8,Y8,Z8] = ptFill(NL121y,NL321y,NL22y,B121y,B321y,B22y) ;
% % fill3(X8,Y8,Z8,'b','FaceAlpha',0.4);
% % hold on
% % %----------------------
% % [X9,Y9,Z9] = ptFill(NL322x,NL332,NL33x,B322x,B332,B33x) ;
% % fill3(X9,Y9,Z9,'r','FaceAlpha',0.4);
% % hold on
% % [X10,Y10,Z10] = ptFill(NL322y,NL332,NL33y,B322y,B332,B33y) ;
% % fill3(X10,Y10,Z10,'r','FaceAlpha',0.4);
% % hold on
% % [X11,Y11,Z11] = ptFill(NL321x,NL331,NL33x,B321x,B331,B33x) ;
% % fill3(X11,Y11,Z11,'b','FaceAlpha',0.4);
% % hold on
% % [X12,Y12,Z12] = ptFill(NL321y,NL331,NL33y,B321y,B331,B33y) ;
% % fill3(X12,Y12,Z12,'b','FaceAlpha',0.4);
% % hold on

%grid on
xlim([-100 max([max(sig1) max(sig2) max(sig3)])+500]);
ylim([-100 max([max(sig1) max(sig2) max(sig3)])+500]);
zlim([-100 max([max(sig1) max(sig2) max(sig3)])+500]) ;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%ax.ZAxisLocation = 'origin';
Xlab = xlabel('$$\sigma_2$$ [MPa]','Interpreter','latex','FontSize',12);
Ylab = ylabel('$$\sigma_3$$ [MPa]','Interpreter','latex','FontSize',12);
Zlab = zlabel('$$\sigma_1$$ [MPa]','Interpreter','latex','FontSize',12);
tit = title(name,'Interpreter','latex','FontSize',14);
%axis off