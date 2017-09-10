function H = imageeq(X,Y,Z,titlestr,xlab,ylab,zlab,ctab,zlimits)
%function H = imageeq(X,Y,Z,titlestr,xlab,ylab,zlab,ctab,zlimits)
% given an array Z, plot it as an image, after equalizing the histogram
% Kurt Feigl with ideas from Xiangfang Zeng
% 20160921
figure;
hold on;
IM1=mat2gray(Z);

%% find extrema
xmin=min(X)
xmax=max(X)
ymin=min(Y)
ymax=max(Y)

switch zlimits
    case 'sig'       
        % make sure extrema are symmetric about zero
        zsig=std(reshape(Z,numel(Z),1));
        zmin=-0.1*zsig;
        zmax=+0.1*zsig;
    case 'abs'
        % make sure extrema are symmetric about zero
        zabs=max(max(abs(Z)));
        zmin=-1*zabs;
        zmax=+1*zabs;        
    case 'minmax'     
        % find actual extrema
        zmin=min(min(Z));
        zmax=max(max(Z));     
    otherwise
        % standard
        zmin=-1;
        zmax=1;   
end



% define coordinates
RI = imref2d(size(IM1));
RI.XWorldLimits = [xmin xmax];
RI.YWorldLimits = [ymin ymax];
xlabel(xlab,'Interpreter','none');
ylabel(ylab,'Interpreter','none');

imshow(histeq(IM1),RI);
axis image;
colormap(ctab);

title(titlestr,'Interpreter','none');
hc=colorbar;
hc.Label.String = zlab;
hc.Limits=[zmin zmax];

H=gcf;
return;



