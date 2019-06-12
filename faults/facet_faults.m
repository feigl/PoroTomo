function FAULTS = facet_faults(FAULTS)
% get indices of vertices and faces
nfaults = numel(FAULTS);

for j=1:nfaults
    %for j=[15,59]
    if sum(FAULTS(j).trigood) ~= 0
        [ntriangles,n3] = size(FAULTS(j).TRI.ConnectivityList);
        %ntriangles = 10;
        fprintf(1,'Working on fault number %d named %s\n',j,char(FAULTS(j).faultname));
        [ntriangles,n3] = size(FAULTS(j).TRI.ConnectivityList);
        X2=zeros(3*ntriangles,1);
        Y2=zeros(3*ntriangles,1);
        Z2=zeros(3*ntriangles,1);
        FF=zeros(ntriangles,3);
        kk=0;
        k1=1;
        k3=3;
        for k=1:ntriangles
            % get indices of vertices
            kA=colvec(FAULTS(j).TRI.ConnectivityList(k,1));
            kB=colvec(FAULTS(j).TRI.ConnectivityList(k,2));
            kC=colvec(FAULTS(j).TRI.ConnectivityList(k,3));
            
            % get coordinates of vertices
            x2=colvec([FAULTS(j).Xtri(kA),FAULTS(j).Xtri(kB),FAULTS(j).Xtri(kC)]);
            y2=colvec([FAULTS(j).Ytri(kA),FAULTS(j).Ytri(kB),FAULTS(j).Ytri(kC)]);
            z2=colvec([FAULTS(j).Ztri(kA),FAULTS(j).Ztri(kB),FAULTS(j).Ztri(kC)]);
            
            
            X2(k1:k3,1)=x2;
            Y2(k1:k3,1)=y2;
            Z2(k1:k3,1)=z2;
            kk = kk+1;
            FF(kk,1:3)=rowvec([k1:k3]);
            k1 = k1+3;
            k3 = k3+3;
        end
        FAULTS(j).vertices = [X2,Y2,Z2];
        FAULTS(j).faces = FF;
    end
end
return
end

