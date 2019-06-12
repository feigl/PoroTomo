function S = read_and_mesh_faults(geologic_model,make_plots)
%% read faults from Siler model and mesh


% 20180209 Kurt Feigl

if contains(geologic_model,'Jolie') == 1
    %% Jolie et al. Models
    %dirname = '/Users/feigl/BoxSync/PoroTomo/BradysEGSProject/DATA v2012 02 14/Fault_Files/';
    dirname = '/Users/feigl/BoxSync/PoroTomo/BradysEGSProject/DATAv20120214/Fault_Files/';
    
    %filname = '1001ft_1';
    %nfaults = 1
    % All the fault files start with the number 1
    flisting=flist(strcat(dirname,'1*'))
    [nr,nc] = size(flisting);
    nfaults = nr
    %nfaults = 5 

    %% Select most important faults
    % From: Nicholas Davatzes <davatzes@temple.edu>
    % Date: 2017- Jun-29 Thursday at 15:23
    % To: Kurt Feigl <feigl@wisc.edu>
    % Subject: fault file ranking
    %
    % Kurt,
    %
    % You should have received a Google Drive invitation so that you can
    % download notes to help you with picking faults; this of course
    % neglects the import geological horizons which could play a role as
    % aquifers or whose boundaries along faults might mark interesting
    % transitions.
    %
    % Cheers,
    % Nick
    % please use this share to download the folder with notes on relative
    % importance of the faults to your hard drive. There is a PPT with some
    % visualization and labeling for reference, as well as tables of the 4-5
    % most important faults. Also there is an annotated xlsx spreadsheet
    % identifying the most important faults.
    % Top 5			Notes
    % fblk	1027	ft_1	important crossing fault, Conjugate to fblk1011ft_1 and fblk1024_1
    % fblk	1012	ft_1
    % fblk	1012	fta_1
    % fblk	1011	ft_1
    % fblk	1015	_2
    % fblk	1024	_1	Northern termination abuts fblk1011ft_1_remesh.ts
    % fblk	1007	ft_1	or use 1009
    % Others to consider
    % fblk1028_1			important crossing fault
    % fblk1010_3
    % To explore, use VTK files and Paraview or Matlab Scripts
    
    if strcmp(geologic_model,'Jolie9')
        kount=0;
        kount=kount+1; flisting{kount} = strcat(dirname,'1027ft_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1012ft_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1012fta_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1011ft_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1015_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1024_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1007ft_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1028_1');
        kount=kount+1; flisting{kount} = strcat(dirname,'1010_1');
        flisting = flisting'; % make into a column vector of strings
        nfaults = kount;
    end
    
    
    for i=1:nfaults
    filename = flisting{i}
        f1 = load(filename);
        [nr,nc] = size(f1);
        npoints = nr;
        
        
        % take only a few of the faults
        %     ikeep = randi(npoints,max([ceil(npoints/10),100]),1);
        %     ikeep = unique(ikeep);
        % keep all the points
        ikeep=1:npoints;
        
        nkeep = numel(ikeep);
        
        
        
        % Read X,Y,Z file
        %     fprintf(1,'Plotting fault number %d %s\n',i,char(F.Fault_names(i)));
        %     S(i).name = char(F.Fault_names(i));
        S(i).filename = filename;
        S(i).faultname = strrep(filename,dirname,'')
        %     f1=cell2mat(F.Faults_xyz(i));
        S(i).UTM.e  = colvec(f1(ikeep,1));
        S(i).UTM.n  = colvec(f1(ikeep,2));
        S(i).UTM.v  = colvec(f1(ikeep,3));
        %     npoints = numel(S(i).UTM.e);
    end
elseif strcmp(geologic_model,'Siler') == 1
    
    %% Siler Model
    % From: Kurt Feigl Sent: Thursday, February 15, 2018 12:03 AM To: Siler,
    % Drew Subject: [porotomo] Brady Earthvision files
    %
    % Hi Drew,
    %
    % Thanks for sharing your files, as listed below. I confirm that we will
    % not share them outside the PoroTomo group and that we will cite the GRC
    % paper by Siler and Faulds (2013).
    %
    % Best regards,
    %
    % Kurt
    %
    % /Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler
    % -rwxrwxrwx  1 feigl  15   7256872 Feb 13 11:58 SilerBradyModel_strat.dat*
    % -rwxrwxrwx@ 1 feigl  15   4852180 Feb 13 11:57 SilerBradyModel_faults.dat*
    % -rwxrwxrwx@ 1 feigl  15   1778624 Jan  1  2014 DrewBradysModel.fault.unsliced.faces
    % -rwxrwxrwx  1 feigl  15  13263112 Sep  6  2013 DrewBradysModel.unsliced.faces*
    %
    %
    % Siler, D. L., and J. E. Faulds (2013), Three-Dimensional Geothermal
    % Fairway Mapping: examples from the Western Great Basin, USA, Geotherm Res
    % Counc Trans, 37.
    
    dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler'
    filename = strcat(dirname,filesep,'SilerBradyModel_faults.dat.xlsx')
    T=readtable(filename,'ReadVariableNames',true);
    SS=table2struct(T,'ToScalar',true)
    faultnames=unique(SS.FaultID);
    nfaults = numel(faultnames)
   %nfaults = 1 % for debugging
    
    for i = 1:nfaults
       ifault = find(strcmp(SS.FaultID,faultnames{i}) == 1);
       npoints = numel(ifault);
       fprintf(1,'i= %5d %20s %5d\n',i,faultnames{i},npoints);
       S(i).UTM.e = colvec(SS.UTMeasting_m(ifault));
       S(i).UTM.n = colvec(SS.UTMnorthing_m(ifault));
       S(i).UTM.v = colvec(SS.UTMelevation_m(ifault));
       S(i).filename = filename;
       S(i).faultname = faultnames{i};
    end
else
    error(sprintf('Unknown geologic_model %s\n',geologic_model));
end


%% fit surfaces to faults  - slow
ngood = 0;
for i=1:nfaults
    %convert from UTM to rotated coordinate system PR
    [S(i).PR.x, S(i).PR.y, S(i).PR.z] = utm2xyz_porotomo(S(i).UTM.e, S(i).UTM.n, S(i).UTM.v);
    
    titlestr = sprintf('%sPR.z_Min%06.0fm_Max%06.0fm',S(i).faultname,nanmin(S(i).PR.z),nanmax(S(i).PR.z))
    %     XYZ = [S(i).PR.x, S(i).PR.y, S(i).PR.z];
    %[S(i).kfaces, S(i).kbounds] = triangulatexyz(S(i).PR.x, S(i).PR.y, S(i).PR.z,titlestr);
    %[S(i).kfaces, S(i).kbounds] = triangulatexyz2(S(i).PR.x, S(i).PR.y, S(i).PR.z,titlestr);
    [S(i).TRI,S(i).trigood,S(i).Xtri,S(i).Ytri,S(i).Ztri] = triangulatexyz3(S(i).PR.x, S(i).PR.y, S(i).PR.z,titlestr,make_plots);
    
    %     %if isfinitearray(S(i).kbounds) == 1 && isfinitearray(S(i).kfaces) == 1
    %     if isfinitearray(S(i).kfaces) == 1
    if numel(S(i).trigood(S(i).trigood)) > 0
       ngood = ngood+1;
    end
end

return
end


